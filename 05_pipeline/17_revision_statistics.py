"""Exact rank-pair summaries of the submitted, frozen external predictions.

The within/cross-family identity decomposes the AUROC statistic, not biological
and confounding effects. Conditional permutation expectations are descriptive;
no sample-level independence assumption or significance test is introduced.
"""
import argparse
import bisect
import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def concordance(positive, negative):
    negative = sorted(negative)
    lower = sum(bisect.bisect_left(negative, score) for score in positive)
    ties = sum(bisect.bisect_right(negative, score) -
               bisect.bisect_left(negative, score) for score in positive)
    return lower + 0.5 * ties, ties


def summarize(rows, selected):
    groups = defaultdict(list)
    for row in rows:
        score = float(row['score'])
        if not math.isfinite(score):
            raise ValueError('Non-finite selected-model prediction')
        groups[row['cohort_family']].append(
            (score, row['true_label'] == 'ATB', row['standard_sample_id']))
    sample_ids = [row['standard_sample_id'] for row in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError('Duplicate selected-model sample IDs')
    positives = [float(row['score']) for row in rows if row['true_label'] == 'ATB']
    negatives = [float(row['score']) for row in rows if row['true_label'] != 'ATB']
    total_pairs = len(positives) * len(negatives)
    total_credit, _ = concordance(positives, negatives)
    base = {key: selected[key] for key in ['task', 'selected_role', 'combination_id']}
    pair_rows, strata = [], defaultdict(lambda: [0, 0., 0])
    null_expected_credit = 0.
    for positive_family, first in sorted(groups.items()):
        pos = [score for score, label, _ in first if label]
        if not pos:
            continue
        for negative_family, second in sorted(groups.items()):
            neg = [score for score, label, _ in second if not label]
            if not neg:
                continue
            n_pairs = len(pos) * len(neg)
            credit, ties = concordance(pos, neg)
            stratum = 'within_family' if positive_family == negative_family else 'cross_family'
            strata[stratum][0] += n_pairs
            strata[stratum][1] += credit
            strata[stratum][2] += ties
            if positive_family == negative_family:
                expected_credit = n_pairs * 0.5
            else:
                all_credit, _ = concordance([r[0] for r in first], [r[0] for r in second])
                expected_credit = n_pairs * all_credit / (len(first) * len(second))
            null_expected_credit += expected_credit
            pair_rows.append(dict(base, positive_family=positive_family,
                negative_family=negative_family, pair_stratum=stratum,
                n_positive=len(pos), n_negative=len(neg), n_pairs=n_pairs,
                pair_weight=n_pairs/total_pairs, AUROC=credit/n_pairs,
                weighted_AUROC_term=credit/total_pairs, n_tied_pairs=ties,
                conditional_permutation_expected_AUROC=expected_credit/n_pairs))
    summary_rows = []
    for stratum in ['within_family', 'cross_family']:
        n_pairs, credit, ties = strata[stratum]
        summary_rows.append(dict(base, pair_stratum=stratum, n_samples=len(rows),
            n_positive=len(positives), n_negative=len(negatives), n_families=len(groups),
            n_pairs=n_pairs, pair_weight=n_pairs/total_pairs, AUROC=credit/n_pairs,
            weighted_AUROC_term=credit/total_pairs, n_tied_pairs=ties))
    reconstructed = sum(row['weighted_AUROC_term'] for row in summary_rows)
    pooled = total_credit/total_pairs
    if not math.isclose(reconstructed, pooled, abs_tol=1e-12):
        raise AssertionError('Pair decomposition does not reproduce pooled AUROC')
    expectation = null_expected_credit/total_pairs
    expectation_row = dict(base, n_samples=len(rows), n_families=len(groups),
        observed_AUROC=pooled, conditional_permutation_expected_AUROC=expectation,
        observed_minus_conditional_expectation=pooled-expectation,
        interpretation='fixed-family score-label randomization expectation; descriptive, not causal attribution',
        formal_p_value='not_computed')
    return summary_rows, pair_rows, expectation_row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    if root.name == 'TB ML':
        raise RuntimeError('Revision analysis must run in the revision copy')
    results = root/'07_results'/'03_ml_final'
    selected = read_tsv(results/'selected_combinations.tsv')
    predictions = read_tsv(results/'external_predictions.tsv.gz')
    summaries, pairs, expectations = [], [], []
    for model in selected:
        rows = [r for r in predictions if r['combination_id'] == model['combination_id']]
        # Some saved files repeat a combination for more than one exported role.
        role_rows = [r for r in rows if r.get('selected_role') == model['selected_role']]
        if role_rows:
            rows = role_rows
        one_summary, one_pairs, expectation = summarize(rows, model)
        summaries.extend(one_summary)
        pairs.extend(one_pairs)
        expectations.append(expectation)
    output = root/'07_results'/'04_revision'/'statistics'
    write_tsv(output/'external_AUROC_pair_decomposition.tsv', summaries)
    write_tsv(output/'external_AUROC_ordered_family_pairs.tsv', pairs)
    write_tsv(output/'conditional_permutation_expectation.tsv', expectations)
    for row in expectations:
        print(row['task'], row['selected_role'],
              f"AUROC={row['observed_AUROC']:.9f}",
              f"conditional_expectation={row['conditional_permutation_expected_AUROC']:.9f}")


if __name__ == '__main__':
    main()
