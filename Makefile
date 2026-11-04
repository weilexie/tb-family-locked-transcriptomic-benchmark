RSCRIPT ?= Rscript

.PHONY: prepare_data optimize_final enhance_final stress_test manuscript_final

prepare_data:
	$(RSCRIPT) run_all.R --mode prepare_data

optimize_final:
	$(RSCRIPT) run_all.R --mode optimize_final

enhance_final:
	$(RSCRIPT) run_all.R --mode enhance_final

stress_test:
	$(RSCRIPT) run_all.R --mode stress_test

manuscript_final:
	$(RSCRIPT) run_all.R --mode manuscript_final
