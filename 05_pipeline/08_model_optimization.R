# ML final optimization orchestrator.
# Work is split into development screening, locked external validation, and reporting.
source(project_path("scripts/08a_development_screen.R"), local = TRUE)
source(project_path("scripts/08b_locked_external_validation.R"), local = TRUE)
source(project_path("scripts/08c_result_summary.R"), local = TRUE)
