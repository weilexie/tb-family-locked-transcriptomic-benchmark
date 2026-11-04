config <- read_config("config/config.yml")
core <- unique(unlist(config$packages$minimal_core, use.names = FALSE))
optional <- unique(unlist(config$packages[setdiff(names(config$packages), "minimal_core")], use.names = FALSE))
installed <- rownames(installed.packages())

availability <- data.frame(
  package = c(core, optional),
  required_level = c(rep("core", length(core)), rep("optional", length(optional))),
  available = c(core, optional) %in% installed,
  stringsAsFactors = FALSE
)
write_tsv_atomic(availability, "results/tables/package_availability.tsv")

versions <- data.frame(
  item = c("R", "platform", "libPaths"),
  value = c(as.character(getRversion()), R.version$platform, paste(.libPaths(), collapse = "; ")),
  stringsAsFactors = FALSE
)
for (pkg in availability$package[availability$available]) {
  versions <- rbind(versions, data.frame(item = pkg, value = as.character(utils::packageVersion(pkg)), stringsAsFactors = FALSE))
}
write_tsv_atomic(versions, "results/tables/software_versions.tsv")

missing_core <- availability$package[availability$required_level == "core" & !availability$available]
if (length(missing_core) > 0) {
  stop(sprintf("Missing core R packages: %s", paste(missing_core, collapse = ", ")), call. = FALSE)
}

write_lines(c(
  sprintf("[%s] 00_setup_packages completed", Sys.time()),
  sprintf("R version: %s", getRversion()),
  sprintf("Library paths: %s", paste(.libPaths(), collapse = "; "))
), "logs/00_setup_packages.log")

