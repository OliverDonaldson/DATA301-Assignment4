# Installs everything needed to render the assignment.
# Run once: Rscript install-deps.R
pkgs <- c("fpp2", "forecast", "ggplot2", "knitr", "rmarkdown")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
