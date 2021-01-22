# Master script for running everything.
# This version loads everything into the global environment.
#
# Author: Yuriy Sverchkov

# TODO: Check for installed packages and raise an informative message

src_folder <- here::here("src", "analysis", "YS")
data_folder <- here::here("data")
results_folder <- here::here("results")

source(file.path(src_folder, "functions.R"))
source(file.path(src_folder, "plan.R"))

# Run make
cmd_targets <- commandArgs(trailingOnly = TRUE)
if (length(cmd_targets) > 0) {
  drake::make(plan, targets = cmd_targets)
} else {
  drake::make(plan)
}