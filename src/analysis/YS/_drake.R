# Main file for drake to configure the analysis pipeline.
# For use with drake::r_make
#
# Author: Yuriy Sverchkov

# TODO: Check for installed packages and raise an informative message

src_folder <- here::here("src", "analysis", "YS")
data_folder <- here::here("data")
results_folder <- here::here("results")

source(file.path(src_folder, "functions.R"))
source(file.path(src_folder, "plan.R"))

drake::drake_config(plan)