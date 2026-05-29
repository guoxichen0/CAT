# Small runnable example for the CAT algorithm.

project_root <- "D:/Testability_of_IV_Set/CAT_Condition_Nonlinear/github_code"

source(file.path(project_root, "R", "data_generation.R"))
source(file.path(project_root, "R", "cat_main.R"))

set.seed(123)

data_obj <- gen_data_case1_no_W(
  function_set = default_function_set(),
  epsilon_set = default_epsilon_set(),
  sample_size = 1000
)

dat <- data_obj$data
iv_set <- dat[, grep("^IV", names(dat)), drop = FALSE]

result <- CAT(
  Treatment = dat$Treatment,
  Outcome = dat$Outcome,
  IV_set = iv_set,
  K = 2,
  effect_type = "nonconstant"
)

print(result)
