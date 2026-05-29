# Small runnable example for the CAT algorithm.

source(file.path(getwd(), "data_generation.R"))
source(file.path(getwd(), "cat_main.R"))

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
