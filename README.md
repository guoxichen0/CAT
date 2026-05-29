# CAT for Nonlinear Instrumental Variable Set Testability

This repository contains R code for simulating nonlinear instrumental variable (IV) data and running the CAT algorithm to select candidate valid IV sets.

## Project structure

```text
.
|-- cat_main.R                 # Core CAT algorithm
|-- data_generation.R          # Simulated data generation functions
|-- run_small_example.R        # Minimal reproducible example
|-- README.md
`-- .gitignore
```

## Installation

Install the required R packages:

```r
install.packages(c("energy", "randomForest"))

# The randomForest package is used when covariates are adjusted by random forest residualization.

# Install controlfunctionIV if it is available from your package source.
# If it is installed from GitHub or a local source, add the exact command here.
# install.packages("controlfunctionIV")
```

## Quick start

Run the small example:

```r
setwd("path/to/CAT")
source("run_small_example.R")
```

The example will:

1. Generate a case1-sin nonlinear IV simulation dataset without covariates and sample size 1000.
2. Use all variables named `IV1`, `IV2`, ..., `IV5` as candidate instruments.
3. Run CAT with `K = 2`.
4. Print the selected IV set and the minimum distance-correlation score.

## Main components

### Data generation

`data_generation.R` contains:

- `function_lib()`: nonlinear and linear transformation library.
- `noise_lib()`: random noise generator.
- `default_function_set()`: default nonlinear relationships used in the example.
- `default_epsilon_set()`: default noise distributions.
- `gen_data_case1_no_W()`: case1-sin simulated data generator without observed covariates.
- `gen_data_case1_with_W()`: case1-sin simulated data generator with observed covariates.

### CAT algorithm

`cat_main.R` contains:

- `CAT()`: core algorithm for evaluating candidate IV subsets of size `K`.

The main inputs are:

```r
CAT(
  Treatment = data$Treatment,
  Outcome = data$Outcome,
  IV_set = data[, grep("^IV", names(data))],
  K = 2,
  effect_type = "nonconstant"
)
```

If `W_set` is provided, CAT first regresses Treatment, Outcome, and each IV on the covariates, then runs CAT on the residuals. If `W_set = NULL`, CAT runs directly on the original Treatment, Outcome, and IV set.

Two covariate-adjustment methods are available:

```r
# Random forest residualization
CAT(Treatment, Outcome, IV_set, K, W_set = W_set, residual_method = "randomForest")

# Linear regression residualization
CAT(Treatment, Outcome, IV_set, K, W_set = W_set, residual_method = "linear")
```

Use `effect_type = "nonconstant"` for the control-function IV estimator branch, or `effect_type = "constant"` for the constant-effect branch:

```r
beta <- cov(Outcome, IV) / cov(Treatment, IV)
V <- Outcome - beta * Treatment
```

The output is a list with the selected IV names and the corresponding CAT score:

```r
$best_IV_names
[1] "IV1" "IV2"

$best_dcor
[1] 0.07759382
```

## Citation

This repository is associated with:

- Testing the Validity of Instrumental Variable Sets in Causal Additive Models with Non-Constant Effects
- Data-Driven Selection of Instrumental Variables for Additive Nonlinear, Constant Effects Models

<!-- ## License

Add a license before making the repository public. Common options are MIT for permissive reuse or GPL-3 if you want derivative code to remain open source. -->
