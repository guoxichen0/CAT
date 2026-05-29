# Data generation utilities for nonlinear IV simulations.

function_lib <- function(U, index = c("non-linear", "random")) {
  U <- as.numeric(U)

  nonlinear_functions <- c(
    "log", "sin", "cos", "exp", "square", "cubic",
    "logistic_sigmoid", "sinh", "atanh", "2dpoly", "3dpoly",
    "log_quadratic", "exp_quadratic", "frequent_sin",
    "frequent_cos", "abs_sqrt", "step"
  )

  relationship_type <- index[1]
  nonlinear_type <- index[2]

  if (relationship_type == "linear") {
    gamma <- 0
    while (gamma == 0) {
      gamma <- runif(1, -1.5, 1.5)
    }
    return(gamma * U)
  }

  if (relationship_type != "non-linear") {
    stop("relationship_type must be either 'linear' or 'non-linear'.")
  }

  if (nonlinear_type == "random") {
    nonlinear_type <- sample(nonlinear_functions, 1)
  }

  out <- switch(
    nonlinear_type,
    log = log(abs(U) + 1),
    exp = exp(U),
    sin = sin(U),
    cos = cos(U),
    square = U^2,
    cubic = U^3,
    sinh = sinh(U),
    logistic_sigmoid = 1 / (1 + exp(-U)),
    atanh = atanh(pmin(pmax(U, -0.99), 0.99)),
    "2dpoly" = 1 - 1.5 * U + U^2,
    "3dpoly" = 1 - 1.5 * U + 0.1 * U^2 + 0.5 * U^3,
    log_quadratic = log(abs(0.5 * U^2 - 0.2 * U) + 1),
    exp_quadratic = exp(0.05 * U^2 - U + 1),
    frequent_sin = sin(3 * U),
    frequent_cos = cos(3 * U),
    abs_sqrt = sqrt(abs(U)),
    step = 1.0 * (U < 0) + 2.5 * (U >= 0),
    stop(paste("Unknown nonlinear type:", nonlinear_type))
  )

  max_abs <- max(abs(out), na.rm = TRUE)
  if (max_abs > 1e2) {
    out <- out / (10^(floor(log10(max_abs)) + 1))
  } else if (max_abs < 1e-2 && max_abs > 0) {
    out <- out * (10^(2 - floor(log10(max_abs))))
  }

  out[!is.finite(out)] <- 0
  as.numeric(out)
}

noise_lib <- function(sample_size, noise_type = "random") {
  noise_types <- c("Uniform", "Gaussian", "Beta", "Gamma", "Lognormal", "Exponential", "T")
  if (noise_type == "random") {
    noise_type <- sample(noise_types, 1)
  }

  out <- switch(
    noise_type,
    Uniform = runif(sample_size, -1, 1),
    Gaussian = rnorm(sample_size),
    Beta = rbeta(sample_size, 2, 2),
    Gamma = rgamma(sample_size, 2, 2),
    Lognormal = rlnorm(sample_size),
    Exponential = rexp(sample_size),
    T = rt(sample_size, 3),
    stop(paste("Unknown noise type:", noise_type))
  )

  as.numeric(out)
}

default_function_set <- function() {
  list(
    XtoY = c("non-linear", "sin"),
    UtoX = c("non-linear", "cubic"),
    UtoY = c("non-linear", "cubic"),
    W1toZ1 = c("linear", NA),
    W2toZ1 = c("linear", NA),
    W1toZ2 = c("linear", NA),
    W2toZ2 = c("linear", NA),
    W1toZ3 = c("linear", NA),
    W2toZ3 = c("linear", NA),
    W1toZ4 = c("linear", NA),
    W2toZ4 = c("linear", NA),
    W1toZ5 = c("linear", NA),
    W2toZ5 = c("linear", NA),
    W1toX = c("linear", NA),
    W2toX = c("linear", NA),
    W1toY = c("linear", NA),
    W2toY = c("linear", NA),
    Z1toX = c("non-linear", "cubic"),
    Z2toX = c("non-linear", "cubic"),
    Z3toX = c("non-linear", "cubic"),
    Z4toX = c("non-linear", "cos"),
    Z5toX = c("non-linear", "cos"),
    Z3toY = c("non-linear", "cubic"),
    Z4toZ5 = c("non-linear", "cubic"),
    Z4toY = c("non-linear", "cubic"),
    Z5toY = c("non-linear", "cubic")
  )
}

default_epsilon_set <- function() {
  list(
    U_epsilon = "Uniform",
    W_epsilon = "Uniform",
    Z_epsilon = "Uniform",
    X_epsilon = "Uniform",
    Y_epsilon = "Uniform"
  )
}

gen_data_case1_with_W <- function(function_set, epsilon_set, sample_size, noise_ratio = 0.5) {
  W1 <- noise_lib(sample_size, epsilon_set$W_epsilon)
  W2 <- noise_lib(sample_size, epsilon_set$W_epsilon)
  U <- noise_lib(sample_size, epsilon_set$U_epsilon)

  Z1 <- function_lib(W1, function_set$W1toZ1) +
    function_lib(W2, function_set$W2toZ1) +
    noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z2 <- function_lib(W1, function_set$W1toZ2) +
    function_lib(W2, function_set$W2toZ2) +
    noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z3 <- function_lib(W1, function_set$W1toZ3) +
    function_lib(W2, function_set$W2toZ3) +
    noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z4 <- function_lib(W1, function_set$W1toZ4) +
    function_lib(W2, function_set$W2toZ4) +
    noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z5 <- function_lib(W1, function_set$W1toZ5) +
    function_lib(W2, function_set$W2toZ5) +
    function_lib(Z4, function_set$Z4toZ5) +
    noise_lib(sample_size, epsilon_set$Z_epsilon)

  X <- function_lib(W1, function_set$W1toX) +
    function_lib(W2, function_set$W2toX) +
    function_lib(Z1, function_set$Z1toX) +
    function_lib(Z2, function_set$Z2toX) +
    function_lib(Z3, function_set$Z3toX) +
    function_lib(Z4, function_set$Z4toX) +
    function_lib(Z5, function_set$Z5toX) +
    (Z3 * Z4 * Z5)^3 +
    function_lib(U, function_set$UtoX) +
    noise_ratio * noise_lib(sample_size, epsilon_set$X_epsilon)

  Y <- function_lib(X, function_set$XtoY) +
    0.5 * function_lib(W1, function_set$W1toY) +
    0.5 * function_lib(W2, function_set$W2toY) +
    function_lib(U, function_set$UtoY) +
    function_lib(Z3, function_set$Z3toY) +
    1.5 * function_lib(Z4, function_set$Z4toY) +
    function_lib(Z5, function_set$Z5toY) +
    noise_ratio * noise_lib(sample_size, epsilon_set$Y_epsilon)

  data <- data.frame(
    IV1 = Z1,
    IV2 = Z2,
    IV3 = Z3,
    IV4 = Z4,
    IV5 = Z5,
    Cov1 = W1,
    Cov2 = W2,
    Treatment = X,
    Outcome = Y
  )

  data <- data[apply(data, 1, function(row) all(is.finite(row))), ]
  if (nrow(data) < sample_size) {
    stop("Generated data contains non-finite values. Change the distributions or function settings.")
  }
  rownames(data) <- NULL

  list(data = data)
}

gen_data_case1_no_W <- function(function_set, epsilon_set, sample_size, noise_ratio = 0.3) {
  U <- noise_lib(sample_size, epsilon_set$U_epsilon)
  Z1 <- noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z2 <- noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z3 <- noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z4 <- noise_lib(sample_size, epsilon_set$Z_epsilon)
  Z5 <- function_lib(Z4, function_set$Z4toZ5) +
    noise_ratio * noise_lib(sample_size, epsilon_set$Z_epsilon)

  X <- function_lib(Z1, function_set$Z1toX) +
    function_lib(Z2, function_set$Z2toX) +
    function_lib(Z3, function_set$Z3toX) +
    function_lib(Z4, function_set$Z4toX) +
    function_lib(Z5, function_set$Z5toX) +
    (Z3 * Z4 * Z5)^3 +
    function_lib(U, function_set$UtoX) +
    noise_ratio * noise_lib(sample_size, epsilon_set$X_epsilon)

  Y <- function_lib(X, function_set$XtoY) +
    function_lib(U, function_set$UtoY) +
    function_lib(Z3, function_set$Z3toY) +
    1.5 * function_lib(Z4, function_set$Z4toY) +
    function_lib(Z5, function_set$Z5toY) +
    noise_ratio * noise_lib(sample_size, epsilon_set$Y_epsilon)

  data <- data.frame(
    IV1 = Z1,
    IV2 = Z2,
    IV3 = Z3,
    IV4 = Z4,
    IV5 = Z5,
    Treatment = X,
    Outcome = Y
  )

  data <- data[apply(data, 1, function(row) all(is.finite(row))), ]
  if (nrow(data) < sample_size) {
    stop("Generated data contains non-finite values. Change the distributions or function settings.")
  }
  rownames(data) <- NULL

  list(data = data)
}


