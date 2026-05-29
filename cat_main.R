# Core CAT algorithm.

required_packages <- c("energy", "controlfunctionIV", "randomForest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install required packages before running CAT: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

rf_residual <- function(target, W_set) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Please install the randomForest package before using W_set.", call. = FALSE)
  }

  fit <- randomForest::randomForest(x = as.data.frame(W_set), y = as.numeric(target))
  as.numeric(target) - as.numeric(predict(fit, as.data.frame(W_set)))
}

linear_residual <- function(target, W_set) {
  df <- data.frame(target = as.numeric(target), as.data.frame(W_set))
  fit <- lm(target ~ ., data = df)
  as.numeric(residuals(fit))
}

get_residual <- function(target, W_set, residual_method) {
  if (residual_method == "randomForest") {
    return(rf_residual(target, W_set))
  }

  linear_residual(target, W_set)
}

residualize_by_W <- function(Treatment, Outcome, IV_set, W_set, residual_method) {
  W_set <- as.data.frame(W_set)
  IV_set <- as.data.frame(IV_set)
  residual_IV_set <- data.frame(
    lapply(IV_set, get_residual, W_set = W_set, residual_method = residual_method),
    check.names = FALSE
  )
  colnames(residual_IV_set) <- colnames(IV_set)

  list(
    Treatment = get_residual(Treatment, W_set, residual_method),
    Outcome = get_residual(Outcome, W_set, residual_method),
    IV_set = residual_IV_set
  )
}

CAT <- function(
  Treatment,
  Outcome,
  IV_set,
  K,
  W_set = NULL,
  effect_type = c("nonconstant", "constant"),
  residual_method = c("randomForest", "linear")
) {
  effect_type <- match.arg(effect_type)
  residual_method <- match.arg(residual_method)
  IV_set <- as.data.frame(IV_set)

  if (!is.null(W_set)) {
    residualized <- residualize_by_W(Treatment, Outcome, IV_set, W_set, residual_method)
    Treatment <- residualized$Treatment
    Outcome <- residualized$Outcome
    IV_set <- residualized$IV_set
  }

  df <- data.frame(Y = Outcome, D = Treatment)
  df <- cbind(df, IV_set)
  df <- df[complete.cases(df), ]
  df <- df[apply(df, 1, function(row) all(is.finite(row))), ]

  Y <- df$Y
  D <- df$D
  IV_set <- df[, colnames(IV_set), drop = FALSE]
  num_IVs <- ncol(IV_set)
  IV_names <- colnames(IV_set)

  if (K < 2 || K > num_IVs) {
    stop("K must be between 2 and the number of candidate IVs.")
  }

  V_list <- list()
  for (j in seq_len(num_IVs)) {
    Z <- IV_set[[IV_names[j]]]

    if (effect_type == "nonconstant") {
      model <- controlfunctionIV::cf(Y ~ D + I(D^2) + I(D^3) + I(D^4) | Z + I(Z^2) + I(Z^3))
      coefs <- coef(model)
      D_terms <- cbind(D, D^2, D^3, D^4)
      h_D <- D_terms %*% coefs[2:5]
    } else {
      denominator <- cov(D, Z)
      if (is.na(denominator) || abs(denominator) < .Machine$double.eps) {
        stop("Cannot estimate constant effect because cov(Treatment, IV) is zero for ", IV_names[j], ".")
      }
      beta <- cov(Y, Z) / denominator
      h_D <- beta * D
    }

    V_list[[IV_names[j]]] <- as.numeric(Y - h_D)
  }

  dcor_matrix <- matrix(0, nrow = num_IVs, ncol = num_IVs)
  all_pairs <- combn(num_IVs, 2)

  for (j in seq_len(ncol(all_pairs))) {
    idx1 <- all_pairs[1, j]
    idx2 <- all_pairs[2, j]

    Z1 <- IV_set[[IV_names[idx1]]]
    Z2 <- IV_set[[IV_names[idx2]]]
    V1 <- V_list[[IV_names[idx1]]]
    V2 <- V_list[[IV_names[idx2]]]

    pair_dcor_sum <- energy::dcor(V1, Z2) + energy::dcor(V2, Z1)
    dcor_matrix[idx1, idx2] <- pair_dcor_sum
    dcor_matrix[idx2, idx1] <- pair_dcor_sum
  }

  candidate_combinations <- combn(num_IVs, K)
  test_results <- numeric(ncol(candidate_combinations))

  for (p in seq_len(ncol(candidate_combinations))) {
    current_IV_indices <- candidate_combinations[, p]
    inner_pairs <- combn(length(current_IV_indices), 2)
    set_dcor_sum <- 0

    for (j in seq_len(ncol(inner_pairs))) {
      real_idx1 <- current_IV_indices[inner_pairs[1, j]]
      real_idx2 <- current_IV_indices[inner_pairs[2, j]]
      set_dcor_sum <- set_dcor_sum + dcor_matrix[real_idx1, real_idx2]
    }

    test_results[p] <- set_dcor_sum
  }

  min_abs_col_index <- which.min(test_results)
  best_dcor <- test_results[min_abs_col_index]
  best_combination_index <- candidate_combinations[, min_abs_col_index]
  best_IV_names <- IV_names[best_combination_index]

  list(
    best_IV_names = best_IV_names,
    best_dcor = best_dcor
  )
}
