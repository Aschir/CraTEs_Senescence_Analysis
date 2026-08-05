# Local backward-compatible override for craTEs::getActivities().
# Source after library(craTEs), then call getActivities() without namespace.

get_lambda_max <- function(Y, X, alpha) {
  mysd <- function(y) sqrt(sum((y - mean(y))^2) / length(y))

  stopifnot(length(Y) == dim(X)[1])

  (1 / (length(Y) * alpha)) *
    max(abs(apply(scale(X, center = TRUE, scale = apply(X, 2, mysd)) * Y, MARGIN = 2, sum)))
}

getActivities <- function(expr, X, alpha = 0.1, k = 5, nlambda = 100,
                          ncores = 1, log_id = "") {
  print("Going parallel !")

  n <- ncol(expr)

  cat("Computing lambda max over ", n, " samples\n")
  if (ncores == 1) {
    lambda_max_list <- vapply(seq_len(ncol(expr)), function(i) {
      get_lambda_max(expr[, i], X, alpha = alpha)
    }, numeric(1))
    names(lambda_max_list) <- colnames(expr)
  } else {
    cl <- parallel::makeCluster(
      ncores,
      type = "FORK",
      outfile = paste(c("./getActivities_log", log_id, ".txt"), collapse = "")
    )
    on.exit(parallel::stopCluster(cl), add = TRUE)
    lambda_max_list <- parallel::parApply(cl, expr, MARGIN = 2, FUN = get_lambda_max, X, alpha = alpha)
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
  }

  cat("Building the lambda grid")
  lambda_max <- max(lambda_max_list)
  log_lambda_grid <- seq(from = log(lambda_max), to = log(lambda_max * 1e-4), length.out = nlambda)
  lambda_values <- exp(log_lambda_grid)

  cat("Regressing activities of TFs for ", n, " samples\n")
  f <- function(y) {
    print("started new sample")

    glmnet::cv.glmnet(
      x = X,
      y = y,
      alpha = alpha,
      nfolds = k,
      lambda = lambda_values,
      type.measure = "mse"
    )
  }

  if (ncores == 1) {
    res <- lapply(seq_len(ncol(expr)), function(i) f(expr[, i]))
    names(res) <- colnames(expr)
  } else {
    cl <- parallel::makeCluster(
      ncores,
      type = "FORK",
      outfile = paste(c("./getActivities_log", log_id, ".txt"), collapse = "")
    )
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res <- parallel::parApply(cl, expr, MARGIN = 2, FUN = f)
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
  }

  cat("Finding the optimal lambda value across samples...\n")

  cat("\t... 1) Computing the per lambda MSE across samples\n")
  errors <- do.call("cbind", Biobase::subListExtract(res, "cvm"))

  lambda_mse <- rowSums(errors) / ncol(errors)
  lambda_mse_sd <- apply(errors, MARGIN = 1, sd)

  lambda_min_idx <- which.min(lambda_mse)
  lambda_min <- lambda_values[lambda_min_idx]

  cat("\t... 2) retrieving the coefficients, crossval MSE, r_squared for the optimal lambda value\n")
  coefs <- lapply(
    X = Biobase::subListExtract(res, "glmnet.fit"),
    FUN = glmnet::predict.glmnet,
    type = "coefficients",
    s = lambda_min
  )

  temp <- Biobase::subListExtract(res, "glmnet.fit")
  r_squared_temp <- Biobase::subListExtract(temp, "dev.ratio")
  r_squared <- lapply(r_squared_temp, FUN = function(x) x[lambda_min_idx])
  rm(list = "temp")

  r_squared_max <- lapply(r_squared_temp, FUN = max)
  best_lambda <- Biobase::subListExtract(res, "lambda.min")

  ans <- list()
  ans$r2 <- unlist(r_squared)
  ans$r2_max <- unlist(r_squared_max)
  ans$lambdas <- lambda_values
  ans$lambda_max_by_sample <- unlist(lambda_max_list)
  ans$lambda_min_overall <- lambda_min
  ans$lambdas_best <- unlist(best_lambda)
  ans$lambda_mse <- unlist(lambda_mse)
  ans$lambda_mse_sd <- unlist(lambda_mse_sd)

  ans$activities <- as.matrix(do.call(cbind, coefs))
  colnames(ans$activities) <- names(coefs)
  ans
}
