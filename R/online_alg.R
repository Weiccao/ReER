#' Calculate the diagnol element of expectile weight matrix
#'
#' @param tau  Expectile level.
#' @param X Covariates matrices.
#' @param Y Independent vector.
#' @param beta Parameters.
#'
#' @return Expectile weight vector.
weightMatrix <- function(tau, X, Y, beta) {
  epsilon <- Y - X %*% beta
  wDiag <- ifelse(epsilon < 0, 1 - tau, tau) %>% as.numeric()
  return(wDiag)
}

#' Calculate the asymptotic covariance matrix of expectile regression
#'
#' @param X Covariates matrices.
#' @param Y Independent vector.
#' @param beta Parameters.
#' @param tau  Expectile level.
#'
#' @return The asymptotic covariance matrix of expectile regression.
qMatrix <- function(X, Y, beta, tau) {
  uMatrix <- function(tau, X, Y, beta) {
    W_vec <- weightMatrix(tau, X, Y, beta)
    return(crossprod(X * sqrt(W_vec)/length(Y)))
  }

  vMatrix <- function(tau, X, Y, beta) {
    W_vec <- weightMatrix(tau, X, Y, beta)
    epsilon <- Y - X %*% beta
    V <- (sum((W_vec * epsilon)^2) *  crossprod(X)) / length(Y)
    return(V)
  }

  U <- uMatrix(tau, X, Y, beta)
  V <- vMatrix(tau, X, Y, beta)

  if (rcond(U) > 1e-10) {
    U_inv <- solve(U)
  } else {
    U_inv <- ginv(U)
  }
  return(U_inv %*% V %*% U_inv)
}

#' Renewable estimation for expectile regression
#'
#' Renewable estimation derived with Taylor approximation.
#'
#' @import expectreg MASS
#'
#' @param X A covaraites matrix.
#' @param Y A vector object.
#' @param beta Initial estimation of parameters.
#' @param Z_cum The cumulative Hessian matrix for history data.
#' @param tau Expectile level.
#'
#' @return A list object, which is basically a list consisting of:
#' \item{beta_new}{Renewable estimator}
#' \item{Z_cum}{Updated cumulative Hessian matrix}
#' \item{intbeta}{The update imputation parameters}

ReER <- function(X, Y, beta, Z_cum, tau) {
  beta_iter <- beta

  W_vec <- weightMatrix(tau, X, Y, beta_iter)
  W_mat <- crossprod(X, X * W_vec)
  Hess <- Z_cum + W_mat
  Grad <- Z_cum %*% beta_iter + crossprod(X, Y * W_vec)
  try({
    beta_new <- solve(Hess, Grad)
  }, silent = TRUE)

  if (!exists("beta_new") || any(is.na(beta_new))) {
    beta_new <- ginv(Hess) %*% Grad
  }

  W_vec <- weightMatrix(tau, X, Y, beta_new)
  W_mat <- crossprod(X, X * W_vec)
  Z_cum <- Z_cum + W_mat
  return(list(beta_new = beta_new, Z_cum = Z_cum))
}


#' Online Algorithm for Renewable Expectile Regression
#'
#' Implements an online (streaming) algorithm for expectile regression using different renewable updating strategies.
#'
#' @param X A list of covariate matrices from the streaming data, where each element corresponds to the covariates of a single data batch.
#' @param Y A list of response vectors from the streaming data, where each element corresponds to the response in a single data batch.
#' @param Xc Centralized version of \code{X}.
#' @param Yc Centralized version of \code{Y}.
#' @param K Integer. Number of data batches.
#' @param tau Numeric. Expectile level (between 0 and 1).
#' @param regFormula A formula specifying the regression model, e.g., \code{Y ~ x1 + x2}.
#' @param dataName Character vector. Variable names in the data, e.g., \code{c("Y", "x1", "x2")}.
#'
#' @return A list with three components:
#' \describe{
#'   \item{DCER}{Distributed estimation via precision matrix averaging.}
#'   \item{PAER}{Parameter aggregation using weighted local estimators.}
#'   \item{ReER}{Renewable estimation using Taylor approximation.}
#' }
#'
#' @details The function iteratively updates regression parameters for streaming data using different online algorithms. Expectile regression is estimated via the \code{expectreg.ls} function in each batch.
#'
#' @import expectreg MASS stats
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' library(expectreg)
#' library(MASS)
#' library(dplyr)
#' library(ReER)
#'
#' # generate data
#' N <- 1000; k <- 10
#' x1 <- runif(N); x2 <- runif(N)
#' Xmat <- cbind(x1, x2)
#' X <- cbind(1, Xmat); Xc <- cbind(1, scale(Xmat, scale = FALSE))
#' u <- rnorm(N)
#' beta <- c(2,1,2); gamma <- c(1, 0, 0)
#' Y <- X %*% beta + X %*% gamma * u
#' Yc <- Xc %*% beta + Xc %*% gamma * u
#'
#' # split data
#' indices <- split(1:N, cut(1:N, breaks = k, labels = FALSE))
#' X_split <- lapply(indices, function(idx) X[idx, , drop = FALSE])
#' Y_split <- lapply(indices, function(idx) Y[idx])
#' Xc_split <- lapply(indices, function(idx) Xc[idx, , drop = FALSE])
#' Yc_split <- lapply(indices, function(idx) Yc[idx])
#'
#' result <- online_alg(X = X_split, Y = Y_split,
#'                      Xc = Xc_split, Yc = Yc_split,
#'                      K = k, tau = 0.25,
#'                      regFormula = Y ~ x1 + x2,
#'                      dataName = c("Y", "x0", "x1", "x2"))

online_alg <- function(X, Y, Xc, Yc,
                       K, tau,
                       regFormula, dataName) {

  #----------- Other renewable algorthms -----------
  DCER <- function(X, Y, beta, Q_cum, Q_beta_cum, tau) {
    start_time <- Sys.time()
    Q_k <- qMatrix(X, Y, beta, tau)
    if (rcond(Q_k) > 1e-10) {
      Q_k_inv <- solve(Q_k)
    } else {
      Q_k_inv <- ginv(Q_k)
    }
    Q_cum <- Q_cum + Q_k_inv
    Q_beta_cum <- Q_beta_cum + Q_k_inv%*% beta

    try({
      beta_new <- solve(Q_cum, Q_beta_cum)
    }, silent = TRUE)

    if (!exists("beta_new") || any(is.na(beta_new))) {
      beta_new <- ginv(Q_cum) %*% Q_beta_cum
    }

    end_time <- Sys.time()
    elapsed <-  as.numeric(difftime(end_time, start_time, units = "secs"))
    return(list(beta_new = beta_new, Q_cum = Q_cum, Q_beta_cum = Q_beta_cum,
                elapsed = elapsed))
  }

  PAER <- function(X, Y, beta, Sigma_cum, Theta_cum) {
    start_time <- Sys.time()
    Sigma_cum <- Sigma_cum + crossprod(X)
    Theta_cum <- Theta_cum + crossprod(X, X %*% beta)

    try({
      beta_new <- solve(Sigma_cum, Theta_cum)
    }, silent = TRUE)

    if (!exists("beta_new") || any(is.na(beta_new))) {
      beta_new <- ginv(Sigma_cum) %*% Theta_cum
    }

    end_time <- Sys.time()
    elapsed <-  as.numeric(difftime(end_time, start_time, units = "secs"))

    return(list(beta_new = beta_new, Sigma_cum = Sigma_cum, Theta_cum = Theta_cum,
                elapsed = elapsed))
  }

  #----------- Initialization -----------
  form <- as.formula(regFormula)

  X_1 <- X[[1]]; Xc_1 <- Xc[[1]]
  Y_1 <- Y[[1]]; Yc_1 <- Yc[[1]]
  Data_1 <- as.data.frame(cbind(Y = Yc_1, Xc_1))
  dims <- ncol(X_1)
  names(Data_1) <- dataName
  m <- expectreg.ls(formula = form, data = Data_1, expectiles = tau)
  beta_1 <- c(m$intercepts, unlist(m$coefficients))

  Q_cum <- Q_beta_cum <- 0
  Q_1 <- qMatrix(X_1, Y_1, beta_1, tau)
  Q_1_inv <- ginv(Q_1)
  Q_cum <- Q_cum + Q_1_inv
  Q_beta_cum <- Q_beta_cum + Q_1_inv%*% beta_1

  Sigma_cum <- Theta_cum <- 0
  Sigma_cum <- Sigma_cum + crossprod(X_1)
  Theta_cum <- Theta_cum + crossprod(X_1, X_1 %*% beta_1)

  Z_cum <- 0
  w_vec <- weightMatrix(tau, X_1, Y_1, beta_1)
  W_mat <- crossprod(X_1, X_1 * w_vec)
  Z_cum <- Z_cum + W_mat

  #----------- Renewable estimation -----------
  ReER_path <- matrix(NA, nrow = K, ncol = dims)
  ReER_path[1, ] <- beta_1

  for (k in 2:K) {
    X_k <- X[[k]]; Xc_k <- Xc[[k]]
    Y_k <- Y[[k]]; Yc_k <- Yc[[k]]
    Data_k <- as.data.frame(cbind(Y = Yc_k, Xc_k))
    names(Data_k) <- dataName
    m <- expectreg.ls(formula = form, data = Data_k, expectiles = tau)
    beta_k <- c(m$intercepts, unlist(m$coefficients))

    # --- DCER ---
    DCER_k <- DCER(X_k, Y_k, beta_k, Q_cum, Q_beta_cum, tau)
    DCER_est <- DCER_k$beta_new %>% t()
    Q_cum <- DCER_k$Q_cum
    Q_beta_cum <- DCER_k$Q_beta_cum

    # --- PAER ---
    PAER_k <- PAER(X_k, Y_k, beta_k, Sigma_cum, Theta_cum)
    PAER_est <- PAER_k$beta_new %>% t()
    Sigma_cum <- PAER_k$Sigma_cum
    Theta_cum <- PAER_k$Theta_cum

    # --- ReER ---
    beta_old <- ReER_path[(k - 1), ]
    ReER_k <- ReER(X_k, Y_k, beta_old, Z_cum, tau = tau)
    ReER_est <- ReER_k$beta_new %>% t()
    Z_cum <- ReER_k$Z_cum
    ReER_path[k, ] <- ReER_est
  }

  beta_est <- list(DCER = DCER_est, PAER = PAER_est, ReER = ReER_est)

  return(beta_est = beta_est)
}
