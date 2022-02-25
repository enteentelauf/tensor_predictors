source('../tensor_predictors/matpow.R')

#' Longitudinal Sliced Inverse Regression
#'
#' @param X matrix of dim \eqn{n \times p t} with each row representing a
#'  vectorized \eqn{p \times t} observation.
#' @param y vector of \eqn{n} elements as factors. (can be coersed to factors)
#' @param p,t,k,r dimensions.
#'
#' @returns a list with components
#'  alpha: matrix of \eqn{t \times r}
#'  beta: matrix of \eqn{p \times k}
#'
#' TODO: finish
#'
LSIR <- function(X, y, p, t, k = 1L, r = 1L) {
    # the code assumes:
    # alpha: T x r, beta: p x k, X_i: p x T, for ith observation
    
    # Check and transform parameters.
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- nrow(X)
    stopifnot(
        ncol(X) == p * t,
        n == length(y)
    )
    if (!is.factor(y)) y <- factor(y)

    # Restructure X into a 3D tensor with axis (observations, predictors, time).
    dim(X) <- c(n, p, t)

    # Estimate predictor/time covariance matrices \hat{Sigma}_1, \hat{Sigma}_2.
    sigma_p <- matrix(rowMeans(apply(X, 3, cov)), p, p)
    sigma_t <- matrix(rowMeans(apply(X, 2, cov)), t, t)

    # Normalize X as vec(Z) = Sigma^-1/2 (vec(X) - E(vec(X)))
    dim(X) <- c(n, p * t)
    sigma_p_isqrt <- matpow(sigma_p, -0.5)
    sigma_t_isqrt <- matpow(sigma_t, -0.5)
    Z <- scale(X, scale = FALSE) %*% kronecker(sigma_t_isqrt, sigma_p_isqrt)
    # Both as 3D tensors.
    dim(X) <- dim(Z) <- c(n, p, t)

    # Estimate the conditional predictor/time covariance matrix Omega = cov(E(Z|Y)).
    omega_p <- matrix(Reduce(`+`, lapply(levels(y), function(l) {
        rowMeans(apply(Z[y == l, , ], 3, function(z) {
            (nrow(z) / n) * tcrossprod(colMeans(z))
        }))
    })), p, p)
    omega_t <- matrix(Reduce(`+`, lapply(levels(y), function(l) {
        rowMeans(apply(Z[y == l, , ], 2, function(z) {
            (nrow(z) / n) * tcrossprod(colMeans(z))
        }))
    })), t, t)
    omega <- kronecker(omega_t, omega_p)

    # Compute seperate SVD of estimated omega's and use that for an estimate of
    # a central subspace basis.
    svd_p <- La.svd(omega_p)
    svd_t <- La.svd(omega_t)
    beta  <- sigma_p_isqrt %*% svd_p$u[, k]
    alpha <- sigma_t_isqrt %*% svd_t$u[, r]

    return(list(sigma_p = sigma_p, sigma_t = sigma_t,
                sigma = kronecker(sigma_t, sigma_p),
                alpha = alpha, beta = beta,
                Delta = omega,
                B = kronecker(alpha, beta)))
}
