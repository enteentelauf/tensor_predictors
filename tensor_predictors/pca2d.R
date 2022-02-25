
#' @param X Matrix of dim (n, p * t) with each row the vectorized (p, t) observation.
#' @param p nr. predictors
#' @param t nr. timepoints
#' @param ppc reduced nr. predictors (p-principal components)
#' @param tpc reduced nr. timepoints (t-principal components)
#'
#' @details The `i`th observation is stored in a row such that its matrix equiv
#'  is given by `matrix(X[i, ], p, t)`.
#' 
PCA2d <- function(X, p, t, ppc, tpc, scale = FALSE) {
    stopifnot(ncol(X) == p * t, ppc <= p, tpc <= t)

    X <- scale(X, center = TRUE, scale = scale)

    # Left/Right aka predictor/time covariance matrices.
    dim(X) <- c(nrow(X), p, t)
    Sigma_p <- matrix(apply(apply(X, 1, tcrossprod), 1, mean), p, p) # Sigma_beta
    Sigma_t <- matrix(apply(apply(X, 1,  crossprod), 1, mean), t, t) # Sigma_alpha
    dim(X) <- c(nrow(X), p * t)

    V_p <- La.svd(Sigma_p, ppc, 0)$u
    V_t <- La.svd(Sigma_t, tpc, 0)$u

    X <- X %*% kronecker(V_t, V_p)

    return(list(reduced = X, alpha = V_t, beta = V_p,
                Sigma_t = Sigma_t, Sigma_p = Sigma_p))
}
