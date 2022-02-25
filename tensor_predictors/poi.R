#' Penalysed Orthogonal Iteration.
#'
#' @param lambda Default: 0.75 * lambda_max for FastPOI-C method.
#'
#' @note use.C required 'poi.so' beeing dynamicaly loaded.
#'  dyn.load('../tensor_predictors/poi.so')
POI <- function(A, B, d,
                lambda = 0.75 * sqrt(max(rowSums(Delta^2))),
                update.tol = 1e-3,
                tol = 100 * .Machine$double.eps,
                maxit = 400L,
                maxit.outer = maxit,
                maxit.inner = maxit,
                use.C = FALSE,
                method = 'FastPOI-C') {

    # TODO:
    stopifnot(method == 'FastPOI-C')

    if (nrow(A) < 100) {
        Delta <- eigen(A, symmetric = TRUE)$vectors[, 1:d, drop = FALSE]
    } else {
        Delta <- try(RSpectra::eigs_sym(A, d)$vectors, silent = TRUE)
        if (is(Delta, 'try-error')) {
            Delta <- eigen(A, symmetric = TRUE)$vectors[1:d, , drop = FALSE]
        }
    }

    # Set initial value.
    Z <- Delta

    # Step 1: Optimization.
    # The "inner" optimization loop, aka repeated coordinate optimization.
    if (use.C) {
        Z <- .Call('FastPOI_C_sub', A, B, Delta, lambda, as.integer(maxit.inner),
                   PACKAGE = 'poi')
    } else {
        p <- nrow(Z)
        for (iter.inner in 1:maxit.inner) {
            Zold <- Z
            for (g in 1:p) {
                a <- Delta[g, ] - B[g, ] %*% Z + B[g, g] * Z[g, ]
                a_norm <- sqrt(sum(a^2))
                if (a_norm > lambda) {
                    Z[g, ] <- a * ((1 - lambda / a_norm) / B[g, g])
                } else {
                    Z[g, ] <- 0
                }
            }
            if (norm(Z - Zold, 'F') < update.tol) {
                break;
            }
        }
    }

    # Step 2: QR decomposition.
    if (d == 1L) {
        Z_norm <- sqrt(sum(Z^2))
        if (Z_norm < tol) {
            Q <- matrix(0, p, d)
        } else {
            Q <- Z / Z_norm
        }
    } else {
        # Detect zero columns.
        zeroColumn <- colSums(abs(Z)) < tol
        if (all(zeroColumn)) {
            Q <- matrix(0, p, d)
        } else if (any(zeroColumn)) {
            Q <- matrix(0, p, d)
            Q[, !zeroColumn] <- qr.Q(qr(Z))
        } else {
            Q <- qr.Q(qr(Z))
        }
    }

    return(list(Z = Z, Q = Q, iter.inner = if (use.C) NA else iter.inner,
                lambda = lambda))
}
