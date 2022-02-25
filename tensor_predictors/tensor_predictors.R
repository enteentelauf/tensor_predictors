source('../tensor_predictors/matpow.R')
source('../tensor_predictors/multi_assign.R')
source('../tensor_predictors/approx_kronecker.R')

log.likelihood <- function(par, X, Fy, Delta.inv, da, db) {
    alpha <- matrix(par[1:prod(da)], da[1L])
    beta <- matrix(par[(prod(da) + 1):length(par)], db[1L])
    error <- X - tcrossprod(Fy, kronecker(alpha, beta))
    sum(error * (error %*% Delta.inv))
}

tensor_predictor <- function(X, Fy, p, t, k = 1L, r = 1L, d1 = 1L, d2 = 1L,
                             method = "KPIR_LS",
                             eps1 = 1e-2, eps2 = 1e-2, maxit = 10L) {
    # Validate method using unexact matching.
    methods <- list(KPIR_LS = "KPIR_LS", KPIR_MLE = "KPIR_MLE",
                    KPFC1 = "KPFC1", KPFC2 = "KPFC2", KPFC3 = "KPFC3")
    method <- methods[[toupper(method), exact = FALSE]]
    if (is.null(method)) {
        stop("Unable to determine method.")
    }

    if (method %in% c("KPIR_LS", "KPIR_MLE")) {
        ## Step 1:
        # OLS estimate of the model `X = F_y B + epsilon`.
        B <- t(solve(crossprod(Fy), crossprod(Fy, X)))

        # Estimate alpha, beta as nearest kronecker approximation.
        c(alpha, beta) %<-% approx.kronecker(B, c(t, r), c(p, k))

        if (method == "KPIR_LS") {
            # Estimate Delta.
            B <- kronecker(alpha, beta)
            rank <- if (ncol(Fy) == 1) 1L else qr(Fy)$rank
            Delta <- crossprod(X - tcrossprod(Fy, B)) / (nrow(X) - rank)

        } else { # KPIR_MLE
            # Estimate initial Delta.
            B <- kronecker(alpha, beta)
            Delta <- crossprod(X - tcrossprod(Fy, B)) / nrow(X)

            for (. in 1:maxit) {
                # Optimize log-likelihood for alpha, beta with fixed Delta.
                opt <- optim(c(alpha, beta), log.likelihood, gr = NULL,
                             X, Fy, matpow(Delta, -1), c(t, r), c(p, k))
                # Store previous alpha, beta and Delta (for break consition).
                Delta.last <- Delta
                B.last     <- B
                # Extract optimized alpha, beta.
                alpha <- matrix(opt$par[1:(t * r)], t, r)
                beta <- matrix(opt$par[(t * r + 1):length(opt$par)], p, k)
                # Calc new Delta with likelihood optimized alpha, beta.
                B <- kronecker(alpha, beta)
                Delta <- crossprod(X - tcrossprod(Fy, B)) / nrow(X)
                # Check break condition 1.
                if (norm(Delta - Delta.last, 'F') < eps1 * norm(Delta, 'F')) {
                    # Check break condition 2.
                    if (norm(B - B.last, 'F') < eps2 * norm(B, 'F')) {
                        break
                    }
                }
            }
        }

        # Construct basis from alpha and beta.
        Gamma_1 <- if(d1 > 1L) La.svd(alpha, d1, 0L)$u
                   else        alpha / norm(alpha, 'F')
        Gamma_2 <- if(d2 > 1L) La.svd(beta,  d2, 0L)$u
                   else        beta / norm(beta, 'F')
        Gamma <- kronecker(Gamma_1, Gamma_2)
    } else if (method %in% c("KPFC1", "KPFC2", "KPFC3")) {
        ## Step 1:
        # OLS extimate of the model `X = F_y B + epsilon`.
        B <- t(solve(crossprod(Fy), crossprod(Fy, X)))

        ## Step 2:
        # Estimate Delta_mle.
        P_Fy <- Fy %*% solve(crossprod(Fy), t(Fy))
        Q_Fy <- diag(nrow(P_Fy)) - P_Fy
        Delta_fit <- crossprod(X, P_Fy %*% X) / nrow(X)
        Delta_res <- crossprod(X, Q_Fy %*% X) / nrow(X)
        # Compute Delta_mle using equation (7).
        D <- matpow(Delta_res, -0.5)
        Delta <- with(La.svd(D %*% Delta_fit %*% D), {
            K <- diag(c(rep(0, d1 * d2), d[-(1:(d1 * d2))]))
            D <- matpow(Delta_res, 0.5)
            Delta_res + (D %*% u %*% tcrossprod(K, u) %*% D)
        })

        ## Step 3:
        # Set Gamma to be the first `d = d1 * d2` eigenvectors of (25).
        D <- matpow(Delta, -0.5)
        Gamma <- with(La.svd(D %*% Delta_fit %*% D, d1 * d2, 0L), {
            La.svd(matpow(Delta, 0.5) %*% u[, 1:(d1 * d2)])$u
        })

        if (method == "KPFC1") {
            # Compute lower_gamma using (26).
            D <- crossprod(Gamma, matpow(Delta, -1))
            lower_gamma <- solve(D %*% Gamma, D %*% B)

            ## Step 4a:
            # Calc MLE estimate of B.
            B <- Gamma %*% lower_gamma
            # Using the VLP approx. for a kronecker product factorization.
            c(alpha, beta) %<-% approx.kronecker(B, c(t, r), c(p, k))

            # Construct basis from alpha and beta.
            Gamma_1 <- if(d1 > 1L) La.svd(alpha, d1, 0L)$u
                       else        alpha / norm(alpha, 'F')
            Gamma_2 <- if(d2 > 1L) La.svd(beta,  d2, 0L)$u
                       else        beta / norm(beta, 'F')
            Gamma <- kronecker(Gamma_1, Gamma_2)

        } else { # KPFC2, KPFC3
            ## Step 4b:
            # Estimate Gamma's as nearest kronecker approximation of Gamma.
            c(Gamma_1, Gamma_2) %<-% approx.kronecker(Gamma, c(t, d1), c(p, d2))
            Gamma <- kronecker(Gamma_1, Gamma_2)
            # Compute lower_gamma using (26).
            D <- crossprod(Gamma, matpow(Delta, -1))
            lower_gamma <- solve(D %*% Gamma, D %*% B)

            if (prod(dim(lower_gamma)) == 1) {
                # If lower_gamma is a scalar, then alpha, beta is only scaled.
                # (shortcut)
                lg1 <- lg2 <- sqrt(abs(as.vector(lower_gamma)))
                alpha <- lg1 * Gamma_1
                beta <- lg2 * Gamma_2
            } else if (method == "KPFC2") {
                ## Step 5c:
                c(alpha, beta) %<-% approx.kronecker(Gamma %*% lower_gamma,
                                                     c(t, r), c(p, k))
            } else { # KPFC3
                ## Step 5d:
                c(lg1, lg2) %<-% approx.kronecker(lower_gamma,
                                                  c(d1, r), c(d2, k))
                alpha <- Gamma_1 %*% lg1
                beta <- Gamma_2 %*% lg2
            }
        }
    }

    return(structure(
        list(alpha = alpha,
             beta = beta,
             Gamma = Gamma,
             Gamma_1 = Gamma_1, Gamma_2 = Gamma_2,
             Delta = Delta),
        class = c("tensor_predictor", method)
    ))
}

#' TODO: Write this properly!
reduce <- function(object, data, use = 'Gamma') {
    if (use == 'Gamma') {
        projection <- object$Gamma
    } else if (use == 'alpha_beta') {
        projection <- kronecker(object$alpha, object$beta)
    } else {
        stop("Unkown 'use' parameter value.")
    }

    # ensure alignement of multiple calls.
    if (projection[1] < 0) {
        projection <- -projection
    }

    return(data %*% projection)
}
