# Source Code.                                    # Loaded functions.
source('../tensor_predictors/multi_assign.R')     # %<-%
source('../tensor_predictors/approx_kronecker.R') # approx_kronecker
source('../tensor_predictors/poi.R')              # POI
source('../tensor_predictors/subspace.R')         # subspace
source('../tensor_predictors/random.R')           # rmvnorm

# Load C impleentation of 'FastPOI-C' subroutine.
# Required for using 'use.C = TRUE' in the POI method.
dyn.load('../tensor_predictors/poi.so')
# When 'use.C = FALSE' the POI method uses a base R implementation.
use.C = TRUE

simulateData.sparse <- function(n, p, t, k, r, scale, degree = 2) {
    # Define true reduction matrices alpha, beta.
    alpha <- diag(1, t, r)
    beta  <- diag(1, p, k)

    # Create true "random" covariance of inverse model.
    R <- matrix(rnorm((p * t)^2), p * t)        # random square matrix.
    sigma <- tcrossprod(R / sqrt(rowSums(R^2))) # sym. pos.def. with diag = 1.

    # Sample responces.
    y <- rnorm(n, 0, 1)
    # equiv to cbind(y^1, y^2, ..., y^degree)
    Fy <- t(vapply(y, `^`, double(degree), seq(degree)))

    # Calc X according the inverse regression model.
    X <- tcrossprod(scale(Fy, scale = FALSE, center = TRUE), kronecker(alpha, beta))
    X <- X + (scale * rmvnorm(n, sigma = sigma))

    return(list(X = X, y = y, Fy = Fy, alpha = alpha, beta = beta))
}

# # True Positives Rate
# tpr <- function(Y, Y_hat) {
#     sum(as.logical(Y_hat) & as.logical(Y)) / sum(as.logical(Y)) # TP / P
# }
# False Positives Rate
fpr <- function(Y, Y_hat) {
    sum(as.logical(Y_hat) & !Y) / sum(!Y) # FP / N
}
# False Negative Rate
fnr <- function(Y, Y_hat) {
    sum(!Y_hat & as.logical(Y)) / sum(as.logical(Y)) # FN / P
}
# False Rate (rate of false positives and negatives)
fr <- function(Y, Y_hat) {
    sum(as.logical(Y) != as.logical(Y_hat)) / length(Y)
}

simulation.sparse <- function(scales, reps, n, p, t, k, r,
                              eps = 100 * .Machine$double.eps) {
    results <- vector('list', length(scales) * reps)

    i <- 0
    for (scale in scales) {
        for (rep in 1:reps) {
            cat(sprintf('\r%4d/%d for scale = %.2f', rep, reps, scale))

            ds <- simulateData.sparse(n, p, t, k, r, scale)
            # Formulate PFC-GEP for given dataset.
            X  <- scale(ds$X,  scale = FALSE, center = TRUE)
            Fy <- scale(ds$Fy, scale = FALSE, center = TRUE)
            Sigma <- crossprod(X) / nrow(X)
            P_Fy <- Fy %*% solve(crossprod(Fy), t(Fy))
            Sigma_fit <- crossprod(X, P_Fy %*% X) / nrow(X)

            poi <- POI(Sigma_fit, Sigma, k * r, use.C = use.C)
            # Calc approx. alpha, beta and drop further drop "zero" from konecker
            # factorization approximation.
            c(alpha, beta) %<-% approx.kronecker(poi$Q, dim(ds$alpha), dim(ds$beta))
            alpha[abs(alpha) < eps] <- 0
            beta[abs(beta) < eps] <- 0

            # Compair estimates against true alpha, beta.
            result <- list(
                scale = scale,
                lambda = poi$lambda,
                # alpha_tpr = tpr(ds$alpha, alpha),
                alpha_fpr = fpr(ds$alpha, alpha),
                alpha_fnr = fnr(ds$alpha, alpha),
                alpha_fr = fr(ds$alpha, alpha),
                # beta_tpr = tpr(ds$beta, beta),
                beta_fpr = fpr(ds$beta, beta),
                beta_fnr = fnr(ds$beta, beta),
                beta_fr = fr(ds$beta, beta)
            )
            # Component-wise validation (_c_ stands for component)
            if (ncol(alpha) > 1) {
                ds_c_alpha <- apply(!!ds$alpha, 1, any)
                   c_alpha <- apply(!!   alpha, 1, any)
                # result$alpha_c_tpr <- tpr(ds_c_alpha, c_alpha)
                result$alpha_c_fpr <- fpr(ds_c_alpha, c_alpha)
                result$alpha_c_fnr <- fnr(ds_c_alpha, c_alpha)
                result$alpha_c_fr <- fr(ds_c_alpha, c_alpha)
            }
            if (ncol(beta) > 1) {
                ds_c_beta <- apply(!!ds$beta, 1, any)
                   c_beta <- apply(!!   beta, 1, any)
                # result$beta_c_tpr <- tpr(ds_c_beta, c_beta)
                result$beta_c_fpr <- fpr(ds_c_beta, c_beta)
                result$beta_c_fnr <- fnr(ds_c_beta, c_beta)
                result$beta_c_fr <- fr(ds_c_beta, c_beta)
            }
            results[[i <- i + 1]] <- result
        }
        cat('\n')
    }

    # Restructure results list of lists as data.frame.
    results <- as.data.frame(t(sapply(results, function(res, cols) {
        unlist(res[cols])
    }, names(results[[1]]))))
    results$scale <- as.factor(results$scale)
    attr(results, 'params') <- list(
        reps = reps, n = n, p = p, t = t, k = k, r = r, eps = eps)

    results
}

reps <- 500
#                   n,  p,  t, k, r
#              --------------------
params <- list( c(100, 10,  5, 1, 2)
              , c(100,  7,  5, 1, 2)
              , c(100,  5,  3, 1, 2)
              , c(500, 10,  5, 1, 2)
              , c(500,  7,  5, 1, 2)
              , c(500,  5,  3, 1, 2)
)
scales <- seq(0.5, 6, 0.25)

for (param in params) {
    c(n, p, t, k, r) %<-% param
    results <- simulation.sparse(scales, reps, n, p, t, k, r)
    sim <- aggregate(results[, 'scale' != names(results)],
                     by = list(scale = results$scale), mean)
    attr(sim, 'params') <- attr(results, 'params')

    file.name <- sprintf("simulation_sparse_%d_%d_%d_%d_%d.rds", n, p, t, k, r)
    saveRDS(sim, file = file.name)

    cat(file.name, '\n')
    print(sim, digits = 2)
}
