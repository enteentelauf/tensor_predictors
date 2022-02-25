source('../tensor_predictors/random.R')
source('../tensor_predictors/multi_assign.R')
source('../tensor_predictors/tensor_predictors.R')
source('../tensor_predictors/lsir.R')
source('../tensor_predictors/pca2d.R')

#' @param n0 number of controls
#' @param n1 number of cases
simulateData.binary <- function(n0, n1, p, t, rho.p, rho.t) {
    # Response vector
    Y <- c(rep(1, n1), rep(0, n0))

    # Section 7.1.2 of Tensor_Predictors-4.pdf
    alpha0 <- as.matrix(rep(0, t))
    alpha1 <- as.matrix(1 / ((t + 1) - 1:t))
    beta   <- as.matrix(rep(1 / sqrt(p), p))
    mu0    <- kronecker(alpha0, beta)
    mu1    <- kronecker(alpha1, beta)

    sigma1 <- rho.t^abs(outer(1:t, 1:t, FUN = `-`))
    sigma2 <- rho.p^abs(outer(1:p, 1:p, FUN = `-`))
    sigma  <- kronecker(sigma1, sigma2)

    # Compute Delta
    # Delta = Sigma + E[vec(X)]E[vec(X)^t] - E{E[vec(X)|Y]E[vec(X)^t|Y]}
    n      <- n0 + n1
    muAvg  <- (n0 * mu0 + n1 * mu1) / n
    mat0   <- mu0 %*% t(mu0)
    mat1   <- mu1 %*% t(mu1)
    matAvg <- (n0 * mat0 + n1 * mat1) / n
    Delta  <- sigma + (muAvg %*% t(muAvg)) - matAvg

    X1     <- rmvnorm(n1, mu1, Delta)
    X0     <- rmvnorm(n0, mu0, Delta)
    X      <- rbind(X1, X0)

    # Center data
    Y <- scale(Y, center = TRUE, scale = FALSE)
    X <- scale(X, center = TRUE, scale = FALSE)

    alpha <- alpha0 - alpha1
    Gamma_1 <- alpha / norm(alpha, 'F')
    Gamma_2 <- beta / norm(beta, 'F')
    list(Y = Y, X = X,
         Gamma_1 = Gamma_1, Gamma_2 = Gamma_2,
         Gamma = kronecker(Gamma_1, Gamma_2),
         alpha = alpha, beta = beta,
         Delta = Delta
    )
}


simulation.binary <- function(methods, reps, n0, n1, p, t, rho.p, rho.t) {
    nsim <- length(methods) * reps
    results <- vector('list', nsim)
    E1      <- vector('list', nsim)
    E2      <- vector('list', nsim)
    vec1    <- vector('list', nsim)
    vec2    <- vector('list', nsim)
    Phi     <- vector('list', nsim)
    phi1    <- vector('list', nsim)
    phi2    <- vector('list', nsim)

    i <- 1
    for (rep in 1:reps) {
        set.seed(rep)
        ds <- simulateData.binary(n0, n1, p, t, rho.p, rho.t)
        for (method.name in names(methods)) {
            cat(sprintf('\r%4d/%d in %s', rep, reps, method.name))

            method <- methods[[method.name]]
            sdr <- method(ds$X, ds$Y, p, t)
            # Store which silumation is at index i.
            results[[i]] <- c(method = method.name, rep = rep)
            # Compute simpulation validation metrics.
            E1[[i]] <-
                norm(kronecker(ds$alpha, ds$beta) - kronecker(sdr$alpha, sdr$beta), 'F') /
                    norm(kronecker(ds$alpha, ds$beta), 'F')
            E2[[i]]   <- norm(ds$Delta - sdr$Delta, 'F') / norm(ds$Delta, 'F')
            vec1[[i]] <- as.double(kronecker(sdr$alpha, sdr$beta))
            vec2[[i]] <- as.double(sdr$Delta)
            # Subspace distances.
            if (!('Gamma' %in% names(sdr))) {
                # Assuming r = k = 1
                sdr$Gamma_1 <- sdr$alpha / norm(sdr$alpha, 'F')
                sdr$Gamma_2 <- sdr$beta / norm(sdr$beta, 'F')
                sdr$Gamma <- kronecker(sdr$Gamma_1, sdr$Gamma_2)
            }
            Phi[[i]]  <- norm(tcrossprod(ds$Gamma) - tcrossprod(sdr$Gamma), 'F')
            phi1[[i]] <- norm(tcrossprod(ds$Gamma_1) - tcrossprod(sdr$Gamma_1), 'F')
            phi2[[i]] <- norm(tcrossprod(ds$Gamma_2) - tcrossprod(sdr$Gamma_2), 'F')
            i <- i + 1
        }
    }
    cat('\n')

    # Aggregate per method statistics.
    statistics <- list()
    for (method.name in names(methods)) {
        m <- which(unlist(lapply(results, `[`, 1)) == method.name)

        # Convert list of vec(alpha %x% beta) to a matrix with vec(alpha %x% beta)
        # in its columns.
        tmp <- matrix(unlist(vec1[m]), ncol = length(m))
        V1 <- sum(apply(tmp, 1, var))

        # Convert list of vec(Delta) to a matrix with vec(Delta) in its columns.
        tmp <- matrix(unlist(vec2[m]), ncol = length(m))
        V2 <- sum(apply(tmp, 1, var))

        statistics[[method.name]] <- list(
            mean.E1 = mean(unlist(E1[m])),
            sd.E1 = sd(unlist(E1[m])),
            mean.E2 = mean(unlist(E2[m])),
            sd.E2 = sd(unlist(E2[m])),
            V1 = V1,
            V2 = V2,
            Phi  = mean(unlist(Phi[m])),
            phi1 = mean(unlist(phi1[m])),
            phi2 = mean(unlist(phi2[m]))
        )
    }
    # transform the statistics list into a data.frame with row and col names.
    stat <- t(matrix(unlist(statistics), ncol = length(statistics)))
    rownames(stat) <- names(statistics)
    colnames(stat) <- names(statistics[[1]])
    stat <- as.data.frame(stat)
    attr(stat, "params") <- c(reps = reps, n0 = n0, n1 = n1, p = p, t = t,
                              rho.p = rho.p, rho.t = rho.t)
    return(stat)
}

methods <- list(
    KPIR_LS   = function(...) tensor_predictor(..., method = "KPIR_LS"),
    KPIR_MLE  = function(...) tensor_predictor(..., method = "KPIR_MLE"),
    KPFC1     = function(...) tensor_predictor(..., method = "KPFC1"),
    KPFC2     = function(...) tensor_predictor(..., method = "KPFC2"),
    KPFC3     = function(...) tensor_predictor(..., method = "KPFC3"),
    LSIR      = function(X, Fy, p, t) LSIR(X, Fy, p, t, k = 1, r = 1),
    PCA2d     = function(X, y = NULL, p, t, k = 1, r = 1, d1 = 1, d2 = 1) {
        pca <- PCA2d(X, p, t, k, r)
        pca$Gamma_1 <- pca$alpha[, 1:d1, drop = FALSE]
        pca$Gamma_2 <- pca$beta[, 1:d2, drop = FALSE]
        pca$Gamma <- kronecker(pca$Gamma_1, pca$Gamma_2)
        pca$Delta <- kronecker(pca$Sigma_t, pca$Sigma_p)
        return(pca)
    }
)

#                   n0,   n1,  p,  t, rho.p, rho.t
#              -----------------------------------
params <- list( c( 250,  250, 10,  5,   0.3,   0.3)
              , c( 500,  500, 10,  5,   0.3,   0.3)
              , c(1000, 1000, 10,  5,   0.3,   0.3)
)

for (param in params) {
    c(n0, n1, p, t, rho.p, rho.t) %<-% param
    sim <- simulation.binary(methods, 500, n0, n1, p, t, rho.p, rho.t)

    print(attr(sim, "params"))
    print(round(sim, 2))

    saveRDS(sim, file = sprintf("simulation_3_desc_%d_%d_%d_%d_%f_%f.rds",
                                n0, n1, p, t, rho.p, rho.t))
}
