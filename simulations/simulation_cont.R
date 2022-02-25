source('../tensor_predictors/random.R')
source('../tensor_predictors/multi_assign.R')
source('../tensor_predictors/tensor_predictors.R')
source('../tensor_predictors/lsir.R')
source('../tensor_predictors/pca2d.R')

simulateData.cont <- function(n, p, t, k, r, d1, d2, delta.identity = FALSE) {

    stopifnot(d1 <= r, d2 <= k)

    y  <- rnorm(n)
    ns <- r * k / 2
    Fy <- do.call(cbind, lapply(1:ns, function(s, z) {
        cbind(cos(s * z), sin(s * z))
    }, z = 2 * pi * y))
    Fy <- scale(Fy, scale = FALSE)

    Gamma_1 <- diag(1, t, d1)
    gamma_1 <- diag(1, d1, r)
    alpha   <- Gamma_1 %*% gamma_1
    Gamma_2 <- diag(1, p, d2)
    gamma_2 <- diag(1, d2, k)
    beta    <- Gamma_2 %*% gamma_2

    if (delta.identity) {
        Delta    <- diag(1, p * t, p * t)
    } else {
        Delta    <- crossprod(matrix(rnorm((p * t)^2), p * t))
        DM_Delta <- diag(sqrt(1 / diag(Delta)))
        Delta    <- DM_Delta %*% Delta %*% DM_Delta
    }

    X <- tcrossprod(Fy, kronecker(alpha, beta)) + rmvnorm(n, sigma = Delta)
    X <- scale(X, scale = FALSE)

    return(list(X = X, y = y, Fy = Fy,
                Gamma = kronecker(Gamma_1, Gamma_2),
                Gamma_1 = Gamma_1, gamma_1 = gamma_1, alpha = alpha,
                Gamma_2 = Gamma_2, gamma_2 = gamma_2, beta = beta,
                Delta = Delta))
}

simulation.cont <- function(methods, reps, n, p, t, k, r, d1, d2) {
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
        ds <- simulateData.cont(n, p, t, k, r, d1, d2)
        for (method.name in names(methods)) {
            cat(sprintf('\r%4d/%d in %s', rep, reps, method.name))

            method <- methods[[method.name]]
            sdr <- method(ds$X, ds$Fy, p, t, k, r, d1, d2)
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
            sd.E1   = sd(unlist(E1[m])),
            mean.E2 = mean(unlist(E2[m])),
            sd.E2   = sd(unlist(E2[m])),
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
    attr(stat, "params") <- c(reps = reps, n = n, p = p, t = t, k = k, r = r,
                              d1 = d1, d2 = d2)
    return(stat)
}

methods <- list(
    KPIR_LS   = function(...) tensor_predictor(..., method = "KPIR_LS"),
    KPIR_MLE  = function(...) tensor_predictor(..., method = "KPIR_MLE"),
    KPFC1     = function(...) tensor_predictor(..., method = "KPFC1"),
    KPFC2     = function(...) tensor_predictor(..., method = "KPFC2"),
    KPFC3     = function(...) tensor_predictor(..., method = "KPFC3"),
    PCA2d     = function(X, y = NULL, p, t, k = 1L, r = 1L, d1 = 1L, d2 = 1L) {
        pca <- PCA2d(X, p, t, k, r)
        # Note: alpha, beta are not realy meaningfull for (d1, d2) != (r, k)
        pca$Gamma_1 <- pca$alpha[, 1:d1, drop = FALSE]
        pca$Gamma_2 <- pca$beta[, 1:d2, drop = FALSE]
        pca$Gamma <- kronecker(pca$Gamma_1, pca$Gamma_2)
        pca$Delta <- kronecker(pca$Sigma_t, pca$Sigma_p)
        return(pca)
    }
)

#                    n,  p,  t,  k,  r, d1, d2
#              -----------------------------
params <- list( c( 500, 10,  8,  6,  6,  6,  6)
              , c( 500, 10,  8,  6,  6,  4,  4)
              , c( 500, 10,  8,  6,  6,  2,  2)
              , c(5000, 10,  8,  6,  6,  6,  6)
              , c(5000, 10,  8,  6,  6,  4,  4)
              , c(5000, 10,  8,  6,  6,  2,  2)
)

for (param in params) {
    c(n, p, t, k, r, d1, d2) %<-% param
    sim <- simulation.cont(methods, 500, n, p, t, k, r, d1, d2)

    print(attr(sim, "params"))
    print(round(sim, 2))

    saveRDS(sim, file = sprintf("simulation_cont_%d_%d_%d_%d_%d_%d_%d.rds",
                                n, p, t, k, r, d1, d2))
}
