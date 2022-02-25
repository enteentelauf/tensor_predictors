# implementation contains fallback if the package is not available but for this
# case required!
library(RSpectra)

# Load POI function and compiled C subroutine.
source('../tensor_predictors/poi.R')
dyn.load('../tensor_predictors/poi.so') # "Shared Object" of POI-Subrountine

# Load data from sent data file (last Email)
dataset <- readRDS('../eeg_analysis/eeg_data.rds')

maxit <- 400L # Upper bound for number of optimization iterations.

for (i in 1:nrow(dataset)) {
    gc() # To be on the save side, call the garbage collector (free memory)

    # Formulate PFC-GEP (Principal Fitted Components - Generalized Eigenvalue
    # Problem) for EEG data.
    X  <- scale(dataset[-i, -(1:2)], scale = FALSE, center = TRUE)
    Fy <- scale(dataset$Case_Control[-i], scale = FALSE, center = TRUE)
    B <- crossprod(X) / nrow(X)                     # Sigma
    P_Fy <- Fy %*% solve(crossprod(Fy), t(Fy))
    A <- crossprod(X, P_Fy %*% X) / nrow(X)         # Sigma_fit

    # Call POI using C subroutine (requires "dyn.load" of subroutine)
    poi_res <- POI(A, B, 1L, maxit = maxit, use.C = TRUE)
    # Again, be nice to memory and delete with an explicit fall to gc.
    rm(A, B)
    gc()

    # Store results, do analysis, ... (addapt to needs) .
    poi_res$maxit = maxit
    poi_res$loo_index = i # Keep track of LOO position.

    # Save i-th LOO result to file for analysis/validation/visualization/...
    saveRDS(poi_res, file = sprintf('eeg_poi_loo_%d.rds', i))
}
