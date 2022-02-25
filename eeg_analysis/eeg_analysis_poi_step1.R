# Source Code.                                    # Loaded functions.
source('../tensor_predictors/poi.R')              # POI

# Load C implentation of 'FastPOI-C' subroutine.
# Required for using 'use.C = TRUE' in the POI method.
# Compiled via.
# $ cd ../tensor_predictors/
# $ R CMD SHLIB poi.c
dyn.load('../tensor_predictors/poi.so')
# dyn.load('../tensor_predictors/poi.dll') # On Windows
# In this case 'use.C = TRUE' is required cause the R implementation is not
# sufficient due to memory exhaustion (and runtime).

# Load Dataset.
# > dataset <- read.table(file = 'egg.extracted.means.txt', header = TRUE,
# >                       stringsAsFactors = FALSE, check.names = FALSE)
# Save as Rdata file for faster loading.
# > saveRDS(dataset, file = 'eeg_data.rds')
dataset <- readRDS('../data_analysis/eeg_data.rds')

# Positive and negative case index.
set.seed(42)
zero <- sample(which(dataset$Case_Control == 0))
one  <- sample(which(dataset$Case_Control == 1))

# 10-fold test groups.
zero <- list(zero[ 1: 4], zero[ 5: 8], zero[ 9:12], zero[13:16],
             zero[17:20], zero[21:25], zero[26:30],
             zero[31:35], zero[36:40], zero[41:45])

one <- list(one[ 1: 8], one[ 9:16], one[17:24], one[25:32],
            one[33:40], one[41:48], one[49:56],
            one[57:63], one[64:70], one[71:77])

# Iterate data folds.
folds <- vector('list', 10)
for (i in seq_along(folds)) {
    cat('\r%d/%d ', i, length(folds))

    # Call garbage collector.
    gc()

    # Formulate PFC-GEP for EEG data.
    index <- c(zero[[i]], one[[i]])
    X  <- scale(dataset[-index, -(1:2)], scale = FALSE, center = TRUE)
    Fy <- scale(dataset$Case_Control[-index], scale = FALSE, center = TRUE)
    B <- crossprod(X) / nrow(X)                     # Sigma
    P_Fy <- Fy %*% solve(crossprod(Fy), t(Fy))
    A <- crossprod(X, P_Fy %*% X) / nrow(X)         # Sigma_fit

    # Before Starting POI on (very big GEP) call the garbage collector.
    gc()
    poi <- POI(A, B, 1L, lambda = lambda, use.C = TRUE)
    rm(A, B)
    gc()

    # Set fold index.
    poi$index = index

    folds[[i]] <- poi
}
cat('\n')

# Save complete 10 fold results.
file <- sprintf('eeg_analysis_poi.rds')
saveRDS(folds, file = file)
