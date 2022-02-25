suppressPackageStartupMessages({
    library(pROC)
})

source('../tensor_predictors/approx_kronecker.R')
source('../tensor_predictors/multi_assign.R')
source('../tensor_predictors/tensor_predictors.R')
source('../tensor_predictors/lsir.R')
source('../tensor_predictors/pca2d.R')

# acc: Accuracy. P(Yhat = Y). Estimated as: (TP+TN)/(P+N).
acc <- function(y_true, y_pred) mean(round(y_pred) == y_true)
# err: Error rate. P(Yhat != Y). Estimated as: (FP+FN)/(P+N).
err <- function(y_true, y_pred) mean(round(y_pred) != y_true)
# fpr: False positive rate. P(Yhat = + | Y = -). aliases: Fallout.
fpr <- function(y_true, y_pred) mean((round(y_pred) == 1)[y_true == 0])
# tpr: True positive rate.  P(Yhat = + | Y = +). aliases: Sensitivity, Recall.
tpr <- function(y_true, y_pred) mean((round(y_pred) == 1)[y_true == 1])
# fnr: False negative rate. P(Yhat = - | Y = +). aliases: Miss.
fnr <- function(y_true, y_pred) mean((round(y_pred) == 0)[y_true == 1])
# tnr: True negative rate.  P(Yhat = - | Y = -).
tnr <- function(y_true, y_pred) mean((round(y_pred) == 0)[y_true == 0])

# Load EEG dataset
dataset <- readRDS('eeg_data.rds')

#' @param ppc Number of "p"redictor "p"rincipal "c"omponents.
#' @param tpc Number of "t"ime "p"rincipal "c"omponents.
egg_analysis_reduced <- function(methods, ppc, tpc) {
    # Set dimenional parameters.
    n <- nrow(dataset)  # sample size (nr. of people)
    p <- 64L            # nr. of predictors (count of sensorce)
    t <- 256L           # nr. of time points (measurements)

    # Extract dimension names from X.
    nNames <- dataset$PersonID
    tNames <- as.character(seq(t))
    pNames <- unlist(strsplit(colnames(dataset)[2 + t * seq(p)], '_'))[c(T, F)]

    # Split into X-y.
    X <- as.matrix(dataset[, -(1:2)])
    y <- dataset$Case_Control
    # Reshape X as 3D tenros of shape (n, t, p) aka. samples, timesteps, predictors.
    # (Each of the n rows in X iterate over the time bevore switching sensorce.)
    X <- array(X, dim = c(n, t, p),
                dimnames = list(nNames, tNames, pNames))
    # Reorder axis to (p, t, n) = (predictors, timesteps, samples).
    X <- aperm(X, c(3, 2, 1))

    # Compute Mean of X.
    X_mean <- apply(X, c(1, 2), mean)
    X_center <- X - as.vector(X_mean)

    # Compute "left" and "right" cov-matrices.
    Sigma_t <- matrix(apply(apply(X_center, 3,  crossprod), 1, mean), t, t)
    Sigma_p <- matrix(apply(apply(X_center, 3, tcrossprod), 1, mean), p, p)
    # Get "left", "right" principal components.
    V_p <- svd(Sigma_p, ppc, 0L)$u
    V_t <- svd(Sigma_t, tpc, 0L)$u

    # Reduce dimension.
    X_reduced <- apply(X_center, 3, function(x) crossprod(V_p, x %*% V_t))
    dim(X_reduced) <- c(ppc, tpc, n)

    # Vectorize to shape of (predictors * timesteps, samples) and transpose to
    # (samples, predictors * timesteps).
    X_vec <- t(matrix(X_reduced, ppc * tpc, n))

    loo.cv <- expand.grid(method = names(methods), fold = 1:n)
    loo.cv$y_true <- y[loo.cv$fold]
    loo.cv$y_pred <- NA

    # Performe LOO cross-validation for each method.
    for (i in 1L:n) {
        # Print progress.
        cat(sprintf("\rCross-Validation (p-PC: %d, t-PC: %d): %4d/%d",
                    ppc, tpc, i, n))
        # Leave Out the i-th element.
        X_train <- X_vec[-i, ]
        X_test <- X_vec[i, ]
        y_train <- y[-i]
        # Center y.
        y_train <- scale(y_train, center = TRUE, scale = FALSE)

        # For each method.
        for (method.name in names(methods)) {
            method <- methods[[method.name]]
            # Compute reduction using current method under common API.
            sdr <- method(X_train, y_train, ppc, tpc)
            B <- kronecker(sdr$alpha, sdr$beta)
            # Fit a linear model (which ensures a common sdr direction if possible).
            model <- glm(y ~ x, family = binomial(link = "logit"),
                        data = data.frame(y = y[-i], x = X_train %*% B))
            # Predict out of sample and store in LOO CV data.frame.
            y_pred <- predict(model, data.frame(x = X_test %*% B), type = "response")
            loo.cv[loo.cv$method == method.name & loo.cv$fold == i, 'y_pred'] <- y_pred
        }
    }

    for (method.name in names(methods)) {
        labels <- loo.cv[loo.cv$method == method.name, 'y_true']
        predictions <- loo.cv[loo.cv$method == method.name, 'y_pred']
        ROC <- roc(unlist(labels), unlist(predictions), quiet = TRUE)
        # Combined accuracy, error, ...
        cat("\nMethod: ", method.name, "\n",
            "acc: ", acc(unlist(labels), unlist(predictions)), "\n",
            "err: ", err(unlist(labels), unlist(predictions)), "\n",
            "fpr: ", fpr(unlist(labels), unlist(predictions)), "\n",
            "tpr: ", tpr(unlist(labels), unlist(predictions)), "\n",
            "fnr: ", fnr(unlist(labels), unlist(predictions)), "\n",
            "tnr: ", tnr(unlist(labels), unlist(predictions)), "\n",
            "auc: ", ROC$auc, "\n",
            "auc sd: ", sqrt(var(ROC)), "\n",
            sep = '')
    }

    loo.cv
}

methods <- list(
    KPIR_LS  = function(...) tensor_predictor(..., method = "KPIR_LS"),
    KPIR_MLE = function(...) tensor_predictor(..., method = "KPIR_MLE"),
    KPFC1    = function(...) tensor_predictor(..., method = "KPFC1"),
    KPFC2    = function(...) tensor_predictor(..., method = "KPFC2"),
    LSIR     = LSIR
)

#                 ppc, tpc
#              ------------
params <- list( c(  4,   3)
              , c( 15,  15)
              , c( 30,  20)
)

for (param in params) {
    c(ppc, tpc) %<-% param
    sim <- egg_analysis_reduced(methods, ppc, tpc)

    attr(sim, 'param') <- c(ppc = ppc, tpc = tpc)

    saveRDS(sim, file = sprintf('eeg_analysis_reduced_%d_%d.rds', ppc, tpc))
}
