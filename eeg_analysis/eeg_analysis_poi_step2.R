suppressPackageStartupMessages({
    library(pROC)
})

source('../tensor_predictors/approx_kronecker.R')
source('../tensor_predictors/multi_assign.R')
# Load EEG dataset
dataset <- readRDS('eeg_data.rds')
# Load EEG k-fold simulation results.
folds <- readRDS('eeg_analysis_poi.rds')
# Set dimenional parameters.
p <- 64L    # nr. of predictors (count of sensorce)
t <- 256L   # nr. of time points (measurements)

labels      <- vector('list', length(folds))
predictions <- vector('list', length(folds))
alphas <- matrix(0, length(folds), t)
betas  <- matrix(0, length(folds), p)
# For each fold.
for (i in seq_along(folds)) {
    fold <- folds[[i]]
    # Factorize POI result in alpha, beta.
    c(alpha, beta) %<-% approx.kronecker(fold$Q, c(t, 1), c(p, 1))
    # Drop small values of alpha, beta.
    alpha[abs(alpha) < 1e-6] <- 0
    beta[abs(beta) < 1e-6] <- 0
    # Reconstruct B from factorization.
    B <- kronecker(alpha, beta)
    # Select folds train/test sets.
    X_train <- as.matrix(dataset[-fold$index, -(1:2)])
    y_train <- as.factor(dataset[-fold$index, 'Case_Control'])
    X_test <- as.matrix(dataset[fold$index, -(1:2)])
    y_test <- as.factor(dataset[fold$index, 'Case_Control'])
    # Predict via a logit model building on the reduced data.
    model <- glm(y ~ x, family = binomial(link = "logit"),
                 data = data.frame(x = X_train %*% B, y = y_train))
    y_hat <- predict(model, data.frame(x = X_test %*% B), type = "response")
    # Set target and prediction values for the ROC curve.
    labels[[i]]      <- y_test
    predictions[[i]] <- y_hat
    alphas[i, ] <- as.vector(alpha)
    betas[i, ]  <- as.vector(beta)
}

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

# Combined accuracy, error, ...
cat("acc: ", acc(unlist(labels), unlist(predictions)), "\n",
    "err: ", err(unlist(labels), unlist(predictions)), "\n",
    "fpr: ", fpr(unlist(labels), unlist(predictions)), "\n",
    "tpr: ", tpr(unlist(labels), unlist(predictions)), "\n",
    "fnr: ", fnr(unlist(labels), unlist(predictions)), "\n",
    "tnr: ", tnr(unlist(labels), unlist(predictions)), "\n",
    "auc: ", roc(unlist(labels), unlist(predictions), quiet = TRUE)$auc, "\n",
    sep = '')
# Confidence interval for AUC.
ci(roc(unlist(labels), unlist(predictions), quiet = TRUE))

# Means of per fold accuracy, error, ...
cat("acc: ", mean(mapply(acc, labels, predictions)), "\n",
    "err: ", mean(mapply(err, labels, predictions)), "\n",
    "fpr: ", mean(mapply(fpr, labels, predictions)), "\n",
    "tpr: ", mean(mapply(tpr, labels, predictions)), "\n",
    "fnr: ", mean(mapply(fnr, labels, predictions)), "\n",
    "tnr: ", mean(mapply(tnr, labels, predictions)), "\n",
    "auc: ", mean(mapply(function(...) roc(...)$auc, labels, predictions,
        MoreArgs = list(direction = '<', quiet = TRUE))), "\n",
    sep = '')
# Means of per fold CI.
rowMeans(mapply(function(...) ci(roc(...)), labels, predictions,
    MoreArgs = list(direction = '<', quiet = TRUE)))
sd(mapply(function(...) roc(...)$auc, labels, predictions,
        MoreArgs = list(direction = '<', quiet = TRUE)))

################################################################################
###                                   plot                                   ###
################################################################################
multiplot <- function(..., plotlist = NULL, cols) {
    library(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # Make the panel
    plotCols = cols
    # Number of rows needed, calculated from cols
    plotRows = ceiling(numPlots / plotCols)
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y) {
        viewport(layout.pos.row = x, layout.pos.col = y)
    }
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i / plotCols)
        curCol = (i - 1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol))
    }
}

pa <- ggplot(data.frame(time = rep(1:ncol(alphas), 2),
                        means = c(colMeans(abs(alphas)), .5 * colMeans(!alphas)),
                        type = factor(rep(c(0, 1), each = ncol(alphas)),
                                      labels = c('mean', 'dropped'))),
             aes(x = time, y = means, fill = type)) +
    geom_col(position = 'dodge') +
    labs(title = 'Components of alpha', x = 'time', y = 'means') +
    coord_cartesian(ylim = c(0, 0.5)) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . * 2,
                                           name = 'dropped',
                                           labels = scales::percent)) +
    theme(legend.position = 'top',
          legend.title = element_blank())

pb <- ggplot(data.frame(time = rep(1:ncol(betas), 2),
                        means = c(colMeans(abs(betas)), .5 * colMeans(!betas)),
                        type = factor(rep(c(0, 1), each = ncol(betas)),
                                      labels = c('mean', 'dropped'))),
             aes(x = time, y = means, fill = type)) +
    geom_col(position = 'dodge') +
    labs(title = 'Components of beta', x = 'sensors', y = 'means') +
    coord_cartesian(ylim = c(0, 0.5)) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . * 2,
                                           name = 'dropped',
                                           labels = scales::percent)) +
    theme(legend.position = 'top',
          legend.title = element_blank())

multiplot(pa, pb, cols = 1)
