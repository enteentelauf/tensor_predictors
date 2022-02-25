# # Generate Sample Data.
# n <- 250
# # see: simulation_binary.R
# data <- simulateData.binary(n / 2, n / 2, (p <- 10), (t <- 5), 0.3, 0.3)
# X <- data$X
# colnames(X) <- paste('X[', outer(1:p, 1:t, paste, sep = ','), ']', sep = '')
# Y <- 2 * data$Y
# write.csv(data.frame(X, Y), file = 'example_data.csv', row.names = FALSE)

suppressPackageStartupMessages({
    library(pROC)
})

source('../tensor_predictors/tensor_predictors.R')

# Read sample data from file and split into predictors and responces.
data <- read.csv('example_data.csv')
X <- as.matrix(data[, names(data) != 'Y'])
Y <- as.matrix(data[, 'Y'])

# Set parameters (and check)
n <- nrow(X)
p <- 10
t <- 5
stopifnot(p * t == ncol(X))

# Setup folds (folds contains indices of the test set).
nr.folds <- n # leave-one-out when number of folds equals the sample size `n`.
folds <- split(sample.int(n), (seq(0, n - 1) * nr.folds) %/% n)
labels      <- vector('list', nr.folds) # True test values (per fold)
predictions <- vector('list', nr.folds) # Predictions on test set.

for (i in seq_along(folds)) {
    fold <- folds[[i]]
    # Split data into train and test sets.
    X.train <- X[-fold, ]
    Y.train <- Y[-fold, , drop = FALSE]
    X.test <- X[fold, ]
    Y.test <- Y[fold, , drop = FALSE]

    # Compute reduction (method = c('KPIR_LS' ,'KPIR_MLE', 'KPFC1', 'KPFC2', 'KPFC3'))
    # or LSIR(X.train, Y.train, p, t) in 'lsir.R'.
    dr <- tensor_predictor(X.train, Y.train, p, t, method = 'KPIR_LS')
    B <- kronecker(dr$alpha, dr$beta) # Also available: Gamma_1, Gamma_2, Gamma, B.
    # Predict via a logit model building on the reduced data.
    model <- glm(y ~ x, family = binomial(link = "logit"),
                 data = data.frame(x = X.train %*% B, y = as.integer(Y.train > 0)))

    labels[[i]] <- as.integer(Y.test > 0)
    predictions[[i]] <- predict(model, data.frame(x = X.test %*% B), type = "response")
}

# Compute classic ROC for predicted samples (mean AUC makes no sense for leave-one-out)
y.true <- unlist(labels)
y.pred <- unlist(predictions)
roc(y.true, y.pred)
