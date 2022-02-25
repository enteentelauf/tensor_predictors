#' Angle between two subspaces
#'
#' Computes the principal angle between two subspaces spaned by the columns of
#' the matrices \code{A} and \code{B}.
#'
#' @param A,B Numeric matrices with column considered as the subspace spanning
#'  vectors. Both must have the same number of rows (a.k.a must live in the
#'  same space).
#' @param is.orth boolean determining if passed matrices A, B are allready
#'  orthogonalized. If set to TRUE, A and B are assumed to have orthogonal
#'  columns (which is not checked).
#'
#' @returns angle in radiants.
#'
subspace <- function(A, B, is.orth = FALSE) {
    if (!is.numeric(A) || !is.numeric(B)) {
        stop("Arguments 'A' and 'B' must be numeric.")
    }
    if (is.vector(A)) A <- as.matrix(A)
    if (is.vector(B)) B <- as.matrix(B)
    if (nrow(A) != nrow(B)) {
        stop("Matrices 'A' and 'B' must have the same number of rows.")
    }
    if (!is.orth) {
        A <- qr.Q(qr(A))
        B <- qr.Q(qr(B))
    }
    if (ncol(A) < ncol(B)) {
        tmp <- A; A <- B; B <- tmp
    }
    for (k in 1:ncol(A)) {
        B <- B - tcrossprod(A[, k]) %*% B
    }
    asin(min(1, La.svd(B, 0L, 0L)$d))
}
