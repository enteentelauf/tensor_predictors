#' Approximates kronecker product decomposition.
#'
#' Approximates the matrices `A` and `B` such that
#'      C = A %x% B
#' with `%x%` the kronecker product of the matrixes `A` and `B`
#' of dimensions `dimA` and `dimB` respectively.
#'
#' @param C desired kronecker product result.
#' @param dimA length 2 vector of dimensions of \code{A}.
#' @param dimB length 2 vector of dimensions of \code{B}.
#'
#' @return list with attributes `A` and `B`.
#'
#' @examples
#' A <- matrix(seq(14), 7, 2)
#' B <- matrix(c(T, F), 3, 4)
#' C <- kronecker(A, B) # the same as 'C <- A %x% B'
#' approx.kronecker(C, dim(A), dim(B))
#'
#' @seealso C.F. Van Loan / Journal of Computational and Applied Mathematics
#'          123 (2000) 85-100 (pp. 93-95)
#'
#' @imports RSpectra
#'
approx.kronecker <- function(C, dimA, dimB) {

    dim(C) <- c(dimB[1L], dimA[1L], dimB[2L], dimA[2L])
    R <- aperm(C, c(2L, 4L, 1L, 3L))
    dim(R) <- c(prod(dimA), prod(dimB))

    svdR <- try(RSpectra::svds(R, 1L), silent = TRUE)
    if (is(svdR, 'try-error')) {
        svdR <- svd(R, 1L, 1L)
    }

    return(list(
        A = array(sqrt(svdR$d[1]) * svdR$u, dimA),
        B = array(sqrt(svdR$d[1]) * svdR$v, dimB)
    ))
}
