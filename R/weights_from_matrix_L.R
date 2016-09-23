


#' Given a matrix L, determine the (unnnormalized) weights to be given to each individual's genotype
#'
#' @param L a matrix, such as one that comes out of matrix_L_from_pedigree.  It should
#' have row and column names and be invertible
#' @export
#' @examples
#' example(matrix_L_from_pedigree)
#' wts <- weights_from_matrix_L(Lmat)
weights_from_matrix_L <- function(Lmat) {
  colSums(solve(Lmat))
}
