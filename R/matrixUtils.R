
#' Compute column means of sparse matrix (dgC format) of nonzero elements
#'
#' @param dgCMat matrix of class dgCMatrix
#'
#' @return vector of column means. If there is no nonzero element in a column, will return NA
#' @import Matrix
#' @export
colMeans_drop0 <- function (dgCMat) {
  nnz_per_col <- diff(dgCMat@p)
  #nnz_per_col[nnz_per_col == 0] <- 1 # Uncomment to return 0 in case of no nonzero element
  return(Matrix::colSums(dgCMat) / nnz_per_col)
}

#' Compute row means of sparse matrix (dgC format) of nonzero elements
#'
#' @param dgCMat matrix of class dgCMatrix
#'
#' @return vector of row means. If there is no nonzero element in a row, will return NA#'
#' @import Matrix
#' @export
rowMeans_drop0 <- function (dgCMat) {
  RowInd <- dgCMat@i + 1
  nnz_per_row <- tabulate(RowInd, nbins = NROW(dgCMat))
  #nnz_per_row[nnz_per_row == 0] <- 1 # Uncomment to return 0 in case of no nonzero element
  return(Matrix::rowSums(dgCMat) / nnz_per_row)
}
