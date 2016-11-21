.datacache <- new.env(hash = TRUE, parent = emptyenv())

getScoreMatrix <- function(score.matrix) {
  data(list = score.matrix, package = "Biostrings", envir = .datacache)
  get(score.matrix, envir = .datacache)
}

#' dcor
#' 
#' Compute Pearson correlation distance on the rows.
#'
#' @param x matrix with numeric values.
#' @param use argument passed down to cor().
#' @param method method use for computing correlation (default: pearson).
#'
#' @return a dist object.
#' @export
#'
#' @examples
#' NULL
dcor <- function(x, use = "pairwise", method = "pearson") {
  as.dist(1 - cor(t(x), use = use, method = method))
}

#' Euclidean distance
#' 
#' Compute euclidean distance of two matrices or vectors.
#' 
#' @param x numeric matrix or vector.
#' @param y numeric matrix or vector.
#' 
#' @return Numeric matrix or vector.
#' @export
eud <- function(x, y) {
  if (any(dim(x) != dim(y)))
    stop("x and y have to have same dimensions.")
  sqrt(sum((x - y) ^ 2))
}
