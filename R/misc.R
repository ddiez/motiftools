#' dcor
#' 
#' Compute Pearson correlation distance on the rows.
#'
#' @param x matrix with numeric values.
#' @param use argument passed down to cor().
#'
#' @return a dist object.
#' @export
#'
#' @examples
#' NULL
dcor <- function(x, use = "pairwise") {
  as.dist(1 - cor(t(x), use = use))
}

df2matrix = function(d, index.cols = c(1,2), value.col = 3, compact.fun = max) {
  if(anyDuplicated(d[, index.cols])) {
    sel = unique(d[, index.cols])
    nd = data.frame(sel, value = apply(sel, 1, function(r) compact.fun(d[d[,1] == r[1] & d[,2] == r[2], value.col])))
  } else {
    nd = data.frame(d[ , c(index.cols, value.col)])
  }
  colnames(nd) = c("source", "target", "value")
  as.matrix(xtabs(value ~ source + target, nd, sparse = TRUE))
}

exportTreedyn = function(x, filename) {
  d = data.frame(Id = rownames(x))
  for(m in 1:ncol(x)) {
    d <- cbind(d, colnames(x)[m], paste("{", x[,m], "}", sep = ""))
  }
  colnames(d) = c("Id", rep(1:ncol(x), each = 2))
  if(!missing(filename))
    write.table(d, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  invisible(d)
}


# find first tree with species.
findSubTree = function(tree, has.any, has.one) {
  has.any = paste(has.any, collapse = "|")
  has.one = paste(has.one, collapse = "|")
  #print(has.any)
  #print(has.one)
  s = subtrees(tree)
  sel = unlist(lapply(s, function(x) any(grepl(has.any, x$tip.label))))
  if(!(any(sel))) stop("no matches for has.any ids")
  #print(which(sel))
  s2 = s[sel]
  found = FALSE
  sel = NA
  foo = sapply(length(s2):1, function(k) {
    #print(k)
    if(any(grepl(has.one, s2[[k]]$tip.label))) {
      if(!found) {
        #print("Found")
        found <<- TRUE
        sel <<- k
      }
    }
  })
  #print(sel)
  s2[[sel]]
}
