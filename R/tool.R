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

dcor <- function(x, use = "pairwise") {
  as.dist(1 - cor(t(x), use = use))
}
