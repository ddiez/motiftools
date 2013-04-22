plotMatrix = function(... , tree, highlight, sep = 2, adj = 10, cex = 0.5, high.col, main) {
  X = list(...)
  M = lapply(X, getMotifMatrix)
  C = lapply(X, function(x) getMotifCount(x, percentage = TRUE))
  
  if (missing(tree)) {
    # the first matrix is used to construct the tree.
    d = dist(M[[1]])
    h = hclust(d)
    require(ape)
    tree = as.phylo(h)
  }
  
  # reorder matrices according to tree.
  M = lapply(M, function(x) {
    x[tree$tip,]
  })
  
  l = layout(matrix(c(0,1,seq(2,length(M)*2+1)), 2, length(M) + 1), heights = c(1, 10), widths = c(5, rep(5, length(M)*2+1)))
  
  op = par(mar = c(0.5,0.5,0.5,0))
  plot(tree, cex = 0.6, show.tip.label = FALSE, root.edge = TRUE, use.edge.length = FALSE, yaxs = "i")
  
  if(!missing(highlight)) {
    l = tree$tip.label
    col = rep("transparent", length(l))
    if(missing(high.col))
      high.col = rainbow(length(highlight))
    if(length(high.col) == 1) {
      highlight = paste(highlight, collapse = "|")
    } else {
      for(k in 1:length(highlight)) {
        col[grepl(highlight[k], l)] = high.col[k]
      }
    }
    tiplabels(pch = 22, col = col, cex = cex, adj = adj, bg = col)
  }
  
  for(k in 1:length(M)) {
    m = M[[k]]
    c = C[[k]]
    if(missing(main))
      mn = NA
    else
      mn = main[k]
    par(mar = c(0, 0, 0.5,sep))
    barplot(c, names.arg = NA, space = 0, xaxs = "i", las = 1, ylim = c(0, 105), main = mn, , xlab = NA, ylab = NA)
    
    par(mar = c(0.5,0,0.5,sep))
    image(1:ncol(m), 1:nrow(m), t(m), col = c("white", "gray"), axes = FALSE, xlab = NA, ylab = NA)
    box()
  }
}
