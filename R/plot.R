plotMotifMatrix = function(object, ..., tree, annot, annot.col) {
  extra_panels=0
  if(!missing(annot)) {
    if(is.matrix(annot))
      extra_panels=ncol(annot)
    else
      stop("annotations must be provided as a matrix.")
  }
  
  X = list(object, ...)
  M = lapply(X,getMotifMatrix)
  
  if (missing(tree)) {
    # the first matrix is used to construct the tree.
    d = dist(M[[1]])
    h = hclust(d)
    tree=as.phylo(h)
  }
  
  # compute layout matrix:
  lm=matrix(c(0,1),ncol=1)
  if(extra_panels>0)
    lm=cbind(lm,rbind(rep(0,extra_panels),seq(2,extra_panels+1)))
  lm=cbind(lm, matrix(seq(2+extra_panels,length(M)*2+extra_panels+1),nrow=2))
  
  
  op = par(mar = c(1,0.5,1,0))
  l = layout(lm, heights = c(1, 10), widths = c(5, rep(.5, extra_panels), rep(5, length(M))))
  # plot tree.
  plot(tree, cex = 0.6, show.tip.label = FALSE, root.edge = TRUE, use.edge.length = FALSE, yaxs = "i")
  
  # reorder matrices according to tree (needs to be done *after* the tree is plotted: conversation with Emmanuel Paradis)
  # fixed by Diego Sep 24/2013
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip <- 1:lastPP$Ntip # or seq_len(lastPP$Ntip)
  YY <- lastPP$yy[tip]
  #o <- order(YY) # This get indexes that cannot reorder the matrix. [DD]
  o = tree$tip.label[YY] # better reorder based on labels. [DD]
  M = lapply(M, function(m) m[o,])
  if(!missing(annot)) annot=annot[o,,drop=FALSE]

  # plot annotations.
  if(!missing(annot)) {
    if(missing(annot.col)) 
      annot.col=rainbow(max(annot))
    for(k in 1:ncol(annot)) {
      par(mar = c(1, 0, 1,0))
      image(t(annot[,k,drop=FALSE]),axes=FALSE,col=annot.col)
      box()
    }
  }
  
  # plot matrix panels.
  for(k in 1:length(M)) {
    m=M[[k]]
    par(mar = c(0, 1, 1,1))
    barplot(100*apply(m,2,sum)/nseq(object),las=1,xaxs="i",ylim=c(0,100),names=FALSE)
    title(X[[k]]@info$tool,line=0)
    par(mar = c(1, 1, 1,1))
    image(t(m),axes=FALSE,col=c("white","grey"))
    box()
  }
  par(op)
}
setGeneric("plotMotifMatrix")

plotCounts = function(x, cut) {
  plot(x$Counts, ylim = c(0, 100), xlab = "Architectures", ylab = "Counts/Percentage", axes = FALSE, type = "l")
  points(x$Counts, col = "black", pch = 21, bg = "gray")
  lines(x$Percentage, col = "darkblue")
  lines(x$Cumulative, col = "darkred")
  points(x$Percentage, col = "darkblue", pch = 21, bg = "steelblue")
  points(x$Cumulative, col = "darkred", pch = 21, bg = "orange")
  if(!missing(cut))
    abline(v = cut, lty = "dotted")
  axis(2, las = 1)
  legend("right", c("Counts", "Percentage", "Cumulative"), pch = 21, col = c("black", "darkblue", "darkred"), pt.bg = c("gray", "steelblue", "orange"), bty = "n")
  box()
}

plotMotifCount = function(M, percentage = FALSE, ...) {
  x = getMotifCount(M, percentage)
  mp = barplot(x, las = 1, axes = FALSE, axisnames = FALSE, ...)
  axis(2, las = 1)
  mtext(1:20, side = 1, at = mp, cex = 0.8)
  box()
}

plotMatrix = function(... , tree, highlight, sep = 2, adj = 10, cex = 0.5, high.col, main) {
  X = list(...)
  #M = lapply(X, getMotifMatrix)
  M = X
  M = lapply(M, function(m) m[order(rownames(m)),])
  
  #C = lapply(X, function(x) getMotifCount(x, percentage = TRUE))
  C = lapply(M, function(x) 100*colSums(x)/nrow(x))
  
  if (missing(tree)) {
    # the first matrix is used to construct the tree.
    d = dist(M[[1]])
    h = hclust(d)
    require(ape)
    tree = as.phylo(h)
  }
  
  l = layout(matrix(c(0,1,seq(2,length(M)*2+1)), 2, length(M) + 1), heights = c(1, 10), widths = c(5, rep(5, length(M)*2+1)))
  
  op = par(mar = c(0.5,0.5,0.5,0))
  plot(tree, cex = 0.6, show.tip.label = FALSE, root.edge = TRUE, use.edge.length = FALSE, yaxs = "i")
  
  # reorder matrices according to tree (needs to be done after the tree is plotted: conversation with Emmanuel Paradis)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip <- 1:lastPP$Ntip # or seq_len(lastPP$Ntip)
  YY <- lastPP$yy[tip]
  o <- order(YY)
  M = lapply(M, function(x) {
    x[o,]
  })
  
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
    legend("topleft", legend = highlight, fill = high.col)
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
