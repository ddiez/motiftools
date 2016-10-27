#' @title conservationMatrix
#' @description computes the conservation score matrix from a MSA
#' @param object sequence alignment object.
conservationMatrix = function(object) {
  x = as.matrix(object)
  rownames(x)=names(object)
  
  m = matrix(1L, ncol=ncol(x),nrow=nrow(x))
  rownames(m)=rownames(x)
  #colnames(m)=1:ncol(m)
  m[x != "-"]=2L
  
  mt=apply(x,2,table,exclude="-")
  mt=lapply(mt,function(z) 100*z/nrow(x))
  for(k in 1:length(mt)) {
    tmp=mt[[k]]
    m[,k][x[,k] %in% names(tmp)[tmp >= 40]] = 3L
    m[,k][x[,k] %in% names(tmp)[tmp >= 60]] = 4L
    m[,k][x[,k] %in% names(tmp)[tmp >= 80]] = 5L
  }
  m
}

#' plotConservation
#'
#'plotConservation
#'
#' @param x a matrix as obtained with conservationMatrix()
#' @param seq.names logical; whether to plot seq names
#' @param cluster.row logical whether to cluster rows
#' 
#' @note This function will be superseeded by plotConservationMatrix()
#'
#' @return NULL
#' @export
#'
plotConservation = function(x, seq.names = FALSE, cluster.row = TRUE) {
  if(cluster.row) {
    h <- hclust(dist(x))
    x <- x[h$order,]
  }
  d = melt(x,varnames = c("sequences","position"), value.name = "conservation")
  d$conservation = factor(d$conservation, levels = 1:5, labels =  c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"))
  g <- ggplot(d, aes(x=position,y=sequences,fill=conservation)) + geom_raster() + scale_fill_manual(values=c("white","grey80","grey50", "steelblue","orange"), limits = c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"), guide="legend") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand=c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
  if(! seq.names)
    g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  print(g)
  invisible(g)
}

#' @title plotConservationMatrix
#' @description plots a heatmap with the conservation score for each residue
#' @param x a matrix as obtained with conservationMatrix()
#' @param tree a tree of class phylo (package ape)
plotConservationMatrix <- function(x, tree) {
  if(!missing(tree))
    x <- x[tree$tip, ]
  
  # tree grob.
  ptree <- ggtree(tree, branch.length = "none", ladderize = FALSE)
  gtree <- ggplotGrob(ptree + theme(plot.margin = unit(c(1,0,1,1), "lines")) + scale_y_discrete())
  
  # heatmap.
  d <- melt(x, varnames = c("sequences","position"), value.name = "conservation")
  d$conservation <- factor(d$conservation, levels = 1:5, labels =  c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"))
  
  gheat <- ggplot(d, aes(x = position, y = sequences, fill = conservation)) +
    geom_tile() + ylab(NULL) + scale_x_discrete() + scale_y_discrete() + theme(panel.margin = unit(c(0,0,0,0), "null"), axis.ticks = element_blank(), axis.text = element_blank(), plot.margin = unit(c(1,1,1,0), "lines")) + scale_fill_manual(values=c("white","grey80","grey50", "steelblue","orange"), limits = c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"), guide="legend")
  
  g <- gtable_add_cols(gtree, unit(5, "null"))
  g <- gtable_add_grob(g, ggplotGrob(gheat), t = 1, b = nrow(g), l = ncol(g))
  grid.newpage()
  grid.draw(g)
  invisible(g)
}