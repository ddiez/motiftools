#' conservationMatrix
#' 
#' Compute the conservation score matrix from a MSA
#'
#' @param object AAMultipleAlignment alignment object.
#' 
#' @export
conservationMatrix <- function(object) {
  x <- as.matrix(object)

  m <- matrix(1L, ncol = ncol(x), nrow = nrow(x))
  rownames(m) <- rownames(x)
  colnames(m) <- apply(x, 2, function(xk) {
    tt <- table(xk)
    names(tt)[which.max(tt)]
    })

  m[x != "-"] = 2L
  
  mt <- apply(x, 2, table, exclude = "-")
  mt <- lapply(mt, function(z) 100 * z / nrow(x))
  for (k in seq_len(length(mt))) {
    tmp <- mt[[k]]
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
  g <- ggplot(d, aes_string(x="position",y="sequences",fill="conservation")) + geom_raster() + scale_fill_manual(values=c("white","grey80","grey50", "steelblue","orange"), limits = c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"), guide="legend") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand=c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
  if(! seq.names)
    g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  print(g)
  invisible(g)
}

#' plotConservationMatrix
#' 
#' plots a heatmap with the conservation score for each residue
#'
#' @param x a matrix as obtained with conservationMatrix()
#' @param tree a tree of class dendrogram or that can be coerced into.
#' @param color color for the tile borders (default: transparent).
#' 
#' @export
plotConservationMatrix <- function(x, tree, color = "transparent") {
  if (!missing(tree)) {
    tree <- as.phylo(tree)
    x <- x[tree$tip.label, ]
  }

  nr <- nrow(x)
  
  # tree grob.
  ptree <- ggtree(tree, branch.length = "none") + 
    scale_y_continuous(limits = c(.5, nr + .5), expand = c(0, 0))
  gtree <- ggplotGrob(ptree)
  
  # heatmap.
  d <- melt(x, varnames = c("sequences","position"), value.name = "conservation")
  d$conservation <- factor(d$conservation, levels = 1:5, labels =  c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"))
  
  pheat <- ggplot(d, aes_string(x = "position", y = "sequences", fill = "conservation")) + 
    geom_tile(color = color) + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c("white", "grey80", "grey50", "steelblue", "orange"),
                      limits = c("gap (-)", "< 40%", ">= 40%", ">= 60%", ">= 80%"), guide = "legend")
  gheat <- ggplotGrob(pheat)
  
  g <- gtable(widths = unit(c(1, 5, 1.5), "null"), heights = unit(c(5, 1), c("null", "lines")))
  g <- gtable_add_grob(g, gtable_filter(gtree, "panel"), t = 1, l = 1)
  g <- gtable_add_grob(g, gtable_filter(gheat, "panel"), t = 1, l = 2)
  g <- gtable_add_grob(g, gtable_filter(gheat, "guide-box"), t = 1, l =3)
  g <- gtable_add_grob(g, gtable_filter(gheat, "xlab"), t = 2, l = 2)
  grid.newpage()
  grid.draw(g)
  invisible(g)
}