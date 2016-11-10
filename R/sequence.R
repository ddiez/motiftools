#' conservationMatrix
#' 
#' Compute the conservation score matrix from a MSA
#'
#' @param object AAMultipleAlignment alignment object.
#' 
#' @export
#' @rdname conservationMatrix-methods
#'
#' @examples
#' NULL
setGeneric("conservationMatrix", function(object) standardGeneric("conservationMatrix"))

#' @rdname conservationMatrix-methods
#' @aliases conservationMatrix,matrix-method
setMethod("conservationMatrix", "matrix", 
function(object) {
  m <- matrix(1L, ncol = ncol(object), nrow = nrow(object))
  rownames(m) <- rownames(object)
  colnames(m) <- apply(object, 2, function(x) {
    tt <- table(x)
    names(tt)[which.max(tt)]
    })

  m[object != "-"] = 2L
  
  mt <- apply(object, 2, table, exclude = "-")
  mt <- lapply(mt, function(z) 100 * z / nrow(object))
  for (k in seq_len(length(mt))) {
    tmp <- mt[[k]]
    m[,k][object[,k] %in% names(tmp)[tmp >= 40]] = 3L
    m[,k][object[,k] %in% names(tmp)[tmp >= 60]] = 4L
    m[,k][object[,k] %in% names(tmp)[tmp >= 80]] = 5L
  }
  m
})

#' @rdname conservationMatrix-methods
#' @aliases conservationMatrix,AAMultipleAlignment-method
setMethod("conservationMatrix", "AAMultipleAlignment", 
function(object) {
  conservationMatrix(as.matrix(object))
})

#' plotConservationMatrix
#' 
#' plots a heatmap with the conservation score for each residue
#'
#' @param x a matrix as obtained with conservationMatrix()
#' @param tree a tree of class dendrogram or that can be coerced into.
#' @param color color for the tile borders (default: transparent).
#' @param show.tips logical; whether to show sequence ids.
#' @param show.consensus logical; whether to show consensus sequence.
#' @param size text size for consensus sequence.
#' 
#' @export
#' @rdname plotConservationMatrix-methods
#'
#' @examples
#' NULL
setGeneric("plotConservationMatrix", function(x, ...) standardGeneric("plotConservationMatrix"))

#' @rdname plotConservationMatrix-methods
#' @aliases plotConservationMatrix,matrix-method
setMethod("plotConservationMatrix", "matrix", 
function(x, tree, color = "transparent", show.tips = FALSE, show.consensus = FALSE, size = 8) {
  if (!missing(tree)) {
    tree <- as.phylo(tree)
    x <- x[tree$tip.label, ]
  }
  
  nr <- nrow(x)
  labels <- NULL
  if (show.consensus)
    labels <- colnames(x)
  tips <- NULL
  if (show.tips)
    tips <- rownames(x)
  colnames(x) <- seq_len(ncol(x))
  
  # tree grob.
  ptree <- ggtree(tree, branch.length = "none") + 
    scale_y_continuous(limits = c(.5, nr + .5), expand = c(0, 0))
  gtree <- ggplotGrob(ptree)
  
  # heatmap.
  d <- melt(x, varnames = c("sequences","position"), value.name = "conservation")
  d$sequences <- factor(d$sequences)
  d$position <- factor(d$position)
  d$conservation <- factor(d$conservation, levels = 1:5, labels =  c("gap (-)","< 40%",">= 40%",">= 60%",">= 80%"))
  
  pheat <- ggplot(d, aes_string(x = "position", y = "sequences", fill = "conservation")) + 
    geom_tile(color = color) + 
    scale_x_discrete(expand = c(0, 0), labels = labels) +
    scale_y_discrete(expand = c(0, 0), labels = tips) +
    scale_fill_manual(values = c("white", "grey80", "grey50", "steelblue", "orange"),
                      limits = c("gap (-)", "< 40%", ">= 40%", ">= 60%", ">= 80%"), guide = "legend") + 
    theme(axis.ticks = element_blank(), axis.text.x = element_text(size = size))
  gheat <- ggplotGrob(pheat)
  
  gt <- gtable(widths = unit(c(.5, 5), "null"), heights = unit(5, "null"))
  gt <- gtable_add_grob(gt, gtable_filter(gtree, "panel"), t = 1, l = 1)
  gt <- gtable_add_grob(gt, gtable_filter(gheat, "panel"), t = 1, l = 2)
  
  if (show.tips) {
    gg <- gtable_filter(gheat, "axis-l")
    gt <- gtable_add_cols(gt, widths = gg$widths, pos = 1)
    gt <- gtable_add_grob(gt, gg, t = 1, l = 2)
  }
  
  if (show.consensus) {
    gg <- gtable_filter(gheat, "axis-b")
    gt <- gtable_add_rows(gt, heights = gg$heights, pos = 0)
    gt <- gtable_add_row_space(gt, unit(0.5, "lines"))
    gt <- gtable_add_grob(gt, gg, t = 1, l = -1)
  }
  
  gt <- gtable_add_rows(gt, heights = unit(1, "lines"), pos = -1)
  gt <- gtable_add_grob(gt, gtable_filter(gheat, "xlab"), t = -1, l = -1)
  
  gg <- gtable_filter(gheat, "guide-box")
  gt <- gtable_add_cols(gt, widths = gg$widths, pos = -1)
  gt <- gtable_add_grob(gt, gg, t = -2, l = -1)
  
  gt <- gtable_add_padding(gt, unit(.5, "lines"))
  grid.newpage()
  grid.draw(gt)
  invisible(gt)
})

#' @rdname plotConservationMatrix-methods
#' @aliases plotConservationMatrix,AAMultipleAlignment-method
setMethod("plotConservationMatrix", "AAMultipleAlignment", 
function(x, tree, color = "transparent", show.tips = FALSE, show.consensus = FALSE, size = 8) {
  plotConservationMatrix(conservationMatrix(x), tree = tree, color = color, show.tips = show.tips, show.consensus = show.consensus, size = size)
})