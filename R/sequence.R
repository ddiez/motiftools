#' Compute the conservation score matrix from a MSA
#' 
#' Compute the conservation score matrix from a multiple sequence alignment (MSA) 
#' according to the method described in the software Jalview (htt://www.jalview.org).
#'
#' @param object AAMultipleAlignment alignment object (or a matrix).
#' 
#' @export
#' @rdname conservationMatrix-methods
#'
#' @examples
#' NULL
setGeneric("conservationMatrix", function(object) standardGeneric("conservationMatrix"))

#' @rdname conservationMatrix-methods
#' @aliases conservationMatrix,AAMultipleAlignment-method
setMethod("conservationMatrix", "AAMultipleAlignment", 
function(object) {
  m <- as.matrix(object)
  nr <- nrow(m)
  nc <- ncol(m)
  
  cons <- consensusMatrix(object) / nr # faster on AAMultipleAlignment object.
  cons <- matrix(cut(cons, breaks = c(0, .4, .6, .8, 1), include.lowest = TRUE, labels = c(2L, 3L, 4L, 5L)), ncol = ncol(cons), nrow = nrow(cons), dimnames = dimnames(cons))
  cons["-",][] <- 1L # replace gaps with "1"
  
  res <- lapply(seq_len(nc), function(j) {
    cons[, j, drop = FALSE][m[, j],] %>% unname
  })
  res <- do.call(cbind, res)
  rownames(res) <- rownames(m)
  mode(res) <- "integer"
  res
})

#' Plots a heatmap with the conservation score for each residue
#' 
#' Plots a heatmap with the conservation score for each residue colored according
#' to the scheme shown in Jalview (htt://www.jalview.org).
#'
#' @param x a matrix as obtained with conservationMatrix() or an AAMultipleAlignment object.
#' @param ... arguments passed down to matrix-method.
#' @param tree a tree of class dendrogram or that can be coerced into.
#' @param raster logical; whether to use geom_raster (best for large alignments).
#' @param tile.color color for the tile borders (default: transparent).
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
function(x, tree, raster = FALSE, tile.color = "transparent", show.tips = FALSE, show.consensus = FALSE, size = 8) {
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
  
  if (raster)
    p <- geom_raster()
  else
    p <- geom_tile(color = tile.color)
  
  pheat <- ggplot(d, aes_string(x = "position", y = "sequences", fill = "conservation")) + 
    p + 
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
function(x, tree, raster = FALSE, tile.color = "transparent", show.tips = FALSE, show.consensus = FALSE, size = 8) {
  plotConservationMatrix(conservationMatrix(x), tree = tree, raster = raster, tile.color = tile.color, show.tips = show.tips, show.consensus = show.consensus, size = size)
})