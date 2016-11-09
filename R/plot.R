#' plotMotifMatrix
#'
#' Plot a heatmap with motif counts per sequence and associated tree and 
#' barplots with motif-wise percentages.
#'
#' @param object a matrix object or list of matrix objects.
#' @param ... arguments passed down to methods.
#' @param tree a tree object.
#' @param fill color used for tiles.
#' @param color color for tile borders.
#' @param annot list of matrices containing annotations.
#' @param annot.fill list of character vectors with fill colors for annotation matrices.
#' @param high vector (factor or character) indicating row groups.
#' @param high.col colors usef for each group in high.
#' @param bar.percentage logical; whether to show percentages in bar plot (default: TRUE).
#' @param plot logical; whether to draw the plot (default: TRUE).
#'
#' @return Returns a grob object invisibly.
#' @export
#' @rdname plotMotifMatrix-methods
#'
#' @examples
#' NULL
setGeneric("plotMotifMatrix", function(object, ...) standardGeneric("plotMotifMatrix"))

#' @rdname plotMotifMatrix-methods
#' @aliases plotMotifMatrix,list-method
setMethod("plotMotifMatrix", "list", 
function(object, tree, fill, color = "transparent", annot = NULL, annot.fill = NULL, high, high.col, bar.percentage = TRUE, plot = TRUE) {
  # check object type.
  type <- unique(sapply(object, class))
  if (length(type) > 1) stop("passing a list of objects of different class are not allowed.")
  
  if (! type %in% c("matrix", "MotifSearchResult"))
    stop("Only objects of class 'matrix' or 'MotifSearchResult' can be used with this function")
  
  if (type == "MotifSearchResult")
    object <- lapply(object, getMotifMatrix)
  
  # initialize commonly used variables.
  n <- length(object) # number of heatmaps.
  nc <- ncol(object[[1]]) # number of columns of matrices.
  nr <- nrow(object[[1]]) # number of rows of matrices.
  
  ## create plots.
  # tree.
  # if tree not provided cluster motif profiles based on Pearson distance.
  if (missing(tree)) {
    m <- do.call(cbind, object)
    tree <- as.phylo(hclust(dcor(m)))
  } else {
    if (class(tree) != "phylo")
      tree <- as.phylo(as.hclust(tree))
  }
  g <- ggtree(tree, branch.length = "none") + scale_y_continuous(limits = c(.5, nr + .5), expand = c(0, 0))
  grob_tree <- ggplotGrob(g)
  
  # reorder data.
  object <- lapply(object, function(o) o[tree$tip.label, , drop = FALSE])
  
  # annotations as matrices.
  if (!is.null(annot)) {
    na <- length(annot)
    nl <- names(annot)
    if (is.null(nl)) nl <- paste0("annot-", seq_len(na))
    # reorder annotations.
    annot <- lapply(annot, function(o) o[tree$tip.label, , drop = FALSE])
    # create plot.
    grob_annot <- lapply(seq_len(na), function(k) {
      l <- nl[k]
      d <- melt(annot[[k]], varnames = c("sequence", "variable"), value.name = l)
      d[[l]] <- factor(d[[l]])
      g <- ggplot(d, aes_string(x = "variable", y = "sequence", fill = l)) +
        geom_tile(color = color) +
        scale_fill_manual(values = annot.fill[[k]]) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0))
      ggplotGrob(g)
    })
  }
  
  # annotations as data.frames.
  # if (!is.null(annot)) {
  #   if (length(annot) == 1 && class(annot) != "list") {
  #     annot <- list(annot = annot)
  #   }
  #   #annot <- lapply(annot, function(o) o[tree$tip.label, , drop = FALSE])
  #   rob_annot <- lapply(seq_len(na), function(k) {
  #     an <- annot[[k]]
  #     an <- an[tree$tip.label, , drop = FALSE] # reorder.
  #     # for data.frame:
  #     an <- an %>% rownames_to_column("sequence")
  #     g <- ggplot(d, aes_string(x = "variable", y = "sequence", fill = l)) +
  #             geom_tile(color = color) +
  #             scale_fill_manual(values = annot.fill[[k]]) +
  #             scale_x_discrete(expand = c(0,0)) +
  #             scale_y_discrete(expand = c(0,0))
  #     ggplotGrob(g)
  #   })
  #   
  # }
  # highlight.
  # if (!missing(high)) {
  #   high <- high[tree$tip.label]
  #   tmp <- melt(matrix(high), nrow = nr)
  #   if (missing(high.col)) {
  #     high.col <- rainbow(nlevels(tmp$value))
  #   }
  #   if (!is.null(names(high)))
  #     high.name <- names(high)
  #   else
  #     high.name <- ""
  #   g <- ggplot(tmp, aes_string(x = "factor(Var2)", y = "factor(Var1)", fill = "value")) +
  #     geom_tile() + 
  #     scale_fill_manual(high.name, values = high.col) + 
  #     theme(legend.key.size = unit(.5, "lines"))  + 
  #     scale_x_discrete(expand = c(0,0)) + 
  #     scale_y_discrete(expand = c(0,0)) + 
  #     guides(fill = guide_legend(direction = "horizontal"))
  #   grob_high <- ggplotGrob(g)
  # }
  
  # heatmaps.
  if (missing(fill)) {
    fill <- cm.colors(128)
  }
  
  grob_heatmap <- lapply(seq_len(n), function(k) {
    d <- melt(object[[k]], varnames = c("sequence", "motif"), value.name = "count")
    g <- ggplot(d, aes_string(x = "motif", y = "sequence", fill = "count")) +
      geom_tile(color = color) +
      scale_fill_gradientn(colours = fill) +
      theme(legend.key.size = unit(.5, "lines")) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) 
    ggplotGrob(g)
  })
  
  # barplots.
  grob_barplot <- lapply(seq_len(n), function(k) {
    d <- melt(object[[k]], varnames = c("sequence", "motif"), value.name = "count")
    d <- d %>% group_by_("motif") %>% summarize_(value = "sum(count)", total = "n()")
    if (bar.percentage)
      d <- d %>% mutate_(value = "100 * value / total")
    g <- ggplot(d, aes_string(x = "motif", y = "value")) +
      geom_bar(stat = "identity", fill = "grey", width = 1) +
      scale_x_discrete(breaks = 1:100, expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    if (bar.percentage)
      g <- g + geom_hline(yintercept = c(0, 25, 50, 75, 100), size = .3, color = "grey20")
    ggplotGrob(g)
  })
  
  ## plot.
  # create main plot structure.
  gt <- gtable(unit(rep(1, n), "null"), unit(c(1, 5), "null"))
  
  # add heatmaps and barplots.
  for (k in seq_len(n)) {
    gg <- gtable_filter(grob_barplot[[k]], "panel")
    gt <- gtable_add_grob(gt, gg, t = 1, l = k)
    gg <- gtable_filter(grob_heatmap[[k]], "panel")
    gt <- gtable_add_grob(gt, gg, t = 2, l = k)
  }
  
  # # add highlights.
  # # if (!missing(high)) {
  # #   gt <- gtable_add_cols(gt, unit(1, "lines"), pos = -1)
  # #   gg <- gtable_filter(grob_high, "panel")
  # #   gt <- gtable_add_grob(gt, gg, t = 3, l = -1)
  # # }
  # 
  # add annotations.
  if (!is.null(annot)) {
    gt <- gtable_add_cols(gt, widths = unit(rep(1, na), "null"))
    for (k in seq_len(na)) {
      # gg <- gtable_filter(grob_annot[[k]], "guide-box")
      # gt <- gtable_add_grob(gt, gg, t = 1, l = n + k)
      gg <- gtable_filter(grob_annot[[k]], "panel")
      gt <- gtable_add_grob(gt, gg, t = 2, l = n + k)
    }
  }

  # add padding between matrices.
  gt <- gtable_add_col_space(gt, unit(.5, "lines"))
  gt <- gtable_add_row_space(gt, unit(.5, "lines"))
  
  # add barplot axis.
  gg <- gtable_filter(grob_barplot[[1]], "axis-l")
  gt <- gtable_add_cols(gt, gg$widths, pos = 0)
  gt <- gtable_add_grob(gt, gg, t = 1, l = 1)
   
  # add barplot title.
  if (!is.null(names(object))) {
    gt <- gtable_add_rows(gt, unit(1, "lines"), pos = 0)
    gg <- lapply(names(object), function(objn) {
      grid::textGrob(objn, x = unit(.5, "npc"), y = unit(.5, "npc"))
    })
    for (k in seq_len(n))
      gt <- gtable_add_grob(gt, gg[[k]], t = 1, l = 2 * k)
  }
  # 
  # # add guides.
  # # FIX: can't have this and the high legend together.
  # # TODO: make sure all heatmaps have same guide range.
  # #gg <- gtable_filter(grob_heatmap[[1]], "guide-box")
  # #gt <- gtable_add_grob(gt, gg, t = 1, l = 1)
  # 
  # if (!missing(high)) {
  #   gg <- gtable_filter(grob_high, "guide-box")
  #   gt <- gtable_add_grob(gt, gg, t = 1, l = 1)
  # }
  
  # add tree.
  gt <- gtable_add_cols(gt, unit(.5, "null"), pos = 0)
  gg <- gtable_filter(grob_tree, "panel")
  gt <- gtable_add_grob(gt, gg, t = -1, l = 1, r = 2)
  
  gt <- gtable_add_padding(gt, unit(.5, "line"))
  if (plot) {
    grid.newpage()
    grid.draw(gt)
  }
  invisible(gt)
})

#' @rdname plotMotifMatrix-methods
#' @aliases plotMotifMatrix,MotifSearchResult-method
setMethod("plotMotifMatrix", "MotifSearchResult",
function(object, ...) {
  plotMotifMatrix(list(getMotifMatrix(object)), ...)
})

#' @rdname plotMotifMatrix-methods
#' @aliases plotMotifMatrix,matrix-method
setMethod("plotMotifMatrix", "matrix",
function(object, ...) {
  plotMotifMatrix(list(object), ...)
})