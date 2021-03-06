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
#' @param bar.percentage logical; whether to show percentages in bar plot (default: TRUE).
#' @param show.tips logical; whether to show the sequence names near the tree tips.
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
function(object, tree, fill, color = "transparent", annot = NULL, annot.fill = NULL, bar.percentage = TRUE, show.tips = FALSE) {
  # check object type.
  type <- unique(lapply(object, class))
  if (length(type) > 1) stop("passing a list of objects of different class is not allowed.") 
  
  type <- unlist(type)

  if (! any(c("matrix", "MotifSearchResult") %in% type))
    stop("Only objects of class 'matrix' or 'MotifSearchResult' can be used with this function")
  
  if ("MotifSearchResult" %in% type)
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
    tree <- as.phylo(hclust(dist(m, method = "binary")))
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
    if (is.null(nl)) nl <- paste0("annot", seq_len(na))
    # reorder annotations.
    annot <- lapply(annot, function(o) o[tree$tip.label, , drop = FALSE])
    # create plot.
    grob_annot <- lapply(seq_len(na), function(k) {
      l <- nl[k]
      d <- melt(annot[[k]], varnames = c("sequence", "variable"), value.name = l)
      d[[l]] <- factor(d[[l]], levels = seq_len(length(annot.fill[[k]])))
      g <- ggplot(d, aes(x = .data$variable, y = .data$sequence, fill = !!l)) +
        geom_tile(color = color) +
        scale_fill_manual(values = annot.fill[[k]]) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0))
      ggplotGrob(g)
    })
  }
  
  # heatmaps.
  if (missing(fill)) {
    fill <- c("white", "grey")
  }
  
  grob_heatmap <- lapply(seq_len(n), function(k) {
    d <- melt(object[[k]], varnames = c("sequence", "motif"), value.name = "count")
    g <- ggplot(d, aes(x = .data$motif, y = .data$sequence, fill = .data$count)) +
      geom_tile(color = color) +
      scale_fill_gradientn(colours = I(fill)) +
      theme(legend.key.size = unit(.5, "lines")) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_discrete(expand = c(0,0)) +
      theme(axis.ticks = element_blank())
    ggplotGrob(g)
  })
  
  # barplots.
  grob_barplot <- lapply(seq_len(n), function(k) {
    d <- melt(object[[k]], varnames = c("sequence", "motif"), value.name = "count")
    d <- d %>% group_by_at("motif") %>% summarize(value = sum(.data$count), total = n())
    if (bar.percentage)
      d <- d %>% mutate(value = 100 * .data$value / .data$total)
    g <- ggplot(d, aes(x = .data$motif, y = .data$value)) +
      geom_bar(stat = "identity", fill = "grey", width = 1) +
      scale_x_discrete(breaks = 1:100, expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    if (bar.percentage)
      g <- g + geom_hline(yintercept = c(0, 25, 50, 75, 100), size = .3, color = "grey20")
    ggplotGrob(g)
  })
  
  ## plot.
  # create main plot structure.
  gt <- gtable(unit(rep(5, n), "null"), unit(c(1, 5), "null"))
  
  # add heatmaps and barplots.
  for (k in seq_len(n)) {
    gg <- gtable_filter(grob_barplot[[k]], "panel")
    gt <- gtable_add_grob(gt, gg, t = 1, l = k)
    gg <- gtable_filter(grob_heatmap[[k]], "panel")
    gt <- gtable_add_grob(gt, gg, t = 2, l = k)
  }
  
  # add annotations.
  if (!is.null(annot)) {
    for (k in seq_len(na)) {
      gt <- gtable_add_cols(gt, widths = unit(ncol(annot[[k]]), "lines"))
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
  
  if (show.tips) {
    gg <- gtable_filter(grob_heatmap[[1]], "axis-l")
    gt <- gtable_add_cols(gt, widths = gg$widths, pos = 0)
    gt <- gtable_add_grob(gt, gg, t = -1, l = 1, r = 2)
  }
  
  # add tree.
  gt <- gtable_add_cols(gt, unit(1, "null"), pos = 0)
  gg <- gtable_filter(grob_tree, "panel")
  if (show.tips)
    gt <- gtable_add_grob(gt, gg, t = -1, l = 1)
  else
    gt <- gtable_add_grob(gt, gg, t = -1, l = 1, r = 2)
  
  
  gt <- gtable_add_padding(gt, unit(.5, "line"))
  grid.newpage()
  grid.draw(gt)
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