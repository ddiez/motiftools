#' plotMatrix
#'
#' plotMatrix
#'
#' @param object a matrix object or list of matrix objects.
#' @param tree a tree object.
#' @param fill color used for tiles.
#' @param color color for tile borders.
#' @param high vector (factor or character) indicating row groups.
#' @param high.col colors usef for each group in high.
#' @param bar.percentage logical; whether to show percentages in bar plot (default: TRUE).
#' @param plot logical; whether to draw the plot (default: TRUE).
#'
#' @return Returns a grob object invisibly.
#' @export
#'
#' @examples
#' NULL
plotMatrix <- function(object, tree, fill, color = "transparent", high, high.col, bar.percentage = TRUE, plot = TRUE) {
  if (class(object) != "list")
    object <- list(object)
  
  n <- length(object) # number of heatmaps.
  nc <- ncol(object[[1]]) # number of columns of matrices.
  nr <- nrow(object[[1]]) # number of rows of matrices.
  
  ## create plots.
  # tree.
  if (missing(tree)) {
    m <- do.call(cbind, object)
    tree <- as.dendrogram(hclust(dist(m)))
  } else {
    tree <- as.dendrogram(as.hclust(tree))
  }
  g <- suppressMessages(ggdendrogram(tree, rotate = FALSE, theme_dendro = FALSE) + coord_flip() + scale_x_continuous(limits = c(.5, nr + .5), expand = c(0, 0)) + scale_y_reverse() + theme_void())
  grob_tree <- ggplotGrob(g)
  
  # reorder data.
  object <- lapply(object, function(o) o[order.dendrogram(tree), ])
  
  # highlight.
  if (!missing(high)) {
    high <- high[order.dendrogram(tree)]
    tmp <- melt(matrix(high), nrow = nr)
    if (missing(high.col)) {
      high.col <- rainbow(nlevels(tmp$value))
    }
    if (!is.null(names(high)))
      high.name <- names(high)
    else
      high.name <- ""
    g <- ggplot(tmp, aes(x = factor(Var2), y = factor(Var1), fill = value)) + geom_tile() + scale_fill_manual(high.name, values = high.col) + theme(legend.key.size = unit(.5, "lines"))  + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + guides(fill = guide_legend(direction = "horizontal"))
    grob_high <- ggplotGrob(g)
  }
  
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
    d <- d %>% group_by(motif) %>% summarize(value = sum(count), total = n())
    if (bar.percentage)
      d <- d %>% mutate(value = 100 * value / total)
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
  
  # add highlights.
  if (!missing(high)) {
    gt <- gtable_add_cols(gt, unit(1, "lines"), pos = -1)
    gg <- gtable_filter(grob_high, "panel")
    gt <- gtable_add_grob(gt, gg, t = 2, l = ncol(gt))
  }
  
  # add padding between matrices.
  gt <- gtable_add_col_space(gt, unit(.5, "lines"))
  gt <- gtable_add_row_space(gt, unit(.5, "lines"))
  
  # add barplot axis.
  gt <- gtable_add_cols(gt, unit(1.5, "lines"), pos = 0)
  gg <- gtable_filter(grob_barplot[[1]], "axis-l")
  gt <- gtable_add_grob(gt, gg, t = 1, l = 1)
  
  # add tree.
  gt <- gtable_add_cols(gt, unit(.5, "null"), pos = 0)
  gg <- gtable_filter(grob_tree, "panel")
  gt <- gtable_add_grob(gt, gg, t = 3, l = 1, r = 2)
  
  # add guides.
  # FIX: can't have this and the high legend together.
  # TODO: make sure all heatmaps have same guide range.
  #gg <- gtable_filter(grob_heatmap[[1]], "guide-box")
  #gt <- gtable_add_grob(gt, gg, t = 1, l = 1)
  
  if (!missing(high)) {
    gg <- gtable_filter(grob_high, "guide-box")
    gt <- gtable_add_grob(gt, gg, t = 1, l = 1)
  }
  
  gt <- gtable_add_padding(gt, unit(.5, "line"))
  if (plot) {
    grid.newpage()
    grid.draw(gt)
  }
  invisible(gt)
}