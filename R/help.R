#' motiftools
#' 
#' Tools to analyze sequence motifs
#' 
#' @name motiftools-package
#' @aliases motiftools
#' @docType package
#' 
#' @import methods Biobase Biostrings xml2 gtable ggplot2
#' @importFrom igraph graph_from_data_frame V as_adj
#' @importFrom ape as.phylo
#' @importFrom ggtree ggtree
#' @importFrom reshape2 melt
#' @importFrom dplyr filter_ select_ mutate_ group_by_ summarize_ bind_rows %>%
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges seqnames reduce
#' @importFrom Rcpp evalCpp
#' @importFrom grid grid.newpage grid.draw
#' @importFrom stats dist as.dist hclust as.hclust cor
#' @importFrom utils data
#' @importFrom grDevices cm.colors rainbow
#' 
#' @useDynLib motiftools
NULL