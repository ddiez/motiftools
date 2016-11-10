#' motiftools
#' 
#' Tools to analyze sequence motifs
#' 
#' @name motiftools-package
#' @aliases motiftools
#' @docType package
#' 
#' @import methods Biobase Biostrings XML xml2 gtable ggplot2
#' @importFrom ape as.phylo
#' @importFrom ggtree ggtree
#' @importFrom reshape2 melt
#' @importFrom dplyr filter_ select_ mutate_ group_by_ summarize_ bind_rows %>%
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom Rcpp evalCpp
#' @importFrom grid grid.newpage grid.draw
#' @importFrom stats as.dist hclust cor
#' 
#' @useDynLib motiftools
NULL