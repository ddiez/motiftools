#' getMotifMatrix
#' 
#' Obtain a congintency matrix with the number of motifs per sequence. Alterna-
#' tivelly return only presence/absence of motif, and/or only motifs/sequences
#' with at least one non-zero entry.
#'
#' @param object a MotifSearchResult object.
#' @param simplify whether to return a matrix of 0 (motif not present) and 1 (motif present).
#' @param filter whether remove columns/rows with all zero entries.
#'
#' @return a matrix.
#' @export
#' @rdname getMotifMatrix-methods
#'
#' @examples
#' NULL
setGeneric("getMotifMatrix", function(object, simplify = TRUE, filter = FALSE) standardGeneric("getMotifMatrix"))

#' @rdname getMotifMatrix-methods
#' @aliases getMotifMatrix,MotifSearchResult-method
setMethod("getMotifMatrix", "MotifSearchResult", 
function(object, simplify = TRUE, filter = FALSE) {
  r <- object@ranges
  tmp <- data.frame(
    sequence = as.character(seqnames(r)),
    motif = r$motif_name,
    stringsAsFactors = FALSE
  )
  tmp$sequence <- factor(tmp$sequence, levels = sequenceNames(object))
  tmp$motif <- factor(tmp$motif, levels = motifNames(object))
  tmp <- unclass(table(tmp))
  mode(tmp) <- "integer"
  
  # simiply?
  if (simplify)
    tmp[tmp != 0] = 1L
  
  # filter?
  if (filter)
    tmp <- tmp[rowSums(tmp) > 0, colSums(tmp) > 0]
  
  tmp
})

#' getMotifCoverage
#' 
#' Obtain the percentage of each sequence that is covered by the motifs matching
#' to it. This enable us to know how much the motifs are representative of the
#' features defining each sequence.
#'
#' @param object a MotifSearchResult object.
#'
#' @return A numeric vector with the percentage of sequence covered by motifs.
#' @export
#' @rdname getMotifCoverage-methods
#'
#' @examples
#' NULL
setGeneric("getMotifCoverage", function(object) standardGeneric("getMotifCoverage"))

#' @rdname getMotifCoverage-methods
#' @aliases getMotifCoverage,MotifSearchResult-method
setMethod("getMotifCoverage", "MotifSearchResult",
function(object) {
  r <- object@ranges
  p <- pData(object@sequences)
  d <- r %>% as.data.frame
  d$tot_length <- p[d$seqnames, "length"]
  
  d %>% group_by_("seqnames") %>% 
    summarize_(length = "sum(width)", tot_length = "unique(tot_length)") %>% 
    mutate_(perc = "length / tot_length") %>% 
    select_("perc") %>% unlist(use.names = FALSE)
})


#' getMotifTotalCoverage
#' 
#' Obtain the total percentage of all sequences that is covered by the motifs matching
#' to it.
#' 
#' @param object a MotifSearchResult object.
#'
#' @return A numeric scalar with the percentage of sequence covered by motifs.
#' @export
#' @rdname getMotifTotalCoverage-methods
#'
#' @examples
#' NULL
setGeneric("getMotifTotalCoverage", function(object) standardGeneric("getMotifTotalCoverage"))

#' @rdname getMotifTotalCoverage-methods
#' @aliases getMotifTotalCoverage,MotifSearchResult-method
setMethod("getMotifTotalCoverage", "MotifSearchResult",
function(object) {
  r <- reduce(object@ranges)
  pdata <- pData(object@sequences)
  100 * sum(width(r)) / sum(pdata[, "length"])
})
