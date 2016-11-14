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
  s <- object@sequences
  
  p <- Biobase::pData(s)
  r %>% as.data.frame %>% 
    group_by_("seqnames") %>% 
    summarize_(length = "sum(width)", tot_length = p[unique(seqnames), "length"]) %>% 
    mutate_(perc = "length / tot_length") %>% 
    select_("perc") %>% unlist(use.names = FALSE)
})


#### TO CHECK:

getTotalCoverage <- function(object) {
  r <- reduce(object@ranges)
  pdata <- pData(object@sequences)
  100 * sum(width(r))/sum(pdata[, "length"])
}
setGeneric("getTotalCoverage")

getMotifBySeq <- function(object) {
  lapply(object@ranges, function(x) x$motif_name)  # for RangedData.
}
setGeneric("getMotifBySeq")

getMotifArchBySeq <- function(object) {
  sapply(getMotifBySeq(object), function(x) paste(x, collapse = "-"))
}
setGeneric("getMotifArchBySeq")

# convert architectures.
.letter2num <- function(x) {
  .architecture.code <- c(LETTERS, letters, as.character(0:9))
  x <- strsplit(x, "")[[1]]
  as.character(sapply(x, function(z) which(.architecture.code %in% z), USE.NAMES = FALSE))
}

.num2letter <- function(x) {
  .architecture.code <- c(LETTERS, letters, as.character(0:9))
  paste(.architecture.code[as.numeric(x)], collapse = "")
}

convertArch <- function(object, to = "string") {
  to <- match.arg(to, c("number", "string"))
  switch(to, number = sapply(object, .letter2num), string = sapply(object, .num2letter))
}

getMotifDistribution <- function(object, by.groups) {
  m <- object
  if (!missing(by.groups)) 
    rownames(m) <- by.groups
  
  nu <- unique(rownames(m))
  tmp <- lapply(nu, function(n) {
    sel.n <- rownames(m) %in% n
    l <- length(which(sel.n))
    100 * apply(m[sel.n, , drop = FALSE], 2, sum)/l
  })
  names(tmp) <- nu
  do.call(rbind, tmp)
}
setGeneric("getMotifDistribution")

setMethod("getMotifDistribution", "MotifSearchResult", function(object, by.groups) {
  m <- getMotifMatrix(object)
  getMotifDistribution(m, by.groups)
})


getMotifCount <- function(object, percentage = FALSE) {
  apply(getMotifMatrix(object), 2, function(m) {
    if (percentage) 
      100 * sum(m)/nseq(object) else sum(m)
  })
}
setGeneric("getMotifCount")