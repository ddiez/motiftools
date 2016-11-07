#' MotifSearchResult
#' 
#' MotifSearchResult
#' 
#' @name MotifSearchResult-class
#' @rdname MotifSearchResult-class
#' @slot info list with result information.
#' @slot sequences AnnotatedDataFrame with sequence information.
#' @slot motifs AnnotatedDataFrame with motif information.
#' @slot probabilities list of motif probabilities. 
#' @slot scores list of motif scores.
#' @slot ranges RangedData with motif hits.
#' @param object MotifSearchResult object.
#' @param x MotifSearchResult object.
#' @param i row number (motif).
#' @param j column number (sequence).
#' @param ... further arguments passed down to method.
#' @param drop logical; whether to coerce to vector when number of columns equal one.
#' 
setClass("MotifSearchResult",
         representation(
           info = "list",
           sequences = "AnnotatedDataFrame",
           motifs = "AnnotatedDataFrame",
           probabilities = "list",
           scores = "list",
           ranges = "RangedData"
         )
)

#' @rdname MotifSearchResult-class
#' @aliases show,MotifSearchResult-method
setMethod("show", "MotifSearchResult",
          function(object) {
            message(class(object), " object.")
            message("Motif tool: ", object@info$tool)
            message("Number of sequences: ", object@info$nseq)
            message("Number of motifs: ", object@info$nmotif)
          })

#' @rdname MotifSearchResult-class
#' @aliases [,MotifSearchResult-method
setMethod("[", "MotifSearchResult",
          function(x, i, j, ..., drop = FALSE) {
            if(!missing(j)) {
              x@motifs <- x@motifs[j, , drop = drop]
              x@probabilities <- x@probabilities[j]
              x@scores <- x@scores[j] 
            }
            if(!missing(i)) {
              x@sequences <- x@sequences[i, , drop = drop]
            }
            x@info$nseq <- unname(nrow(x@sequences))
            x@info$nmotif <- unname(nrow(x@motifs))
            x@ranges <- x@ranges[space(x@ranges) %in% sequenceNames(x) &
                                 x@ranges$motif_name %in% motifNames(x), ]
            x
          })

#' @rdname MotifSearchResult-class
#' @aliases nmotif,MotifSearchResult-method
#' @export
nmotif <- function(object) object@info$nmotif
setGeneric("nmotif")

#' @rdname MotifSearchResult-class
#' @aliases nseq,MotifSearchResult-method
#' @export
nseq <- function(object) object@info$nseq
setGeneric("nseq")

#' @rdname MotifSearchResult-class
#' @aliases motifNames,MotifSearchResult-method
#' @export
motifNames <- function(object) featureNames(object@motifs)
setGeneric("motifNames")

#' @rdname MotifSearchResult-class
#' @aliases sequenceNames,MotifSearchResult-method
#' @export
sequenceNames <- function(object) featureNames(object@sequences)
setGeneric("sequenceNames")

#' @rdname MotifSearchResult-class
#' @aliases scores,MotifSearchResult-method
#' @export
scores <- function(object) object@scores
setGeneric("scores")

#' @rdname MotifSearchResult-class
#' @aliases pwm,MotifSearchResult-method
#' @export
pwm <- function(object) object@probabilities
setGeneric("pwm")

##### MotifCompareResult
#' MotifCompareResult
#' 
#' MotifCompareResult
#' 
#' @name MotifCompareResult-class
#' @rdname MotifCompareResult-class
#' @slot info list with result information.
#' @slot alphabet data.frame with alphabet information.
#' @slot query data.frame with query motifs information.
#' @slot target data.frame with target motifs information.
#' @slot probabilities list of motif probabilities. 
#' @slot matches list of motif scores.
#' 
setClass("MotifCompareResult",
         slots = c(
           info = "list",
           alphabet = "data.frame",
           query = "data.frame",
           target = "data.frame",
           probabilities = "list",
           matches = "data.frame"
         )
)

#' @rdname MotifCompareResult-class
#' @aliases show,MotifCompareResult-method
setMethod("show", "MotifCompareResult",
          function(object) {
            message(class(object), " object.")
            message("Motif tool: ", object@info$tool)
            message("Number of query motifs: ", object@info$nquery)
            message("Number of target motifs: ", object@info$ntarget)
          })

