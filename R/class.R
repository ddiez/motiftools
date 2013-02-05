setClass("MotifSet",
         representation(
           nmotif = "numeric",
           motif = "list",
           nseq = "numeric",
           sequence = "character"
         )
)

setClass("MotifLib", representation = list(
  library = "list",
  annotation = "data.frame"
))

setMethod("[", "MotifLib",
          function(x, i, j, drop = TRUE) {
            x@library = x@library[i]
            x@annotation = x@annotation[i,]
            x
          })

setMethod("show", "MotifLib",
          function(object) {
            message("MotifLib object with: ", length(object@library), " motifs.")
          })

setMethod("names", "MotifLib",
          function(x) {
            names(x@library)
          })

setOldClass("dendrogram")
setClass("TomTom", representation(
  matrix = "matrix",
  cutoff.type = "character",
  cutoff = "numeric",
  matrix_key = "matrix",
  color_key = "character",
  dendrogram = "dendrogram"
))
setGeneric("getPWMnames", function(object) standardGeneric("getPWMnames"))
setMethod("getPWMnames", "TomTom",
          function(object) {
            colnames(object@matrix)
          })