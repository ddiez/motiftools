
setClass("MotifSearchResult",
         representation(
           info = "list",
           sequences = "AnnotatedDataFrame",
           motifs = "AnnotatedDataFrame",
           ranges = "RangedData"
           )
         )

setMethod("show", "MotifSearchResult",
          function(object) {
            message(class(object), " object.")
            message("Motif tool: ", object@info$tool)
            message("Number of sequences: ", object@info$nseq)
            message("Number of motifs: ", object@info$nmotif)
          })


nmotif = function(object) object@info$nmotif
setGeneric("nmotif")

nseq = function(object) object@info$nseq
setGeneric("nseq")

motifNames=function(object) featureNames(object@motifs)
setGeneric("motifNames")

sequenceNames=function(object) sequenceNames(object@sequences)
setGeneric("sequenceNames")


getMotifMatchMatrix = function(object, motif, extend=TRUE) {
  r=object@ranges
  r_motif=r[r$motif_name==motif,]
  
  best_f=as.character(object@motifs[motif,]$best_f)
  bs=strsplit(best_f,"")[[1]]
  
  ms=do.call(rbind,strsplit(as.character(r_motif$sequence_hit), ""))
  colnames(ms)=bs
  rownames(ms)=space(r_motif)
  
  m=matrix(0,ncol=nchar(best_f),nrow=nrow(r_motif))
  colnames(m)=bs
  rownames(m)=space(r_motif)
  
  data(BLOSUM62)
  for(k in 1:length(bs)) {
    m[,k]=BLOSUM62[bs[k],ms[,k]]
  }
  
  #for(k in 1:length(bs)) {
  #  m[ms[,k]==bs[k],k]=1
  #}
  m
}
setGeneric("getMotifMatchMatrix")

plotMotifMatchMatrix = function(object, motif, ...) {
  op=par(mar=c(2,1,2,1))
  m=getMotifMatchMatrix(object, motif)
  image(x=1:ncol(m),1:nrow(m),t(m[order(rowSums(m)),]),main=paste("motif",motif),axes=FALSE,xlab="",ylab="")
  mtext(colnames(m),at=1:ncol(m),side=1, ...)
  box()
  par(op)
}
setGeneric("plotMotifMatchMatrix")

##### DEPRECATED.
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