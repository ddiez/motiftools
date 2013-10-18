
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

sequenceNames=function(object) featureNames(object@sequences)
setGeneric("sequenceNames")


getMotifMatchMatrix = function(object, motif, pssm="BLOSUM62") {
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
  
  PSSM=get(data(list=pssm))
  for(k in 1:length(bs)) {
    m[,k]=PSSM[bs[k],ms[,k]]
  }
  m
}
setGeneric("getMotifMatchMatrix")

plotMotifMatchMatrix = function(object, motif, pssm="BLOSUM62", ylim=NULL,...) {
  l=layout(matrix(c(1,2),ncol=2),width=c(10,10))
  op=par(mar=c(2,1,2,1))
  m=getMotifMatchMatrix(object, motif, pssm=pssm)
  s=rowSums(m)
  so=order(s)
  s=s[so]
  m=m[so,]
  image(t(m),main=paste("motif",motif),axes=FALSE,xlab="",ylab="",ylim=ylim)
  mtext(colnames(m),at=1:ncol(m),side=1, ...)
  box()
  barplot(s,horiz=TRUE,names=NA)
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