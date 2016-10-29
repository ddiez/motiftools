#' MotifSearchResult
#' 
#' MotifSearchResult
#' 
#' @name MotifSearchResult-class
#' @rdname MotifSearchResult-class
#' @slot info list. 
#' @slot sequences AnnotatedDataFrame. 
#' @slot motifs AnnotatedDataFrame. 
#' @slot models list. 
#' @slot ranges RangedData.
#' @param object MotifSearchResult object.
#' @param x MotifSearchResult object.
#' @param i row number.
#' @param j column number.
#' @param ... further arguments passed down to method.
#' @param drop logical; whether to coerce to vector when number of columns equal one.
#' 
setClass("MotifSearchResult",
         representation(
           info = "list",
           sequences = "AnnotatedDataFrame",
           motifs = "AnnotatedDataFrame",
           models = "list",
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
              x@motifs=x@motifs[j,,drop=drop]
              x@models=x@models[j] 
            }
            if(!missing(i)) {
              x@sequences=x@sequences[i,,drop=drop]
            }
            x@info$nseq = unname(nrow(x@sequences))
            x@info$nmotif = unname(nrow(x@motifs))
            x@ranges = x@ranges[space(x@ranges) %in% sequenceNames(x) & x@ranges$motif_name %in% motifNames(x),]
            x
          })

#' @rdname MotifSearchResult-class
#' @aliases nmotif,MotifSearchResult-method
nmotif <- function(object) object@info$nmotif
setGeneric("nmotif")

#' @rdname MotifSearchResult-class
#' @aliases nseq,MotifSearchResult-method
nseq <- function(object) object@info$nseq
setGeneric("nseq")

#' @rdname MotifSearchResult-class
#' @aliases motifNames,MotifSearchResult-method
motifNames <- function(object) featureNames(object@motifs)
setGeneric("motifNames")

#' @rdname MotifSearchResult-class
#' @aliases sequenceNames,MotifSearchResult-method
sequenceNames <- function(object) featureNames(object@sequences)
setGeneric("sequenceNames")


# getMotifMatchMatrix <- function(object, motif, pssm="BLOSUM62") {
#   r=object@ranges
#   r_motif=r[r$motif_name==motif,]
#   
#   best_f=as.character(object@motifs[motif,]$best_f)
#   bs=strsplit(best_f,"")[[1]]
#   
#   ms=do.call(rbind,strsplit(as.character(r_motif$sequence_hit), ""))
#   colnames(ms)=bs
#   rownames(ms)=space(r_motif)
#   
#   m=matrix(0,ncol=nchar(best_f),nrow=nrow(r_motif))
#   colnames(m)=bs
#   rownames(m)=space(r_motif)
#   
#   PSSM=get(data(list=pssm))
#   for(k in 1:length(bs)) {
#     m[,k]=PSSM[bs[k],ms[,k]]
#   }
#   m
# }
# setGeneric("getMotifMatchMatrix")
# 
# .plotBoxPlot <- function(x, ...) {
#   xs=.4
#   xss=.1
#   
#   plot(NA,xlim=c(1-.5,ncol(x)+.5),xaxs="i", ...)
#   abline(h=0,col="grey")
#   for(k in 1:ncol(x)) {
#     tmp=fivenum(x[,k])
#     lines(c(k-xs,k+xs),c(tmp[3],tmp[3]),lwd=2)
#     rect(k-xs,tmp[2],k+xs,tmp[4])
#     lines(c(k,k),c(tmp[2],tmp[1]),lty="dotted")
#     lines(c(k-xss,k+xss),c(tmp[1],tmp[1]),lwd=1)
#     lines(c(k,k),c(tmp[4],tmp[5]),lty="dotted")
#     lines(c(k-xss,k+xss),c(tmp[5],tmp[5]),lwd=1)
#   }
# }
# 
# plotMotifMatchMatrix <- function(object, motif, pssm="BLOSUM62", ...) {
#   l=layout(matrix(c(1,0,2,3),ncol=2,byrow=TRUE),width=c(15,5),height=c(5,15))
#   m=getMotifMatchMatrix(object, motif, pssm=pssm)
#   
#   bmotif=object@motifs[motif,]$best_f
#   bmotif=strsplit(as.character(object@motifs[motif,]$best_f),"")[[1]]
#   max_score=diag(BLOSUM62[bmotif,bmotif])
#   ylim=range(BLOSUM62[bmotif,bmotif]) + c(-1,1)
#   limit=range(BLOSUM62[bmotif,bmotif])
#   
#   cs=colMeans(m)
# #  cs=100*cs/max_score
#     
#   rs=rowMeans(m)
#   
#   mm=t(t(m)/max_score)
#   
#   rso=order(rs)
#   rs=rs[rso]
#   m=m[rso,]
#   mm=mm[rso,]
#   
#   op=par(mar=c(0,2,2,0))
#   .plotBoxPlot(m, xlab="Position",ylab="PSSM score",axes=FALSE,ylim=ylim)
#   axis(2,las=1,lwd=0,line=-.5)
#   points(1:length(max_score),max_score, col="red",pch=19,bg="red")
#   abline(h=limit,col="blue")
#   box()
#   
#   op=par(mar=c(2,2,.5,0))
#   col=colorRampPalette(c("white","blue","orange"))(128)
#   image(1:ncol(mm), 1:nrow(mm),t(mm),main="",axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i",col=col)
#   mtext(colnames(mm),at=1:ncol(mm),side=1,line=.5, ...)
#   box()
#   par(mar=c(2,.5,.5,2))
#   plot(rs,1:length(rs),xaxs="i",yaxs="i",axes=FALSE,type="l",col="red")
#   axis(3,las=1,lwd=0,line=-.5)
#   box()
#   par(op)
# }
# setGeneric("plotMotifMatchMatrix")
# 
# ##### DEPRECATED.
# setClass("MotifSet",
#          representation(
#            nmotif = "numeric",
#            motif = "list",
#            nseq = "numeric",
#            sequence = "character"
#          )
# )
# 
# setClass("MotifLib", representation = list(
#   library = "list",
#   annotation = "data.frame"
# ))
# 
# setMethod("[", "MotifLib",
#           function(x, i, j, drop = TRUE) {
#             x@library = x@library[i]
#             x@annotation = x@annotation[i,]
#             x
#           })
# 
# setMethod("show", "MotifLib",
#           function(object) {
#             message("MotifLib object with: ", length(object@library), " motifs.")
#           })
# 
# setMethod("names", "MotifLib",
#           function(x) {
#             names(x@library)
#           })
# 
# setOldClass("dendrogram")
# setClass("TomTom", representation(
#   matrix = "matrix",
#   cutoff.type = "character",
#   cutoff = "numeric",
#   matrix_key = "matrix",
#   color_key = "character",
#   dendrogram = "dendrogram"
# ))
# setGeneric("getPWMnames", function(object) standardGeneric("getPWMnames"))
# setMethod("getPWMnames", "TomTom",
#           function(object) {
#             colnames(object@matrix)
#           })