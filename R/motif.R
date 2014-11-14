#' @title
#' @param object
getCoverage = function(object) {
  r=reduce(object@ranges)
  pdata=pData(object@sequences)
  sapply(names(r), function(n) {
    100*sum(width(r[n]))/pdata[n,"length"]
  })
}
setGeneric("getCoverage")

#'
#'
getTotalCoverage = function(object) {
  r=reduce(object@ranges)
  pdata=pData(object@sequences)
  100*sum(width(r))/sum(pdata[,"length"])
}
setGeneric("getTotalCoverage")

#'
getMotifBySeq = function(object) {
  lapply(object@ranges, function(x) x$motif_name) # for RangedData.
}
setGeneric("getMotifBySeq")

getMotifArchBySeq = function(object) {
  sapply(getMotifBySeq(object),function(x) paste(x,collapse="-"))
}
setGeneric("getMotifArchBySeq")

# convert architectures.
.letter2num = function(x){
  .architecture.code=c(LETTERS,letters,as.character(0:9))
  x=strsplit(x,"")[[1]]
  as.character(sapply(x,function(z) which(.architecture.code %in% z),USE.NAMES=FALSE))
}

.num2letter = function(x){
  .architecture.code=c(LETTERS,letters,as.character(0:9))
  paste(.architecture.code[as.numeric(x)],collapse="")
}

convertArch = function(object,to="string") {
  to=match.arg(to,c("number","string"))
  switch(to,
         "number"=sapply(object,.letter2num),
         "string"=sapply(object,.num2letter)
  )
}

#' @title 
#' @param object a MotifSearchResult object.
#' @param simplify whether to return an presence/absence matrix.
#' @param filter whether remove columns/rows with all zero entries.
getMotifMatrix = function(object, simplify=TRUE, filter=FALSE) {
#   r=object@ranges
#   tmp=data.frame(seq_id=space(r), motif_name=r$motif_name, stringsAsFactors = FALSE)    
#   tmp=unclass(table(tmp))
#   if(simplify)
#     tmp[tmp!=0]=1
#   tmp
  r=object@ranges
  tmp=data.frame(sequence=as.character(space(r)),
                 motif=r$motif_name,stringsAsFactors = FALSE)
  tmp$sequence <- factor(tmp$sequence, levels = sequenceNames(object))
  tmp$motif <- factor(tmp$motif, levels = motifNames(object))
  tmp <- table(tmp)
  if(simplify)
    tmp[tmp!=0]=1
  if(filter)
    tmp[rowSums(tmp)>0, colSums(tmp)>0]
  else
    tmp
}
setGeneric("getMotifMatrix")

getMotifDistribution=function(object,by.groups) {
  m=object
  if(!missing(by.groups))
    rownames(m)=by.groups
  
  nu=unique(rownames(m))
  tmp=lapply(nu,function(n) {
    sel.n=rownames(m) %in% n
    l=length(which(sel.n))
    100*apply(m[sel.n,,drop=FALSE],2,sum)/l
  })
  names(tmp)=nu
  do.call(rbind,tmp)
}
setGeneric("getMotifDistribution")

setMethod("getMotifDistribution","MotifSearchResult",
function(object,by.groups) {
  m=getMotifMatrix(object)
  getMotifDistribution(m,by.groups)
})


getMotifCount = function(object, percentage = FALSE) {
  apply(getMotifMatrix(object),2,function(m) {
    if(percentage)
      100*sum(m)/nseq(object)
    else
      sum(m)
  })
}
setGeneric("getMotifCount")