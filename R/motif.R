getMotifBySeq = function(object) {
  lapply(object@ranges, function(x) x$motif_name) # for RangedData.
  #lapply(object@ranges,names) # for old IRanges.
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

#getMotifMatrix = function(object,seqnames) {
getMotifMatrix = function(object) {
  #if(missing(seqnames))
  #  seqnames=object@info$sequence_info
  #nseq=length(seqnames)
  
  m = matrix(0,nrow=nseq(object),ncol=nmotif(object))
  rownames(m)=object@info$sequence_info
  colnames(m)=object@info$motif_info$motif_name
  #m = matrix(0,nrow=nseq,ncol=nmotif(object))
  #rownames(m)=seqnames
  #colnames(m)=object@info$motif_info$motif_name
  for(n in names(object@ranges)) {
    #r=object@ranges[[n]] # for old IRanges.
    #h=unique(names(r))
    h=unique(object@ranges[n]$motif_name) # for RangedData.
    m[n,h]=1
  }
  m
}
setGeneric("getMotifMatrix")

getMotifDistribution=function(object,convert.func) {
  m=object
  if(!missing(convert.func))
    rownames(m)=convert.func(rownames(m))
  
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
function(object,convert.func) {
  m=getMotifMatrix(object)
  getMotifDistribution(m,convert.func)
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