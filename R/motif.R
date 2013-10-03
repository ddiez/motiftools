getMotifBySeq = function(object) {
  lapply(object@ranges,names)
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
    r=object@ranges[[n]]
    h=unique(names(r))
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

## for old classes.

setMethod("getMotifBySeq","MotifSet",
function(object) {
  res=list()
  for(n in names(object@motif)) {
    tmp=object@motif[[n]]
    #apply(tmp,1,function(x) {
    all_ids = unique(tmp$Id)
    for(id in all_ids) {
      tmp2 = tmp[tmp$Id == id,]
      if(is.null(res[[id]]))
        res[[id]]=data.frame(Motif=rep(n,nrow(tmp2)),Start=tmp2$Start)
      else
        res[[id]]=rbind(res[[id]], data.frame(Motif=rep(n,nrow(tmp2)),Start=tmp2$Start))
    }
  }
  reorderMotifs(res)
})

reorderMotifs = function(M) {
  for(s in names(M)) {
    d = M[[s]]
    #rownames(d) = d[, "Motif"]
    M[[s]] <- d[order(d[, "Start"]), ]
  }
  M
}


setMethod("getMotifMatrix","MotifSet",
function(object) {
  M = object
  r = matrix(0, nrow = M@nseq, ncol = M@nmotif)
  rownames(r) = M@sequence
  colnames(r) = 1:M@nmotif
  for(n in 1:M@nmotif) {
    r[as.character(M@motif[[n]]$Id), n] = 1
  }
  r
})

# getMotifDistribution = function(m) {
#   n = sub("(..).*", "\\1", rownames(m))
#   r = t(apply(m, 2, function(x) {
#     sel = which(x == 1)
#     nn = n[sel]
#     nn = factor(nn, levels = unique(n))
#     table(nn)
#   }))
#   
#   nt = table(n)
#   nc = as.vector(nt)
#   names(nc) = names(nt)
#   rr = sapply(colnames(r), function(x) {
#     100*r[,x]/nc[x]
#   })
#   rr
# }

setMethod("getMotifCount", "MotifSet",
function(object, percentage = FALSE) {
  x = unlist(lapply(object@motif, function(m) length(unique(m$Id))))
  if (percentage)
    100 * x / object@nseq
  else
    x
})
