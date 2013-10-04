getMotifArchString = function(object, convert.to.letter=FALSE,return.unique=FALSE) {
  .architecture.code=c(LETTERS,letters,as.character(0:9))
  tmp = lapply(getMotifBySeq(object),function(x)x[,1])
  if(convert.to.letter) {
    if(max(as.numeric(unlist(tmp)))>length(.architecture.code)) stop("more motifs than architectures! cannot convert to letters.")
    tmp=sapply(tmp, .num2letter)
  }
  else
    tmp=sapply(tmp, function(x) paste(x, collapse="-"))
  if(return.unique) unique(tmp) else tmp
}