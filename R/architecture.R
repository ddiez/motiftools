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




### OLD
countArchs = function(M, min.motif = 3, check.coherence = TRUE, model, remove.rep = TRUE) {
  if (check.coherence && missing(model)) stop("model is needed to check coherence")
  A = list()
  for(s in names(M)) {
    #print(s)
    ms = as.character(M[[s]][, "Motif"])
    if (check.coherence)
      ms = ms[ms %in% model]
    #print(ms)
    if (remove.rep)
      msl = unique(ms)
    else
      msl = ms
    if (length(msl) >= min.motif) {
      doit = TRUE
      if (check.coherence) {
        doit = FALSE
        o = sapply(ms, function(m) which(model %in% m))
        #print(o)
        if (! is.unsorted(o))
          doit = TRUE
        if (! is.unsorted(rev(o)))
          doit = TRUE
      }
      
      if (doit) {
        a = paste(ms, collapse = "-")
        #print(a)
        if(length(A[[a]]) == 0)
          A[[a]] = 1
        else
          A[[a]] = A[[a]] + 1
      }
    }
  }
  d = data.frame(Architecture = names(A), Counts = as.numeric(unlist(A)), Percentage = 100 * as.numeric(unlist(A)) / as.numeric(length(names(M))))
  d = d[order(d[, "Counts"], decreasing = TRUE), ]
  d = data.frame(d, Cumulative = cumsum(d$Percentage))
  d
}