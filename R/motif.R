getMotifBySeq = function(object) {
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
}

reorderMotifs = function(M) {
  for(s in names(M)) {
    d = M[[s]]
    #rownames(d) = d[, "Motif"]
    M[[s]] <- d[order(d[, "Start"]), ]
  }
  M
}



## OLD

getMotifDistribution = function(m) {
  n = sub("(..).*", "\\1", rownames(m))
  r = t(apply(m, 2, function(x) {
    sel = which(x == 1)
    nn = n[sel]
    nn = factor(nn, levels = unique(n))
    table(nn)
  }))
  
  nt = table(n)
  nc = as.vector(nt)
  names(nc) = names(nt)
  rr = sapply(colnames(r), function(x) {
    100*r[,x]/nc[x]
  })
  rr
}

getMotifMatrix = function(M) {
  r = matrix(0, nrow = M@nseq, ncol = M@nmotif)
  rownames(r) = M@sequence
  colnames(r) = 1:M@nmotif
  for(n in 1:M@nmotif) {
    r[as.character(M@motif[[n]]$Id), n] = 1
  }
  r
}

getMotifCount = function(M, percentage = FALSE) {
  #x = unlist(lapply(M@motif, nrow))
  x = unlist(lapply(M@motif, function(m) length(unique(m$Id))))
  if (percentage)
    x = 100 * x / M@nseq
  x
}
