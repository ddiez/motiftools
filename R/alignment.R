#' @title align
#' @description aligns two sequences x1 and x2.
#' @param x1 sequence (a character vector).
#' @param x2 sequence (a character vector).
#' @param score.matrix substitution matrix.
#' @param gap.score gap score.
#' @param type type of alignment (global, local), default to local.
align <- function(x1, x2, score.matrix, gap.score=-1, type="local", debug = TRUE) {
  if(length(x1) == 1 && length(x2) == 1) {
    x1 = strsplit(x1, "")[[1]]
    x2 = strsplit(x2, "")[[1]]
  }
  
  type <- match.arg(type, c("global", "local"))
  
  switch(type,
         "local" = {
           .sw(x1, x2, score.matrix=score.matrix, gap.score=gap.score, debug = debug) 
         },
         "global" = {
           .nw(x1, x2, score.matrix=score.matrix, gap.score=gap.score, debug=debug)
         })
}

.checkMax <- function(x) x[which.max(x)]

.generateScoreMatrix <- function(x,match.score=1,mismatch.score=-1) {
  x <- unique(x)
  m <- matrix(mismatch.score, ncol=length(x),nrow=length(x),dimnames = list(x,x))
  diag(m) <- match.score
  m
}

.sw <- function(x1, x2, score.matrix, gap.score=-1,debug=FALSE) {
  l1 = length(x1)
  l2 = length(x2)
  
  if(missing(score.matrix))
    score.matrix <- .generateScoreMatrix(c(x1,x2))
  
  # initialize.
  m <- matrix(NA, ncol=l1+1, nrow=l2+1)
  colnames(m) = c("", x1)
  rownames(m) = c("", x2)
  m[1,]=0
  m[,1]=0
  
  # fill
  res = list()
  res = c(res, lapply(2:ncol(m), function(j) data.frame(row=1, col=j, ii=1, jj=j-1, type="h")))
  res = c(res, lapply(2:nrow(m), function(i) data.frame(row=i, col=1, ii=i-1, jj=1, type="v")))
  for(i in 2:nrow(m)) {
    for(j in 2:ncol(m)) {
      # diagonal.
      d <- m[i-1,j-1] + score.matrix[x1[j-1],x2[i-1]]
      # horizonal
      h <- m[i,j-1] + gap.score
      # vertical
      v <- m[i-1,j] + gap.score
      # assign.
      tmp = .checkMax(c(diag=d, hori=h, vert=v))
      if(tmp<0) tmp[]<-0 # no negative scores.
      if(tmp>0) {
        switch(names(tmp),
               diag = {
                 res = c(res, list(data.frame(row=i, col=j, ii=i-1, jj=j-1, type="d")))
               },
               hori = {
                 res = c(res, list(data.frame(row=i, col=j, ii=i, jj=j-1, type="h")))
               },
               vert = {
                 res = c(res, list(data.frame(row=i, col=j, ii=i-1, jj=j, type="v")))
               }
        )
      }
      m[i,j] <- tmp
    }
  }
  res = do.call(rbind, res)
  
  # backtrace
  index = which(m == max(m), arr.ind = TRUE)
  index<-index[nrow(index),] # there may be several max-- get the last one!
  i <- index[1]
  j <- index[2]
  r1 = c()
  r2 = c()
  st <- 0
  repeat {
    if(m[i,j]==0) break
    
    tmp <- subset(res, row==i & col==j)
    switch(as.character(tmp$type),
           d = {
             r1 = c(r1, x1[j-1])
             r2 = c(r2, x2[i-1])
             st = st + m[i,j]
             i=i-1
             j=j-1
           },
           v = {
             r1 = c(r1, "-")
             r2 = c(r2, x2[i-1])
             i=i-1
           },
           h = {
             r1 = c(r1, x1[j-1])
             r2 = c(r2, "-")
             j=j-1
           }
    )
  }
  
  # build alignment.
  if(length(r1)>0) {
    al=data.frame(rbind(rev(r1),rev(r2)),row.names = c("s1","s2"),check.names = FALSE)
    colnames(al) = 1:ncol(al)
  }
  else
    al=data.frame()
  
  # results.
  if(debug)
    list(score.matrix=score.matrix, gap.score=gap.score, scores=m, traceback=res, sequences=list(s1=x1,s2=x2), score=st, alignment=al)
  else
    list(alignment=al, score=st)
}


.nw <- function(x1, x2, score.matrix, gap.score=-1, debug=FALSE) {
  l1 = length(x1)
  l2 = length(x2)
  
  if(missing(score.matrix))
    score.matrix <- .generateScoreMatrix(c(x1,x2))
  
  # initialize.
  m <- matrix(NA, ncol=l1+1, nrow=l2+1)
  colnames(m) = c("", x1)
  rownames(m) = c("", x2)
  m[1,]=-1*(0:l1)
  m[,1]=-1*(0:l2)
  
  # fill
  res = list()
  res = c(res, lapply(2:ncol(m), function(j) data.frame(row=1, col=j, ii=1, jj=j-1, type="h")))
  res = c(res, lapply(2:nrow(m), function(i) data.frame(row=i, col=1, ii=i-1, jj=1, type="v")))
  for(i in 2:nrow(m)) {
    for(j in 2:ncol(m)) {
      # diagonal.
      d <- m[i-1,j-1] + score.matrix[x1[j-1],x2[i-1]]
      # horizonal
      h <- m[i,j-1] + gap.score
      # vertical
      v <- m[i-1,j] + gap.score
      # assign.
      tmp = .checkMax(c(diag=d, hori=h, vert=v))
      switch(names(tmp),
             diag = {
               res = c(res, list(data.frame(row=i, col=j, ii=i-1, jj=j-1, type="d")))
             },
             hori = {
               res = c(res, list(data.frame(row=i, col=j, ii=i, jj=j-1, type="h")))
             },
             vert = {
               res = c(res, list(data.frame(row=i, col=j, ii=i-1, jj=j, type="v")))
             }
      )
      m[i,j] <- tmp
    }
  }
  res = do.call(rbind, res)
  
  # backtrace
  i <- nrow(m)
  j <- ncol(m)
  r1 = c()
  r2 = c()
  st <- 0
  repeat {
    tmp <- subset(res, row==i & col==j)
    switch(as.character(tmp$type),
           d = {
             r1 = c(r1, x1[j-1])
             r2 = c(r2, x2[i-1])
             st = st + m[i,j]
             i=i-1
             j=j-1
           },
           v = {
             r1 = c(r1, "-")
             r2 = c(r2, x2[i-1])
             i=i-1
           },
           h = {
             r1 = c(r1, x1[j-1])
             r2 = c(r2, "-")
             j=j-1
           }
    )
    if(i==1 && j==1) break
  }
  
  # build alignment.
  al=data.frame(rbind(rev(r1),rev(r2)),row.names = c("s1","s2"),check.names = FALSE)
  
  # results.
  if(debug)
    list(score.matrix=score.matrix, gap.score=gap.score, scores=m, traceback=res, sequences=list(s1=x1,s2=x2), score=st, alignment=al)
  else
    list(alignment=al, score=st)
}

motifAlign <- function(x1, x2, align.method="local", score.method="PCC", debug=TRUE) {
  score.method <- match.arg(score.method, c("PCC","euclid"))
  
  # rename columns for motifs A and B.
  m1 <- x1
  m2 <- x2
  colnames(m1) <- paste0("A",1:ncol(m1))
  colnames(m2) <- paste0("B",1:ncol(m2))
  
  # generate motifs column scoring matrix.
  switch(score.method,
         "PCC" = {
           s <- cor(m1,m2)
         },
         "euclid" = {
           m <- cbind(m1,m2)
           s <- as.matrix(dist(t(m), diag = TRUE, upper = TRUE))
         })
  
  # align columns.
  align(colnames(m1), colnames(m2), score.matrix = s, debug=debug, type = align.method)
}

motifDistance <- function(x1, x2, score.method="PCC", align.method="local", debug=TRUE) {
  motifAlign(x1,x2,score.method=score.method, align.method=align.method, debug=debug)$score
}
