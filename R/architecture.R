getMotifArchString <- function(object, convert.to.letter = FALSE, return.unique = FALSE) {
  .architecture.code <- c(LETTERS, letters, as.character(0:9))
  tmp <- getMotifsBySeq(object)
  if (convert.to.letter) {
    if (max(as.numeric(unlist(tmp))) > length(.architecture.code)) 
      stop("more motifs than architectures! cannot convert to letters.")
    tmp <- sapply(tmp, .num2letter)
  } else tmp <- sapply(tmp, function(x) paste(x, collapse = "-"))
  if (return.unique) 
    unique(tmp) else tmp
}


#' getMotifsBySeq
#' 
#' Obtain list with a vector of motifs per sequence.
#' 
#' @param object a MotifSearchResult object.
#' 
#' @export
#' @return A list of sequences.
getMotifsBySeq <- function(object) {
  tmp <- as.data.frame(object@ranges)
  split(tmp$motif_name, tmp$seqnames)
}

# simpler version of getMotifArchString
getMotifArchBySeq <- function(object) {
  sapply(getMotifsBySeq(object), function(x) paste(x, collapse = "-"))
}

# convert architectures.
.letter2num <- function(x) {
  .architecture.code <- c(LETTERS, letters, as.character(0:9))
  x <- strsplit(x, "")[[1]]
  as.character(sapply(x, function(z) which(.architecture.code %in% z), USE.NAMES = FALSE))
}

.num2letter <- function(x) {
  .architecture.code <- c(LETTERS, letters, as.character(0:9))
  paste(.architecture.code[as.numeric(x)], collapse = "")
}

convertArch <- function(object, to = "string") {
  to <- match.arg(to, c("number", "string"))
  switch(to, number = sapply(object, .letter2num), string = sapply(object, .num2letter))
}