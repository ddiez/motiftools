#' getMotifArchString
#' 
#' Obtain  a vector of motif architectures per sequence, as a string. Optionally
#' encode the motif numbers using letter, LETTER and 0:9 symbols (if the number
#' of motifs is lower than the number of symbols).
#' 
#' @param object a MotifSearchResult object.
#' @param convert.to.letter logical; whether to encode motif numbers as character
#' symbols.
#' @param return.unique logical; whether to return unique architectures.
#' 
#' @export
#' @return A list of sequences.
getMotifArchString <- function(object, convert.to.letter = FALSE, return.unique = FALSE) {
  .codedb <- c(LETTERS, letters, as.character(0:9))
  tmp <- getMotifsBySeq(object)
  if (convert.to.letter) {
    if (nmotif(object) > length(.codedb)) 
      stop("more motifs than architectures! cannot convert to index letters.")
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

# convert architecture from letter to index.
.letter2num <- function(x) {
  .codedb <- c(LETTERS, letters, as.character(0:9))
  x <- strsplit(x, "")[[1]]
  as.character(sapply(x, function(z) which(.codedb %in% z), USE.NAMES = FALSE))
}

# convert architecture from index to letter.
.num2letter <- function(x) {
  .codedb <- c(LETTERS, letters, as.character(0:9))
  paste(.codedb[as.numeric(x)], collapse = "")
}

# convert architectures between types.
convertArch <- function(object, to = "string") {
  to <- match.arg(to, c("number", "string"))
  switch(to, number = sapply(object, .letter2num), string = sapply(object, .num2letter))
}