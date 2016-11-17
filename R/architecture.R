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

# get architechtures encoded as data.frame.
getMotifArchDataFrame <- function(object) {
  tmp <- getMotifsBySeq(object)
  lapply(tmp, function(x) {
    x <- c("H", x, "T")
    z <- lapply(seq_along(head(x, -1)), function(k) {
      data.frame(from = x[k], to = x[k + 1], stringsAsFactors = FALSE)
    })
    z <- bind_rows(z)
  })
}

# get architechtures encoded as graphs (igraph)
getMotifArchGraph <- function(object) {
  lapply(getMotifArchDataFrame(object), graph_from_data_frame)
}

# compute euclidean distance of two matrices or vectors.
eud <- function(x, y) {
  if (any(dim(x) != dim(x)))
    stop("x and y have to have same dimensions.")
  sqrt(sum((x - y) ^ 2))
}


# compute jaccard similarity of two matrices or vectors.
jaccard_sim <- function(x, y) {
  sum(x & y) / sum(x | y)
}

# compute jaccard similarity of two graphs.
graph_sim <- function(x, y) {
  n <- unique(c(V(x)$name, c(V(y)$name)))
  l <- length(n)
  m <- matrix(0, ncol = l, nrow = l, dimnames = list(n, n))
  
  m1 <- as_adj(x, sparse = FALSE)
  m2 <- as_adj(y, sparse = FALSE)
  m1e <- m2e <- m
  
  m1e[rownames(m1), colnames(m1)] <- m1
  m2e[rownames(m2), colnames(m2)] <- m2
  
  jaccard_sim(m1e, m2e)
}

# compute jaccard similarity matrix of a list of graphs.
graph_sim_list <- function(x) {
  l <- length(x)
  comb <- t(combn(names(x), 2))
  m <- matrix(NA, ncol = l, nrow = l, dimnames = list(names(x), names(x)))
  diag(m) <- 1
  for (k in seq_len(nrow(comb))) {
    i <- comb[k, 1]
    j <- comb[k, 2]
    m[i, j] <- graph_sim(x[[i]], x[[j]])
    m[j, i] <- m[i, j]
  }
  m
}

#' getMotifArchSimilarity
#' 
#' Computes similarity between sequence motif architectures.
#' 
#' @param object MotifSearchResult object.
#' 
#' @return A matrix with Jaccard index similarity.
#' @export
getMotifArchSimilarity <- function(object) {
  gl <- getMotifArchGraph(object)
  graph_sim_list(gl)
}

#' plotMotifArchSimilarity
#' 
#' Plots heatmap with similarity between sequence motif architectures.
#' 
#' @param object MotifSearchResult object.
#' 
#' @return NULL
#' @export
plotMotifArchSimilarity <- function(object) {
  m <- object
  h <- hclust(dcor(m))
  m <- m[h$order, h$order]
  d <- reshape2::melt(m) # need to fix this.
  ggplot(d, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile(color = "black") + 
    viridis::scale_fill_viridis("similarity", limits = c(0, 1), guide = guide_legend(reverse = TRUE)) + 
    theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1), plot.title = element_text(hjust = .5)) + 
    labs(x = "", y = "", title = "Architecture similarity") + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))
}