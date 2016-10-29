# handy function.
getAttr <- function(id, attr, info) {
  info[info["id"] == id, attr]
}
getAttr <- Vectorize(getAttr, vectorize.args = "id", USE.NAMES = FALSE)

# get alphabet.
getAlphabet <- function(x) {
  p <- xml_find_first(x, ".//training_set/alphabet")
  d <- lapply(1:xml_length(p), function(i) {
    tmp <- xml_attrs(xml_child(p, i))
    data.frame(id = tmp["id"], symbol = tmp["symbol"], name = tmp["name"], stringsAsFactors = FALSE, row.names = tmp["id"])
  })
  do.call(rbind, d)
}

# get motifs.
getMotifs <- function(x) {
  p <- xml_find_first(x, ".//motifs")
  d <- lapply(1:xml_length(p), function(i) {
    tmp <- xml_attrs(xml_child(p, i))
    data.frame(
      id = tmp["id"],
      name = tmp["name"],
      width = as.integer(tmp["width"]),
      sites = as.integer(tmp["sites"]),
      ic = as.numeric(tmp["ic"]),
      re = as.numeric(tmp["re"]),
      llr = as.numeric(tmp["llr"]),
      e_value = as.numeric(tmp["e_value"]),
      stringsAsFactors = FALSE,
      row.names = tmp["name"]
    )
  })
  do.call(rbind, d)  
}

# get sequences.
getSequences <- function(x) {
  p <- xml_find_all(x, ".//sequence")
  tmp <- do.call(rbind, xml_attrs(p))
  data.frame(
    id = tmp[, "id"],
    name = tmp[, "name"],
    length = as.numeric(tmp[, "length"]),
    weight = as.numeric(tmp[, "weight"]),
    stringsAsFactors = FALSE,
    row.names = tmp[, "name"]
  )
}

# further sequence info (nsites and p.value).
getSeqStats <- function(x) {
  p <- xml_find_all(x, ".//scanned_sites")
  
  tmp <- do.call(rbind, xml_attrs(p))
  data.frame(
    sequence_id = tmp[, "sequence_id"],
    p.value = as.numeric(tmp[, "pvalue"]),
    nsites = as.integer(tmp[, "num_sites"]),
    stringsAsFactors = FALSE,
    row.names = tmp[, "sequence_id"]
  )
}

# get motif hits.
getMotifHits <- function(x, motif_info) {
  p <- xml_find_all(x, ".//scanned_sites")
  
  hits <- lapply(p, function(node) {
    seq_info <- xml_attr(node, "sequence_id")
    tmp <- do.call(rbind, xml_attrs(xml_children(node)))
    data.frame(
      sequence_id = seq_info,
      motif_id = tmp[, "motif_id"],
      position = as.integer(tmp[, "position"]) + 1,
      strand = tmp[, "strand"],
      p.value = as.numeric(tmp[, "pvalue"]),
      stringsAsFactors = FALSE
    )
  })
  hits <- do.call(rbind, hits)
  hits$width <- getAttr(hits$motif_id, attr = "width", info = motif_info)
  #hits$width <- motif_info[hits$motif_id, "width"]
  hits
}

# get motif scores.
getMotifScores <- function(x) {
  motifs <- xml_find_first(x, ".//motifs")
  scores <- lapply(seq_len(xml_length(motifs)), function(k) {
    width <- as.integer(xml_attr(xml_child(motifs, k), "width"))
    nodes <- xml_find_all(xml_child(motifs, k), "scores/alphabet_matrix/alphabet_array/value")
    aa <- xml_attrs(nodes) %>% unlist(use.names = FALSE)
    aa <- apply(matrix(aa, ncol = width), 1, unique)
    scores <- xml_integer(nodes)
    scores <- matrix(scores, ncol = width, dimnames = list(aa, seq_len(width)))
    scores
  })
  names(scores) <- xml_attr(xml_find_all(x, ".//motif"), "id")
  scores
}

# get motif probabilities.
getMotifProbabilities <- function(x) {
  motifs <- xml_find_first(x, ".//motifs")
  prob <- lapply(seq_len(xml_length(motifs)), function(k) {
    width <- as.integer(xml_attr(xml_child(motifs, k), "width"))
    nodes <- xml_find_all(xml_child(motifs, k), "probabilities/alphabet_matrix/alphabet_array/value")
    aa <- xml_attrs(nodes) %>% unlist(use.names = FALSE)
    aa <- apply(matrix(aa, ncol = width), 1, unique)
    prob <- xml_double(nodes)
    prob <- matrix(prob, ncol = width, dimnames = list(aa, seq_len(width)))
    prob
  })
  names(prob) <- xml_attr(xml_find_all(x, ".//motif"), "id")
  prob
}

#' readMEME
#' 
#' readMEME
#'
#' @param file XML output file from MEME.
#'
#' @return MotifSearchResult object.
#' @export
#'
#' @examples
#' NULL
readMEME <- function(file) {
  doc <- read_xml(file)
  root <- xml_root(doc)
  
  alf_info <- getAlphabet(root)
  seq_info <- getSequences(root)
  motif_info <- getMotifs(root)
  motif_hit <- getMotifHits(root, motif_info)
  motif_score <- getMotifScores(root)
  motif_prob <- getMotifProbabilities(root)
  
  # alphabetData.
  alphabetData <- AnnotatedDataFrame(alf_info)
  
  # sequenceData.
  sequenceData <- AnnotatedDataFrame(seq_info)
  
  # motifData.
  motifData <- AnnotatedDataFrame(motif_info)
  
  # rangeData.
  rangeData <- RangedData(
    IRanges(start = motif_hit$position, width = motif_hit$width),
    motif_name = getAttr(motif_hit$motif_id, attr = "name", info = motif_info),
    p.value = motif_hit$p.value,
    space = getAttr(motif_hit$sequence_id, attr = "name", info = seq_info)
  )
  
  new(
    "MotifSearchResult",
    info = list(
      tool = "MEME",
      description = description,
      nseq = nrow(seq_info),
      nmotif = nrow(motif_info),
      alphabet = alphabetData
    ),
    sequences = sequenceData,
    motifs = motifData,
    probabilities = motif_prob,
    scores = motif_score,
    ranges = rangeData
  )
}
