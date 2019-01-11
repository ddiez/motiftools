# handy function.
getMemeAttr <- function(id, attr, info) {
  info[info["id"] == id, attr]
}
getMemeAttr <- Vectorize(getMemeAttr, vectorize.args = "id", USE.NAMES = FALSE)

# get alphabet.
getMemeAlphabet <- function(x) {
  p <- xml_find_first(x, ".//training_set/alphabet")
  d <- lapply(1:xml_length(p), function(i) {
    tmp <- xml_attrs(xml_child(p, i))
    data.frame(id = tmp["id"], symbol = tmp["symbol"], name = tmp["name"], stringsAsFactors = FALSE, row.names = tmp["id"])
  })
  do.call(rbind, d)
}

# get motifs.
getMemeMotifs <- function(x) {
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
getMemeSequences <- function(x) {
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
getMemeSeqStats <- function(x) {
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

# get scanned sites hits.
getMemeScannedSiteHits <- function(x, motif_info) {
  p <- xml_find_all(x, ".//scanned_sites")
  
  hits <- lapply(p, function(node) {
    seq_info <- xml_attr(node, "sequence_id")
    tmp <- do.call(rbind, xml_attrs(xml_children(node)))
    if (!is.null(tmp)) {
      data.frame(
        sequence_id = seq_info,
        motif_id = tmp[, "motif_id"],
        position = as.integer(tmp[, "position"]) + 1,
        strand = tmp[, "strand"],
        pvalue = as.numeric(tmp[, "pvalue"]),
        stringsAsFactors = FALSE
      )  
    } else {
      data.frame()
    }
    
  })
  hits <- do.call(rbind, hits)
  hits$width <- getMemeAttr(hits$motif_id, attr = "width", info = motif_info)
  hits
}

# get motif hits.
getMemeMotifHits <- function(x, motif_info) {
  p <- xml_find_all(x, ".//contributing_sites")
  
  hits <- lapply(seq_along(p), function(k) {
    node <- p[[k]]
    motif_id <- motif_info[k, "id"]
    tmp <- do.call(rbind, xml_attrs(xml_children(node)))
    if (!is.null(tmp)) {
      data.frame(
        sequence_id = tmp[, "sequence_id"],
        motif_id = motif_id,
        position = as.integer(tmp[, "position"]) + 1,
        strand = tmp[, "strand"],
        pvalue = as.numeric(tmp[, "pvalue"]),
        stringsAsFactors = FALSE
      )  
    } else {
      data.frame()
    }
    
  })
  hits <- do.call(rbind, hits)
  hits$width <- getMemeAttr(hits$motif_id, attr = "width", info = motif_info)
  hits
}


# get motif scores.
getMemeMotifScores <- function(x, alphabet = NULL) {
  motifs <- xml_find_first(x, ".//motifs")
  scores <- lapply(seq_len(xml_length(motifs)), function(k) {
    width <- as.integer(xml_attr(xml_child(motifs, k), "width"))
    nodes <- xml_find_all(xml_child(motifs, k), "scores/alphabet_matrix/alphabet_array/value")
    aa <- xml_attrs(nodes) %>% unlist(use.names = FALSE)
    aa <- apply(matrix(aa, ncol = width), 1, unique)
    if (!is.null(alphabet))
      aa <- alphabet[aa, "symbol"]
    scores <- xml_integer(nodes) # not working? (returns "numeric").
    scores <- matrix(scores, ncol = width, dimnames = list(aa, seq_len(width)))
    mode(scores) <- "integer"
    scores
  })
  names(scores) <- xml_attr(xml_find_all(x, ".//motif"), "id")
  scores
}

# get motif probabilities.
getMemeMotifProbabilities <- function(x, alphabet = NULL) {
  motifs <- xml_find_first(x, ".//motifs")
  prob <- lapply(seq_len(xml_length(motifs)), function(k) {
    width <- as.integer(xml_attr(xml_child(motifs, k), "width"))
    nodes <- xml_find_all(xml_child(motifs, k), "probabilities/alphabet_matrix/alphabet_array/value")
    aa <- xml_attrs(nodes) %>% unlist(use.names = FALSE)
    aa <- apply(matrix(aa, ncol = width), 1, unique)
    if (!is.null(alphabet))
      aa <- alphabet[aa, "symbol"]
    prob <- xml_double(nodes)
    prob <- matrix(prob, ncol = width, dimnames = list(aa, seq_len(width)))
    prob
  })
  names(prob) <- xml_attr(xml_find_all(x, ".//motif"), "id")
  prob
}

#' Read results from MEME motif search stored in XML format.
#' 
#' Read results from MEME motif search stored in XML format and returns a
#' MotifSearchResult object.
#'
#' @param file XML output file from MEME.
#' @param description character string describing the dataset.
#'
#' @return MotifSearchResult object.
#' @export
#'
#' @examples
#' NULL
readMEME <- function(file, description = NULL) {
  doc <- read_xml(file)
  root <- xml_root(doc)
  
  alf_info <- getMemeAlphabet(root)
  seq_data <- getMemeSequences(root)
  seq_stats <- getMemeSeqStats(root)
  seq_info <- merge(seq_data, seq_stats, by.x = "id", by.y = "sequence_id")
  rownames(seq_info) <- rownames(seq_data)
  motif_info <- getMemeMotifs(root)
  motif_hit <- getMemeMotifHits(root, motif_info)
  motif_score <- getMemeMotifScores(root, alf_info)
  motif_prob <- getMemeMotifProbabilities(root, alf_info)
  
  # alphabetData.
  alphabetData <- AnnotatedDataFrame(alf_info)
  
  # sequenceData.
  sequenceData <- AnnotatedDataFrame(seq_info)
  
  # motifData.
  motifData <- AnnotatedDataFrame(motif_info)
  
  # rangeData.
  rangeData <- GRanges(
    seqnames = getMemeAttr(motif_hit$sequence_id, attr = "name", info = seq_info),
    strand = "*",
    ranges = IRanges(start = motif_hit$position, width = motif_hit$width),
    motif_name = getMemeAttr(motif_hit$motif_id, attr = "name", info = motif_info),
    p_value = motif_hit$pvalue
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
