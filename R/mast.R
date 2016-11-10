#' Read results from MAST motif search stored in XML format.
#' 
#' Read results from MAST motif search stored in XML format and returns a
#' MotifSearchResult object.
#'
#' @param file XML output file from MAST
#' @param description character string describing the dataset.
#'
#' @return MotifSearchResult object.
#' @export
#'
#' @examples
#' NULL
readMAST <- function(file, description = NULL) {
  doc <- read_xml(file)
  root <- xml_root(doc)
  
  motif_info <- getMastMotifInfo(root)
  seq_info <- getMastSequenceInfo(root)
  hit_info <- getMastMotifHits(root)
  
  motifData <- AnnotatedDataFrame(motif_info)
  sequenceData <- AnnotatedDataFrame(seq_info)
  
  nmotif <- nrow(motif_info)
  nseq <- nrow(seq_info)
  
  rangeData <- GRanges(
    seqnames = hit_info$seq_id,
    strand = "*",
    ranges = IRanges(start = hit_info$start, width = motif_info[hit_info$motif_id, "width"]),
    motif_name = hit_info$motif_id,
    score = rep(NA, nrow(hit_info)),
    p_value = hit_info$pvalue
  )
  
  new(
    "MotifSearchResult",
    info = list(
      tool = "MAST",
      description = description,
      nseq = nseq,
      nmotif = nmotif
    ),
    sequences = sequenceData,
    motifs = motifData,
    ranges = rangeData
  )
}

# motif info.
getMastMotifInfo <- function(x) {
  tmp <- xml_attrs(xml_find_all(x, ".//motifs/motif"))
  
  motif_info <- lapply(tmp, function(x) {
    tmp <- data.frame(
      db = x["db"],
      id = x["id"],
      width = as.integer(x["length"]),
      nsites = as.integer(x["nsites"]),
      evalue = as.numeric(x["evalue"]),
      bad = FALSE,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    
    if ("bad" %in% names(x))
      tmp$bad <- TRUE
    
    tmp
  })
  motif_info <- do.call(rbind, motif_info)
  rownames(motif_info) <- motif_info$id
  motif_info
}

# sequence info.
getMastSequenceInfo <- function(x) {
  tmp <- do.call(rbind, xml_attrs(xml_find_all(x, ".//sequences/sequence")))
  mode(tmp) <- "character"
  tmp2 <- do.call(rbind, xml_attrs(xml_find_all(x, ".//sequences/sequence/score")))
  mode(tmp2) <- "character"
  data.frame(
    db = tmp[, "db"],
    name = tmp[, "name"],
    comment = tmp[, "comment"],
    length = as.integer(tmp[, "length"]),
    strand = tmp2[, "strand"],
    combined_pvalue = as.numeric(tmp2[, "combined_pvalue"]),
    evalue = as.numeric(tmp2[, "evalue"]),
    row.names = "name",
    stringsAsFactors = FALSE
  )
}

# motif hits.
getMastMotifHits <- function(x) {
  seq_data <- lapply(xml_find_all(x, ".//sequences/sequence"), function(sequence) {
    hit_data <- lapply(xml_find_all(sequence, ".//seg/hit"), function(hit) {
      tmp <- xml_attrs(hit)
      data.frame(
        seq_id = xml_attr(sequence, "name"),
        motif_id = as.character(as.integer(tmp["idx"]) + 1),
        start = as.integer(tmp["pos"]),
        pvalue = as.numeric(tmp["pvalue"]),
        stringsAsFactors = FALSE
      )
    })
    bind_rows(hit_data)
  })
  bind_rows(seq_data)
  
}
