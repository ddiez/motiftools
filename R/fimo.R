#' Read results from FIMO motif search stored in XML format.
#' 
#' Read results from FIMO motif search stored in XML format and returns a
#' MotifSearchResult object.
#'
#' @param file file XML output file from FIMO.
#' @param description character string describing the dataset.
#'
#' @return MotifSearchResult object.
#' @export
#'
#' @examples
#' NULL
readFIMO <- function(file, description = NULL) {
  # parse fimo.xml.
  doc <- read_xml(file)
  root <- xml_root(doc)
  
  alfa_data <- getFimoAlphabet(root)
  alfa_bkgd <- getFimoBackground(root)
  alfa_info <- merge(alfa_data, alfa_bkgd, by.x = "id", by.y = "letter")
  
  nseq <- getFimoSequenceInfo(root)$nseq
  motif_info <- getFimoMotifInfo(root)
  nmotif <- nrow(motif_info)
  cisml_file <- getFimoCisml(root)
  cisml_file <- file.path(dirname(file), cisml_file)
  
  # parse cisml.xml.
  doc <- read_xml(cisml_file)
  root <- xml_root(doc)
  
  motif_hits <- getFimoMotifHits(root)
  
  seq_names <- unique(motif_hits$seq_id)
  seq_info <- data.frame(seq_names, sequence_id = seq_names, row.names = 1, stringsAsFactors = FALSE)
  
  if (nrow(motif_hits) > 0) {
    motif_hits <-
      GRanges(
        seqnames = motif_hits$seq_id,
        ranges = IRanges(start = motif_hits$start, width = motif_info[motif_hits$motif_name, "width"]),
        strand = "*",
        motif_name = motif_hits$motif_name,
        score = motif_hits$score,
        pvalue = motif_hits$pvalue,
        qvalue = motif_hits$qvalue,
        evalue = rep(NA, nrow(motif_hits)),
        sequence = motif_hits$sequence,
        space = motif_hits$seq_id
        )
  } else {
    motif_hits <- GRanges()
  }
  
  motifData <- AnnotatedDataFrame(motif_info)
  sequenceData <- AnnotatedDataFrame(seq_info)
  
  new(
    "MotifSearchResult",
    info = list(
      tool = "FIMO",
      description = description,
      nseq = nseq,
      nmotif = nmotif,
      alphabet = alfa_info
    ),
    sequences = sequenceData,
    motifs = motifData,
    ranges = motif_hits
  )
}

# motif info.
getFimoMotifInfo <- function(x) {
  tmp <- do.call(rbind, xml_attrs(xml_find_all(x, "motif")))
  data.frame(
    name = tmp[, "name"],
    alt = tmp[, "alt"],
    width = as.integer(tmp[, "width"]),
    best_match = tmp[, "best-possible-match"],
    stringsAsFactors = FALSE
  )
}

# sequence background.
getFimoBackground <- function(x) {
  tmp <- xml_find_all(x, ".//background/value")
  data.frame(
    letter = xml_attr(tmp, "letter"),
    background = xml_double(tmp),
    stringsAsFactors = FALSE
  )
}

# cisml file.
getFimoCisml <- function(x) {
  xml_text(xml_find_first(x, "cisml-file")) 
}

# alphabet.
getFimoAlphabet <- function(x) {
  tmp <- xml_attrs(xml_find_all(x, ".//letter"))
  tmp <- lapply(tmp, function(z) {
    data.frame(id = z["id"],
               symbol = z["symbol"],
               name = z["name"],
               color = z["colour"],
               aliases = z["aliases"],
               equals = z["equals"],
               row.names = NULL,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, tmp)
}

# sequence info.
getFimoSequenceInfo <- function(x) {
  tmp <- xml_attrs(xml_find_first(x, "sequence-data"))
  data.frame(
    nseq = as.integer(tmp["num-sequences"]),
    nresidues = as.integer(tmp["num-residues"])
  )
}

# motif hits.
getFimoMotifHits <- function(x) {
  tmp_motif <- lapply(xml_children(x), function(motif) { # iterate over motifs ("pattern").
    # check we are in a pattern.
    if (xml_name(motif) == "pattern") {
      motif_name <- xml_attr(motif, "name")
      tmp_seq <- lapply(xml_children(motif), function(seq) { # iterate over sequences.
        # check we are in a sequence.
        if (xml_name(seq) == "scanned-sequence") {
          seq_id <- xml_attr(seq, "accession")
          tmp_match <- lapply(xml_children(seq), function(hit) { # iterate over hits.
            # check we are in a hit.
            if (xml_name(hit) == "matched-element") {
              seq_str <- xml_text(xml_contents(hit)[[1]])
              qvalue <- xml_double(xml_contents(hit)[[2]])
              data.frame(
                motif_name = motif_name,
                seq_id = seq_id,
                start = as.integer(xml_attr(hit, "start")),
                score = as.numeric(xml_attr(hit, "score")),
                pvalue = as.numeric(xml_attr(hit, "pvalue")),
                qvalue = qvalue,
                sequence = seq_str,
                stringsAsFactors = FALSE
              )
            }
          })
          bind_rows(tmp_match)
        }
      })
      bind_rows(tmp_seq)
    }
  })
  bind_rows(tmp_motif)
}
