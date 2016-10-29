#' Read MEME XML file format
#' 
#' Read MEME XML file format
#'
#' @param filename MEME XML file.
#' @param sequenceData sequence data object (optional).
#' @param description dataset description.
#'
#' @return MotifSearchResult object.
#'
#' @examples
#' NULL
readMEME_old <- function(filename, sequenceData, description = NULL) {
  .readAlphabet <- function(top) {
    alphabet <- xmlApply(top[["training_set"]][["alphabet"]], function(s) {
      data.frame(id = xmlGetAttr(s, "id"), symbol = xmlGetAttr(s, "symbol"), stringsAsFactors = TRUE)
    })
    alphabet <- do.call(rbind,alphabet)
    rownames(alphabet) <- alphabet$id
    alphabet
  }
  
  .readProbabilities <- function(x, alphabet) {
    res <- xmlApply(x[["motifs"]], function(m) {
      w <- as.numeric(xmlGetAttr(m, "width"))
      tmp <- t(as.matrix(xmlToDataFrame(m[["probabilities"]][["alphabet_matrix"]], colClasses = rep("numeric", 20), collectNames = FALSE)))
      #rownames(tmp) <- alphabet$symbol
      colnames(tmp) <- 1:w
      tmp
    })
    names(res) <- unlist(xmlApply(top[["motifs"]], function(x) xmlGetAttr(x, "name")))
    res
  }
  
  doc <- xmlParse(filename)
  top <- xmlRoot(doc)
  
  # get sequence ids and number:
  seq_info <- xmlApply(top[["training_set"]], function(s) {
    if (xmlName(s) == "sequence") {
      data.frame(seq_id = xmlGetAttr(s, "id"), seq_name = xmlGetAttr(s, "name"), length = as.integer(xmlGetAttr(s, "length")), weight = as.numeric(xmlGetAttr(s, "weight")), stringsAsFactors = FALSE)
    }
  })
  seq_info <- seq_info[! sapply(seq_info, is.null)]
  seq_info <- do.call(rbind, seq_info)
  rownames(seq_info) <- seq_info$seq_id
  nseq <- nrow(seq_info)
  
  # motifs.
  nmotif <- as.integer(xmlValue(top[["model"]][["nmotifs"]]))
  
  motif_info <- xmlApply(top[["motifs"]], function(m) {
    if (xmlName(m) == "motif")
      data.frame(
        motif_id = xmlGetAttr(m, "id"),
        motif_name = xmlGetAttr(m, "name"),
        width = as.integer(xmlGetAttr(m, "width")),
        stringsAsFactors = FALSE
      )
  })
  motif_info <- do.call(rbind, motif_info)
  rownames(motif_info) <- motif_info$motif_id
  
  motif_alf <- .readAlphabet(top)
  
  motif_list <- .readProbabilities(top, motif_alf)
  
  # sequence hits:
  res <- xmlApply(top[["scanned_sites_summary"]], function(s) {
    if (xmlName(s) == "scanned_sites") {
      seq_id <- xmlGetAttr(s, "sequence_id")
      tmp <- xmlApply(s, function(ss) {
        if (xmlName(ss) == "scanned_site") {
          data.frame(
            seq_id = seq_info[seq_id, "seq_name"],
            motif_name = motif_info[xmlGetAttr(ss, "motif_id"), "motif_name"],
            pos = as.numeric(xmlGetAttr(ss, "position")),
            pvalue = as.numeric(xmlGetAttr(ss, "pvalue")),
            stringsAsFactors = FALSE
          )
        }
      })
      do.call(rbind, tmp)
    }
  })
  res <- do.call(rbind,res)
  rownames(res) = NULL
  
  # sequenceData.
  rownames(seq_info) <- seq_info$seq_name
  sequenceData <- AnnotatedDataFrame(seq_info)
  
  # motifData.
  rownames(motif_info) <- motif_info$motif_name
  motifData <- AnnotatedDataFrame(motif_info)
  
  if (is.null(res)) {
    ranges <- RangedData()
  } else {
    res <- res[with(res, order(seq_id,pos)), ] # reorder by sequence and then position.
    ranges <-
      RangedData(
        IRanges(start = res$pos, width = motif_info[res$motif_name, "width"]),
        motif_name = res$motif_name,
        score = rep(NA, nrow(res)),
        pvalue = res$pvalue,
        qvalue = rep(NA, nrow(res)),
        evalue = rep(NA, nrow(res)),
        space = res$seq_id
      )
  }
  
  new(
    "MotifSearchResult",
    info = list(
      tool = "MEME",
      description = description,
      nseq = nseq,
      nmotif = nmotif
    ),
    sequences = sequenceData,
    motifs = motifData,
    models = motif_list,
    ranges = ranges
  )
}