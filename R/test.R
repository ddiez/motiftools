readMEME2 <- function(file, description = NULL) {
  library(xml2)
  #f <- "~/projects/motiftools/inst/files/meme_ras/meme.xml"
  doc <- read_xml(file)
  root <- xml_root(doc)
  
  # alphabet
  p <- xml_find_first(root, ".//training_set/alphabet")
  d <- lapply(1:xml_length(p), function(i) {
    tmp <- xml_attrs(xml_child(p, i))
    data.frame(id = tmp["id"], symbol = tmp["symbol"], name = tmp["name"], stringsAsFactors = FALSE, row.names = NULL)
  })
  alphabet <- do.call(rbind, d)
  #alphabet
  
  # motifs.
  p <- xml_find_first(root, ".//motifs")
  d <- lapply(1:xml_length(p), function(i) {
    tmp <- xml_attrs(xml_child(p, i))
    data.frame(id = tmp["id"], name = tmp["name"], width = as.numeric(tmp["width"]), stringsAsFactors = FALSE, row.names = tmp["id"])
  })
  motif.info <- do.call(rbind, d)
  #motif.info
  
  # probabilities.
  p <- xml_find_all(root, ".//probabilities")
  length(p) # one per motif.
  motif.prob <- lapply(p, function(n) {
    tmp <- sapply(xml_find_all(n, ".//alphabet_array"), function(l) {
      l <- xml_children(l)
      txt <- xml_text(l)
      # print(txt)
      names(txt) <- xml_attr(l, "letter_id")
      txt
    })
    mode(tmp) <- "numeric"
    colnames(tmp) <- 1:ncol(tmp)
    tmp
  })
  names(motif.prob) <- motif.info$id
  #motif.prob
  
  # scores.
  p <- xml_find_all(root, ".//scores")
  length(p) # one per motif.
  motif.score <- lapply(p, function(n) {
    tmp <- sapply(xml_find_all(n, ".//alphabet_array"), function(l) {
      l <- xml_children(l)
      txt <- xml_text(l)
      # print(txt)
      names(txt) <- xml_attr(l, "letter_id")
      txt
    })
    mode(tmp) <- "numeric"
    colnames(tmp) <- 1:ncol(tmp)
    tmp
  })
  #motif.score
  
  # scanned_sites_summary.
  p <- xml_find_all(root, ".//scanned_sites")
  
  tmp <- do.call(rbind, xml_attrs(p))
  seq.info <-
    data.frame(
      sequence_id = tmp[, "sequence_id"],
      p.value = as.numeric(tmp[, "pvalue"]),
      nsites = as.integer(tmp[, "num_sites"]),
      stringsAsFactors = FALSE
    )
  rownames(seq.info) <- seq.info$sequence_id
  #seq.info
  
  
  motif.site <- lapply(p, function(n) {
    tmp <- sapply(xml_find_all(p[[1]], ".//scanned_site"), function(s) {
      xml_attrs(s)
    })
    data.frame(
      motif_id = tmp["motif_id",],
      position = as.integer(tmp["position",]) + 1,
      width = motif.info[tmp["motif_id",], "width"],
      p.value = as.numeric(tmp["pvalue",]),
      stringsAsFactors = FALSE)
  })
  names(motif.site) <- seq.info$sequence_id
  motif.site <- do.call(rbind, motif.site)
  #motif.site
  
  # range.
  r <- RangedData(
    IRanges(start = motif.site$position, width = motif.site$width),
    motif_name = motif.site$motif_id,
    p.value = motif.site$p.value,
    space = motif.site$sequence_id
  )
  
  # sequenceData.
  sequenceData <- AnnotatedDataFrame(seq.info)
  
  # motifData.
  motifData <- AnnotatedDataFrame(motif.info)
  
  nseq <- nrow(seq.info)
  nmotif <- nrow(motif.info)
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
    models = motif.prob,
    ranges = r
  )
}