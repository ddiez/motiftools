#' Read results from FIMO motif search stored in XML format.
#' 
#' Read results from FIMO motif search stored in XML format and returns a
#' MotifSearchResult object.
#'
#' @param file file XML output file from FIMO.
#' @param description 
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
  
  nseq <- getFimoSequenceInfo(root)$nseq
  motif_info <- getFimoMotifInfo(root)
  nmotif <- nrow(motif_info)
  alpha_bkgd <- getFimoBackground(root)
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
      RangedData(
        IRanges(start = motif_hits$start, width = motif_info[motif_hits$motif_name, "width"]),
        motif_name = motif_hits$motif_name,
        score = motif_hits$score,
        pvalue = motif_hits$pvalue,
        qvalue = motif_hits$qvalue,
        evalue = rep(NA, nrow(motif_hits)),
        sequence = motif_hits$sequence,
        space = motif_hits$seq_id)
  } else {
    motif_hits <- RangedData()
  }
  
  motifData <- AnnotatedDataFrame(motif_info)
  sequenceData <- AnnotatedDataFrame(seq_info)
  
  new(
    "MotifSearchResult",
    info = list(
      tool = "FIMO",
      description = description,
      nseq = nseq,
      nmotif = nmotif
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
  xml_text(xml_find_first(x, "alphabet"))
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

### Legacy code.
## reads XML FIMO output.
readFIMO_old = function(filename, sequenceData, description=NULL) {
  doc = xmlParse(filename)
  top = xmlRoot(doc)
  
  # sequences.
  nseq=as.numeric(xmlGetAttr(top[["sequence-data"]],"num-sequences"))
  
  # motifs.
  motif_info=xmlApply(top, function(m) {
    if(xmlName(m)=="motif") {
      data.frame(motif_name=xmlGetAttr(m,"name"),width=as.numeric(xmlGetAttr(m,"width")),best_f=xmlGetAttr(m,"best-possible-match"),stringsAsFactors = FALSE)
    }
  })
  motif_info=motif_info[!sapply(motif_info,is.null)]
  motif_info=do.call(rbind,motif_info)
  rownames(motif_info)=motif_info$motif_name
  
  nmotif=nrow(motif_info)
  
  # get cisml.
  cismlfile=xmlValue(top[["cisml-file"]])
  free(doc)
  doc=xmlParse(file.path(dirname(filename),cismlfile))
  top=xmlRoot(doc)
  
  all_seq=c()
  res=xmlApply(top, function(p) {
    if(xmlName(p)=="pattern") {
      motif_name=xmlGetAttr(p,"name")
      tmp_seq=xmlApply(p, function(s) {
        if(xmlName(s)=="scanned-sequence") {
          seq_id=xmlGetAttr(s,"accession")
          all_seq<<-c(all_seq,seq_id)
          tmp_match=xmlApply(s, function(m) {
            if(xmlName(m)=="matched-element") {
              qvalue=xmlValue(m[["qvalue"]])
              sequence=xmlValue(m[["sequence"]])
              data.frame(motif_name=motif_name,seq_id=seq_id,pos=as.numeric(xmlGetAttr(m,"start")),score=as.numeric(xmlGetAttr(m,"score")), pvalue=as.numeric(xmlGetAttr(m,"pvalue")),qvalue=as.numeric(qvalue),sequence_hit=sequence,stringsAsFactors = FALSE)
            }
          })
          do.call(rbind,tmp_match)
        }
      })
      do.call(rbind,tmp_seq)
    }
  })
  res=res[!sapply(res,is.null)]
  res=do.call(rbind,res)
  
  # build sequenceData.
  all_seq=sort(unique(all_seq))
  if(missing(sequenceData))
    sequenceData=AnnotatedDataFrame(data.frame(all_seq,sequence_id=all_seq,row.names=1,stringsAsFactors = FALSE))
  # motifData.
  motifData=AnnotatedDataFrame(data.frame(motif_info, row.names="motif_name",stringsAsFactors = FALSE))
  
  if(is.null(res)) {
    ranges=RangedData()
  } else {
    res=res[with(res, order(seq_id,pos)),] # reorder by sequence and then position.
    ranges=RangedData(IRanges(start=res$pos,width=motif_info[res$motif_name,"width"]), motif_name=res$motif_name, score=res$score, pvalue=res$pvalue, qvalue=res$qvalue, evalue=rep(NA, nrow(res)), sequence_hit=res$sequence_hit, space=res$seq_id)
  }
  
  new("MotifSearchResult", info=list(tool="FIMO", description=description, nseq=nseq,nmotif=nmotif), sequences=sequenceData, motifs=motifData, ranges=ranges)
}