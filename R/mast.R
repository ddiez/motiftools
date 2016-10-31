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
  
  rangeData <- RangedData(
    IRanges(start = hit_info$start, width = motif_info[hit_info$motif_id, "width"]),
    motif_name = hit_info$motif_id,
    score = rep(NA, nrow(hit_info)),
    p_value = hit_info$pvalue,
    space = hit_info$seq_id
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
### Legacy code.

## reads XML MAST output.
readMAST_old = function(filename, sequenceData, description=NULL) {
  doc=xmlParse(filename)
  top=xmlRoot(doc)
  
  # get motifs.
  motif_info=xmlApply(top[["motifs"]],function(m) {
    if(xmlName(m)=="motif") {
      attr=xmlAttrs(m)
      bad=FALSE
      if("bad" %in% names(attr)) bad=TRUE
      data.frame(motif_id=attr[["id"]],motif_name=attr[["name"]],width=as.numeric(attr[["width"]]),best_f=attr[["best_f"]],bad=bad)
    }
  })
  motif_info=do.call(rbind,motif_info)
  rownames(motif_info)=motif_info$motif_id
  
  nmotif=nrow(motif_info)
  nseq=as.numeric(xmlAttrs(top[["sequences"]][["database"]])[["seq_count"]])
  
  # get sequences.
  all_seqs=c()
  res=xmlApply(top[["sequences"]],function(s) {
    if(xmlName(s)=="sequence") {
      seq_id=xmlAttrs(s)[["name"]]
      all_seqs<<-c(all_seqs,seq_id)
      xmlApply(s, function(ss) {
        if(xmlName(ss)=="seg") {
          xmlApply(ss, function(h) {
            if(xmlName(h)=="hit") {
              data.frame(seq_id=seq_id, motif_name=motif_info[xmlAttrs(h)[["motif"]],"motif_name"], pvalue=as.numeric(xmlAttrs(h)[["pvalue"]]),pos=as.numeric(xmlAttrs(h)[["pos"]]))
            }
          })
        }
      })
    }
  })
  
  res=unlist(res,recursive=FALSE,use.names=FALSE)
  res=unlist(res,recursive=FALSE,use.names=FALSE)
  res=res[!sapply(res,is.null)]
  res=do.call(rbind,res)
  
  # build sequenceData.
  all_seqs=sort(unique(all_seqs))
  if(missing(sequenceData))
    sequenceData=AnnotatedDataFrame(data.frame(all_seqs,sequence_id=all_seqs,row.names=1))
  # motifData.
  rownames(motif_info)=motif_info$motif_name
  motifData=AnnotatedDataFrame(data.frame(motif_info, row.names="motif_name"))
  
  if(is.null(res)) {
    ranges=RangedData()
  } else {
    res=res[with(res, order(seq_id,pos)),] # reorder by sequence and then position.
    
    ranges=RangedData(IRanges(start=res$pos,width=motif_info[res$motif_name,"width"]), motif_name=res$motif_name, score=rep(NA,nrow(res)), pvalue=res$pvalue, qvalue=rep(NA, nrow(res)), evalue=rep(NA, nrow(res)), space=res$seq_id)
  }
  
  new("MotifSearchResult", info=list(tool="MAST", description=description, nseq=nseq,nmotif=nmotif), sequences=sequenceData,motifs=motifData,ranges= ranges)
}