readTOMTOM = function(file, do.cut = FALSE, cut.col = "q-value", cutoff = 0.05, filter = "TRANSFAC", do.clean = TRUE) {
  foo = read.table(file, comment.char = "", sep = "\t", header = TRUE, check.names = FALSE, as.is = TRUE)
  colnames(foo)[1] = sub("#", "", colnames(foo)[1])
  
  pwm = unique(c(foo[,1], foo[,2]))
  
  #k = 0
  pwmm = matrix(0, nrow = length(pwm), ncol = length(pwm), dimnames = list(pwm, pwm))
  n <- apply(foo, 1, function(x) {
    if (do.cut) {
      #			print(x[cut.col])
      if (as.numeric(x[cut.col]) <= cutoff) {
        #k <<- k + 1
        pwmm[x[1], x[2]] <<- 1
      }
    } else {
      #k <<- k + 1
      pwmm[x[1], x[2]] <<- 1
    }
  })
  
  print(paste("Accepted comparisons:",length(which(pwmm == 1))))
  if(do.clean) {
    pwmm = cleanMatrix(pwmm)
    print(paste("After cleaning:",length(which(pwmm == 1))))
  }
  
  
  # add default annotations.
  mat.key = matrix(0, nrow = nrow(pwmm), 1, dimnames = list(rownames(pwmm), "Dataset"))
  #mat.key[grep("sci09", rownames(mat.key))] = 1
  #mat.key[grep("cell08", rownames(mat.key))] = 2
  #mat.key[grep("embo10", rownames(mat.key))] = 3
  #mat.key[grep("^MA", rownames(mat.key))] = 4
  
  # colors.
  col.key = c("Default" = "gray")
  #col.key = c("grey", "darkblue", "darkgreen", "yellow2", "darkred")
  #names(col.key) = c("TRANSFAC", "SCI09", "CELL08", "EMBO10", "JASPAR")
  
  # filter?
  if (!missing(filter)) {
    #print(col.key)
    sel.f = which(! names(col.key) %in% filter) - 1
    #print(sel.f)
    sel.m = rownames(mat.key[mat.key[,1] %in% sel.f,,drop = FALSE])
    #print(sel.m)
    mat.key = mat.key[sel.m,,drop = FALSE]
    col.key = col.key[sel.f + 1]
    pwmm = pwmm[rownames(pwmm) %in% sel.m, colnames(pwmm) %in% sel.m]
    n = length(which(pwmm == 1))
    print(paste("Filtering:", n, "comparisons remained after filtering."))
  }
  
  # cluster
  d = clusterMatrix(pwmm)
  
  
  new("TomTom", matrix = pwmm, cutoff.type = cut.col, cutoff = cutoff, matrix_key = mat.key, color_key = col.key, dendrogram = d)
}

## reads XML FIMO output.
readFIMO = function(filename, sequenceData, description=NULL) {
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

#' Read MEME XML file format
#' 
#' Read MEME XML file format
#'
#' @param filename MEME XML file.
#' @param sequenceData sequence data object (optional).
#' @param description dataset description.
#'
#' @return MotifSearchResult object.
#' @export
#'
#' @examples
#' NULL
readMEME <- function(filename, sequenceData, description = NULL) {
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

## reads XML MAST output.
readMAST = function(filename, sequenceData, description=NULL) {
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