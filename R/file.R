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
readFIMO = function(filename, description=NULL) {
  doc = xmlParse(filename)
  top = xmlRoot(doc)
  
  # sequences.
  nseq=as.numeric(xmlGetAttr(top[["sequence-data"]],"num-sequences"))
  
  # motifs.
  motif_info=xmlApply(top, function(m) {
    if(xmlName(m)=="motif") {
      data.frame(motif_name=xmlGetAttr(m,"name"),width=as.numeric(xmlGetAttr(m,"width")),best_f=xmlGetAttr(m,"best-possible-match"))
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
              data.frame(motif_name=motif_name,seq_id=seq_id,pos=as.numeric(xmlGetAttr(m,"start")),score=as.numeric(xmlGetAttr(m,"score")), pvalue=as.numeric(xmlGetAttr(m,"pvalue")),qvalue=as.numeric(qvalue),sequence_hit=sequence)
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
  
  res=res[with(res, order(seq_id,pos)),] # reorder by sequence and then position.
  
  ranges=RangedData(IRanges(start=res$pos,width=motif_info[res$motif_name,"width"]), motif_name=res$motif_name, score=res$score, pvalue=res$pvalue, qvalue=res$qvalue, evalue=rep(NA, nrow(res)), sequence_hit=res$sequence_hit, space=res$seq_id)

  new("MotifSearchResult", info=list(tool="FIMO", description=description, nseq=nseq,nmotif=nmotif,motif_info=motif_info,sequence_info=sort(unique(all_seq))), ranges= ranges)
}

## reads XML MEME output.
readMEME = function(filename, description=NULL, return.old=FALSE) {
  doc = xmlParse(filename)
  top = xmlRoot(doc)
  
  # get sequence ids and number:
  seqset=xmlApply(top[["training_set"]], function(s) {
    if(xmlName(s)=="sequence") {
      data.frame(seq_id=xmlGetAttr(s, "id"),seq_name=xmlGetAttr(s,"name"))
    }
  })
  seqset=seqset[! sapply(seqset,is.null)]
  seqset=do.call(rbind,seqset)
  rownames(seqset)=seqset$seq_id
  nseq=nrow(seqset)
  
  # motifs.
  nmotif=as.numeric(xmlValue(top[["model"]][["nmotifs"]]))
  
  motif_info=xmlApply(top[["motifs"]], function(m) {
    if(xmlName(m)=="motif")
      data.frame(motif_id=xmlGetAttr(m,"id"),motif_name=xmlGetAttr(m,"name"),width=as.numeric(xmlGetAttr(m,"width")))
  })
  motif_info=do.call(rbind,motif_info)
  rownames(motif_info)=motif_info$motif_id
  
  # sequence hits:
  res=xmlApply(top[["scanned_sites_summary"]], function(s) {
    if(xmlName(s)=="scanned_sites") {
      seq_id=xmlGetAttr(s,"sequence_id")
      tmp=xmlApply(s, function(ss) {
        if(xmlName(ss)=="scanned_site") {
          data.frame(seq_id=seqset[seq_id,"seq_name"],motif_name=motif_info[xmlGetAttr(ss,"motif_id"),"motif_name"],pos=as.numeric(xmlGetAttr(ss,"position")),pvalue=as.numeric(xmlGetAttr(ss,"pvalue")))
        }
      })
      do.call(rbind,tmp)
    }
  })
  res=do.call(rbind,res)
  
  rownames(motif_info)=motif_info$motif_name
  res=res[with(res, order(seq_id,pos)),] # reorder by sequence and then position.
  
  ranges=RangedData(IRanges(start=res$pos,width=motif_info[res$motif_name,"width"]), motif_name=res$motif_name, score=rep(NA,nrow(res)), pvalue=res$pvalue, qvalue=rep(NA, nrow(res)), evalue=rep(NA, nrow(res)), space=res$seq_id)
  
  new("MotifSearchResult", info=list(tool="MEME", description=description, nseq=nseq,nmotif=nmotif,motif_info=motif_info,sequence_info=sort(seqset$seq_name)), ranges= ranges)
}

## reads XML MAST output.
readMAST = function(filename, description=NULL, return.old=FALSE) {
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
  
  rownames(motif_info)=motif_info$motif_name
  res=res[with(res, order(seq_id,pos)),] # reorder by sequence and then position.
  
  ranges=RangedData(IRanges(start=res$pos,width=motif_info[res$motif_name,"width"]), motif_name=res$motif_name, score=rep(NA,nrow(res)), pvalue=res$pvalue, qvalue=rep(NA, nrow(res)), evalue=rep(NA, nrow(res)), space=res$seq_id)
  
  new("MotifSearchResult", info=list(tool="MAST", description=description, nseq=nseq,nmotif=nmotif,motif_info=motif_info,sequence_info=sort(all_seqs)), ranges= ranges)
}