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