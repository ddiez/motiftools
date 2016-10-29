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