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


readFIMO = function(file, meme) {
  d = read.table(file, as.is = TRUE)
  tm = unique(d[,1])
  seqs = unique(d[,2])
  
  motifs = list()
  for(m in tm) {
    tmp = d[d[,1] == m,]
    motifs[[as.character(m)]] = data.frame(Id = tmp[,2], Start = as.numeric(tmp[,3]), P = as.numeric(tmp[,6]))
  }
  
  #new("MotifSet", nmotif = length(motifs), motif = motifs, nseq = length(seqs), sequence = seqs)
  new("MotifSet", nmotif = length(motifs), motif = motifs, nseq = meme@nseq, sequence = meme@sequence)
}


readMEMEold = function(filename) {
	con = file(filename, "r")
    lines = readLines(con)
    close(con)
    gms = 0
    gtm = 0
    motifs = list()
    mn = 0
    rts = 0
    rs = 0
    ts = c()
    for(line in lines) {
		if (grepl("^TRAINING SET", line)) { rts = 1; next }
		
		if (rts) {
			if (grepl("^Sequence name", line)) { rs = 1; next }
			if (rs & grepl("^\\*\\*\\*", line)) { rts = 0; rs = 0; next }
			if (rs & grepl("^-------------", line)) next
			if (rs) {
				cks = unlist(strsplit(line, " +"))
				if (length(cks) == 6)
					cks = cks[c(1, 4)]
				else
					cks = cks[1]
				#print(cks)
				ts = c(ts, cks)
			}
		}
		
		if (grepl("^MOTIF", line)) { gms = 1; mn = mn + 1; next }
		if (gms) {
			if (grepl("^Sequence name", line)) { gtm = 1; next }
			if (gtm & grepl("^-----------------------------------------", line)) { gms = 0; gtm = 0; next }
			if (gtm & grepl("^-------------", line)) next
			if (gtm) {
				cks = unlist(strsplit(line, " +"))
				mnc = as.character(mn)
				#print(cks[1:3])
				if (length(motifs[[mnc]]) == 0)
					motifs[[mnc]] = data.frame(Id = cks[1], Start = as.numeric(cks[2]), P = as.numeric(cks[3]))
				else {
					d = data.frame(Id = cks[1], Start = as.numeric(cks[2]), P = as.numeric(cks[3]))
					motifs[[mnc]] = rbind(motifs[[mnc]], d)
				}
			}
		}
    }
    #lines
    #print(paste("Number of sequences:", length(ts)))
    #motifs
    new("MotifSet", nmotif = length(motifs), motif = motifs, nseq = length(ts), sequence = ts)
}

readMEME = function(filename) {
  message("reading XML file ... ",appendLF=FALSE)
  doc = xmlToList(filename)
  message("DONE")
  
  # get sequence ids and number.
  seqset=doc[["training_set"]]
  seqset=seqset[names(seqset)=="sequence"]
  seqset
  seqs=c()
  seqs_id=c()
  for(seq in seqset) {
    seqs=c(seqs, seq["name"])
    seqs_id=c(seqs_id, seq["id"])
  }
  names(seqs)=seqs_id
  nseq=length(seqs)
  
  # model.
  model=doc[["model"]]
  nmotif=as.numeric(model$nmotifs) # get motif number.
  
  # motifs.
  motifset=doc[["motifs"]] # TODO get PSSM?
  motifset=motifset[seq(5,100,5)]
  motifs=c()
  for(m in motifset) {
    motifs=c(motifs,m["id"])
  }
  
  # sequences.
  seqset=doc[["scanned_sites_summary"]]
  seqset=seqset[names(seqset)=="scanned_sites"]
  res=list()
  for(seq in seqset) {
    if(is.list(seq)) { # if no motifs this results in a vector with the values in .attrs
      seq_id=seq$.attrs[["sequence_id"]]
      seq_name=seqs[seq_id]
      l=length(seq)-1
      for(s in seq(1,l,4)) {
        motif_id=seq[[s]]
        pos=seq[[s+2]]
        pvalue=seq[[s+3]]
        res=c(res,list(data.frame(seq_id=seq_name, motif_id=motif_id, pos=as.numeric(pos), pvalue=as.numeric(pvalue))))
      }
    }
  }
  res=do.call(rbind,res)
  
  ## convert to old style (for now!)
  motifs=lapply(unique(res$motif), function(m) {
    tmp=res[res$motif_id==m,]
    data.frame(Id=tmp$seq_id,Start=tmp$pos,P=tmp$pvalue)
  })
  names(motifs)=sub("motif_","",unique(res$motif))
  names(seqs)=NULL
  
  #list(seqs=seqs, nseq=nseq, nmotif=nmotif, motifs=motifs, results=res)
  new("MotifSet", nmotif = nmotif, motif = motifs, nseq = nseq, sequence = seqs)
}

readMAST = function(filename) {
  message("reading XML file ... ",appendLF=FALSE)
  doc = xmlToList(filename)
  message("DONE")
  
  # get motifs.
  motifs=doc[["motifs"]]
  motifs=motifs[names(motifs) == "motif"]
  nmotif=length(motifs)
  motif_info=lapply(motifs, function(x)
    data.frame(motif_id=x["id"],motif_name=x["name"],width=x["width"],consensus=x["best_f"])
  )
  motif_info=do.call(rbind, motif_info)
  rownames(motif_info)=motif_info$motif_id
  
  # get sequences.
  seq=doc[["sequences"]]
  nseq=as.numeric(seq[["database"]]["seq_count"]) # number of sequences.
  
  seq=seq[names(seq)=="sequence"]
  seqs=c() # store all sequence ids.
  res=list()
  for(seg in seq) {
    seq_id=seg$.attrs["name"]
    seqs=c(seqs,seq_id)
    seg=seg[names(seg)=="seg"]
    for(hit in seg) {
      hit=hit[names(hit)=="hit"]
      for(h in hit) {
        res=c(res, list(data.frame(seq_id=seq_id,motif_id=h["motif"],motif_name=motif_info[h["motif"],"motif_name"],pos=as.numeric(h["pos"]),pvalue=as.numeric(h["pvalue"]))))
      }
    }
  }
  res=do.call(rbind, res)

  ## convert to old style (for now!)
  motifs=lapply(motif_info$motif_name, function(m) {
    tmp=res[res$motif_name==m,]
    if(nrow(tmp)==0)
      data.frame()
    else
      data.frame(Id=tmp$seq_id,Start=tmp$pos,P=tmp$pvalue)
  })
  names(motifs)=motif_info$motif_name
  
  names(seqs)=NULL
  
  new("MotifSet", nmotif = nmotif, motif = motifs, nseq = nseq, sequence = seqs)
}

readMASTold = function(filename, meme) {
	con = file(filename, "r")
    lines = readLines(con)
    close(con)

    motifs = list()
    sn = 0
    ts = c()
    rm = 0
    rms = 0
	rs = 0
	cs = ""
	tm = NULL
	tw = NULL
	pval = list()
    for(line in lines) {
		if (grepl("Database contains", line)) {
			sn = as.numeric(gsub(".+Database contains (.+?) .+", "\\1", line))
			#print(sn)
			next
		}

		if (grepl("MOTIF WIDTH", line)) {
			rm = 1
			next
		}

		if (grepl("^SECTION II", line)) {
			rms = 1
			next
		}

		if (rms) {
			if (grepl("^-------------", line)) { rs = 1; next }
		}

		if (rs) {
			if (grepl("^\\*\\*\\*\\*", line)) { rs = 0; rms = 0; next }
			cks = unlist(strsplit(line, " +"))
			if (is.na(cks[1])) { #print("end?");
			next }
			if (cks[1] == "") {
				#print("previous seq")
				#print(cs)
				#print(cks[2])
				motifs[[cs]] = paste(motifs[[cs]], cks[2], collapse = "")
			} else {
				cs = cks[1]
				#print(cks[1])
				#print(cks[2])
				#print(cks[3])
				motifs[[cks[1]]] = cks[3]
				pval[[cks[1]]] = cks[2]
			}
			#print(cks[1])
			#if (length(motifs[[mnc]]) == 0)
            #	motifs[[mnc]] = data.frame(Id = cks[1], Start = as.numeric(cks[2]), P = as.numeric(cks[3]))
            #else {
            #	d = data.frame(Id = cks[1], Start = as.numeric(cks[2]), P = as.numeric(cks[3]))
            #	motifs[[mnc]] = rbind(motifs[[mnc]], d)
            #}
		}

		if (rm) {
			if (grepl("-----", line)) next
			#if (grepl("", line)) { rm = 0; next }
			m = unlist(strsplit(line, " +"))
			if (is.na(m[1])) { #print("end?");
				rm = 0;
				next
			}
			#print(m)
			if (is.null(tm)) tm = m[2]
			else tm = c(tm, m[2])
			if (is.null(tw)) tw = m[3]
			else tw = c(tw, m[3])
		}
		#if (grepl("^SECTION I", line)) { rts = 1; next }
	}
	ts = names(motifs)
	parse_motifs = function(x, p, m, w) {
		#print(x)
		#print(w)
		#print(m)
		ml = list()
		for(n in 1:length(m)) {
			ml[[m[n]]] = data.frame()
		}
		for(n in 1:length(x)) {
			nn = names(x)[n]
			#print(paste(">Doing sequence", nn))
			tmp = gsub(" ", "", x[n])
			tmp = gsub("^\\[", "0_\\[", tmp)
			tmp = gsub("\\]_\\[", "\\]_0_\\[", tmp)
			#print(tmp)
			cks = unlist(strsplit(tmp, "_"))
			sel = grepl("\\[", cks)
			ss = cks[!sel]
			mm = cks[sel]
			mm = gsub("\\[|\\]", "", mm)
			#print(mm)
			if (length(mm) > 0) {
			for(k in 1:length(mm)) {
			#print(paste("  +Found motif: ", mm[k]))
			#print(paste("     -Eval: ", p[k]))
			#print(paste("     -Start:", ss[k]))
				if (is.null(ml[[mm[k]]]))
					ml[[mm[nn]]] = data.frame(Id = nn, Start = ss[k], P = p[n])
				else {
					d = data.frame(Id = nn, Start = ss[k], P = p[n])
					ml[[mm[k]]] = rbind(ml[[mm[k]]], d)
				}
			}
			}
		}
		ml
	}
	motifs = parse_motifs(unlist(motifs), unlist(pval), tm, tw)
	if(missing(meme))
		new("MotifSet", nmotif = length(motifs), motif = motifs, nseq = sn, sequence = ts)
	else
		new("MotifSet", nmotif = length(motifs), motif = motifs, nseq = meme@nseq, sequence = meme@sequence)
	#motifs
}

