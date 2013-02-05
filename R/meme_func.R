setClass("MotifSet",
	representation(
		nmotif = "numeric",
		motif = "list",
		nseq = "numeric",
		sequence = "character"
	)
)

readMEME = function(filename) {
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


readMAST = function(filename, meme) {
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

#foo = readMAST("mast_cby_20.txt")
#foo = readMAST("mast2/mast.txt")
#goo = readMEME("protein/meme2/meme.txt")

getMotifCount = function(M, percentage = FALSE) {
	#x = unlist(lapply(M@motif, nrow))
	x = unlist(lapply(M@motif, function(m) length(unique(m$Id))))
	if (percentage)
		x = 100 * x / M@nseq
	x
}


plotMotifCount = function(M, percentage = FALSE, ...) {
	x = getMotifCount(M, percentage)
	mp = barplot(x, las = 1, axes = FALSE, axisnames = FALSE, ...)
	axis(2, las = 1)
	mtext(1:20, side = 1, at = mp, cex = 0.8)
	box()
}

getMotifMatrix = function(M) {
	r = matrix(0, nrow = M@nseq, ncol = M@nmotif)
	rownames(r) = M@sequence
	colnames(r) = 1:M@nmotif
	for(n in 1:M@nmotif) {
		r[as.character(M@motif[[n]]$Id), n] = 1
	}
	r
}

exportTreedyn = function(x, filename) {
	d = data.frame(Id = rownames(x))
	for(m in 1:ncol(x)) {
		d <- cbind(d, colnames(x)[m], paste("{", x[,m], "}", sep = ""))
	}
	colnames(d) = c("Id", rep(1:ncol(x), each = 2))
	if(!missing(filename))
		write.table(d, file = filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	invisible(d)
}

getMotifArch = function(M, motifs) {
	if (missing(motifs)) motifs = 1:M@nmotif
	M = M@motif
	
	seq = list()
	
	for(m in names(M[motifs])) {
		apply(M[[m]], 1, function(x) {
			if (length(seq[[x[1]]]) == 0)
				seq[[x[1]]] <<- data.frame(Motif = m, Start = as.numeric(x[2]), P = as.numeric(x[3]))
			else {
				d = data.frame(Motif = m, Start = as.numeric(x[2]), P = as.numeric(x[3]))
				seq[[x[1]]] <<- rbind(seq[[x[1]]], d)
			}
		})
	}

	# reorder based on start
	reorderMotifs(seq)
	#seq
}

reorderMotifs = function(M) {
	for(s in names(M)) {
		d = M[[s]]
		#rownames(d) = d[, "Motif"]
		M[[s]] <- d[order(d[, "Start"]), ]
	}
	M
}

countArchs = function(M, min.motif = 3, check.coherence = TRUE, model, remove.rep = TRUE) {
	if (check.coherence && missing(model)) stop("model is needed to check coherence")
	A = list()
	for(s in names(M)) {
		#print(s)
		ms = as.character(M[[s]][, "Motif"])
		if (check.coherence)
			ms = ms[ms %in% model]
		#print(ms)
		if (remove.rep)
			msl = unique(ms)
		else
			msl = ms
		if (length(msl) >= min.motif) {
			doit = TRUE
			if (check.coherence) {
				doit = FALSE
				o = sapply(ms, function(m) which(model %in% m))
				#print(o)
				if (! is.unsorted(o))
					doit = TRUE
				if (! is.unsorted(rev(o)))
					doit = TRUE
			}
			
			if (doit) {
				a = paste(ms, collapse = "-")
				#print(a)
				if(length(A[[a]]) == 0)
					A[[a]] = 1
				else
					A[[a]] = A[[a]] + 1
			}
		}
	}
	d = data.frame(Architecture = names(A), Counts = as.numeric(unlist(A)), Percentage = 100 * as.numeric(unlist(A)) / as.numeric(length(names(M))))
	d = d[order(d[, "Counts"], decreasing = TRUE), ]
	d = data.frame(d, Cumulative = cumsum(d$Percentage))
	d
}

plotCounts = function(x, cut) {
	plot(x$Counts, ylim = c(0, 100), xlab = "Architectures", ylab = "Counts/Percentage", axes = FALSE, type = "l")
	points(x$Counts, col = "black", pch = 21, bg = "gray")
	lines(x$Percentage, col = "darkblue")
	lines(x$Cumulative, col = "darkred")
	points(x$Percentage, col = "darkblue", pch = 21, bg = "steelblue")
	points(x$Cumulative, col = "darkred", pch = 21, bg = "orange")
	if(!missing(cut))
		abline(v = cut, lty = "dotted")
	axis(2, las = 1)
	legend("right", c("Counts", "Percentage", "Cumulative"), pch = 21, col = c("black", "darkblue", "darkred"), pt.bg = c("gray", "steelblue", "orange"), bty = "n")
	box()
}

plotMatrix = function(... , tree, highlight, sep = 2, adj = 10, cex = 0.5, high.col) {
	X = list(...)
	M = lapply(X, getMotifMatrix)
	C = lapply(X, function(x) getMotifCount(x, percentage = TRUE))
	
	if (missing(tree)) {
		# the first matrix is used to construct the tree.
		d = dist(M[[1]])
		h = hclust(d)
		require(ape)
		tree = as.phylo(h)
	}
	
	# reorder matrices according to tree.
	M = lapply(M, function(x) {
		x[tree$tip,]
	})
		
	l = layout(matrix(c(0,1,seq(2,length(M)*2+1)), 2, length(M) + 1), heights = c(1, 10), widths = c(5, rep(5, length(M)*2+1)))

	op = par(mar = c(0.5,0.5,0.5,0))
	plot(tree, cex = 0.6, show.tip.label = FALSE, root.edge = TRUE, use.edge.length = FALSE, yaxs = "i")

	if(!missing(highlight)) {
		l = tree$tip.label
		col = rep("transparent", length(l))
		if(missing(high.col))
			high.col = rainbow(length(highlight))
		if(length(high.col) == 1) {
			highlight = paste(highlight, collapse = "|")
		} else {
			for(k in 1:length(highlight)) {
				col[grepl(highlight[k], l)] = high.col[k]
			}
		}
		tiplabels(pch = 22, col = col, cex = cex, adj = adj, bg = col)
	}
	
	for(k in 1:length(M)) {
		m = M[[k]]
		c = C[[k]]
		par(mar = c(0, 0, 0.5,sep))
		barplot(c, names.arg = NA, space = 0, xaxs = "i", las = 1, ylim = c(0, 105))

		par(mar = c(0.5,0,0.5,sep))
		image(1:ncol(m), 1:nrow(m), t(m), col = c("white", "gray"), axes = FALSE)
		box()
	}
}


# find first tree with species.
findSubTree = function(tree, has.any, has.one) {
	has.any = paste(has.any, collapse = "|")
	has.one = paste(has.one, collapse = "|")
	#print(has.any)
	#print(has.one)
	s = subtrees(tree)
	sel = unlist(lapply(s, function(x) any(grepl(has.any, x$tip.label))))
	if(!(any(sel))) stop("no matches for has.any ids")
	#print(which(sel))
	s2 = s[sel]
	found = FALSE
	sel = NA
	foo = sapply(length(s2):1, function(k) {
		#print(k)
		if(any(grepl(has.one, s2[[k]]$tip.label))) {
			if(!found) {
				#print("Found")
				found <<- TRUE
				sel <<- k
			}
		}
	})
	#print(sel)
	s2[[sel]]
}

getMotifDistribution = function(m) {
	n = sub("(..).*", "\\1", rownames(m))
	r = t(apply(m, 2, function(x) {
		sel = which(x == 1)
		nn = n[sel]
		nn = factor(nn, levels = unique(n))
		table(nn)
	}))

	nt = table(n)
	nc = as.vector(nt)
	names(nc) = names(nt)
	rr = sapply(colnames(r), function(x) {
		100*r[,x]/nc[x]
	})
	rr
}

