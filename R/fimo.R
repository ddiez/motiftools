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
