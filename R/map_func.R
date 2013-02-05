
readIdsTransfac = function(file, source = "transfac") {
	foo <- read.table(file = file, fill = TRUE, na.string = "", as.is = TRUE)

	sel.id = which(foo[,1] == "AC")
	sel.ann = sel.id+8
	
	id = foo[sel.id, 2]
	ann = sapply(foo[sel.ann, 2], function(x) if(is.na(x)) "" else x)
	
	data.frame(ID = foo[sel.id,2], SOURCE = source,  ANN_ORI = ann, ANN = ann, stringsAsFactors = FALSE)
}

# uniprobe ids use Symbols. Moreover, each dataset is focused on an organism. Therefore we reverse map the symbols
# to Entrez Gene using the mouse annotation package in Bioconductor.
readIdsUniprobe = function(file, source = "unknown") {
	foo = read.table(file, fill = TRUE, as.is = TRUE)
	sel.id = which(foo[,1] == "A:") - 1
	ids = foo[sel.id,1]
	# remove secondary ones:
	ids = ids[! grepl("_secondary", ids)]
	ann_ori = sub("_.*", "", ids)
	ann = ann_ori
	
	# hack for some ids.
	# Dobox4, Dobox5 and IRC900814 TFs cannot be identified.
	ann[ann == "E2F2"] = "E2f2"
	ann[ann == "E2F3"] = "E2f3"
#	ann[ann == "Zscan4"] = "Zscan4c" # this? in mouse there are several pseudogenes and then another gene call "Zscan4c". Is it the same?
	
	data.frame(ID = ids, SOURCE = source, ANN_ORI = ann_ori, ANN = ann, ACC_ORI = ann, ACC = ann, ACC_MAPPED = ann, stringsAsFactors = FALSE)
}

# JASPAR matrices are generally mapped into UNIPROT ids, although some cases where RefSeq and EMBL ids were used
#  can be found. The strategy is use the UNIPROT API for id mapping to map those IDs into Entrez Gene Ids.
#
# readIdsJaspar: read ids from matrix_list.txt
# --------------------------------------------
# file: file name
#
readIdsJaspar = function(file, source = "jaspar") {
	foo = read.table(file, sep = "\t", as.is = TRUE)
	ids = foo[,1]
	sym = foo[,3] # symbol.
	ann = foo[,5] # annotations.
	
	ann = strsplit(ann, "; ")
	names(ann) = ids
	
	# get accession number.
	acc = lapply(ann, function(x) grep("^acc ", x, value = TRUE))
	acc = lapply(acc, function(x) gsub(" ", "", sub("acc ", "", x)))
	acc[acc == "-"] = ""

	# get family names.
	fam = lapply(ann, function(x) grep("^family ", x, value = TRUE))
	fam = lapply(fam, function(x) gsub(" ", "", sub("family ", "", x)))
	fam[fam == "-"] = ""
	
	d = data.frame(ID = ids, SOURCE = source, ANN = sym, ACC = unlist(acc), FAMILY = unlist(fam), stringsAsFactors = FALSE)
	
	# manage duplicates (complexes):
	sel.d = grepl(",", d$ACC)
	sel = sel.d
	tmp1 = d[!sel,]
	tmp2 = d[sel,]
	
	# hack for AP1
	tmp2$ANN[tmp2$ANN == "AP1"] = "JUN::FOS"
	
	tmp3 = data.frame(ID = rep(tmp2$ID, each = 2), SOURCE = rep(tmp2$SOURCE, each = 2), ANN_ORI = rep(tmp2$ANN, each = 2), ANN = unlist(strsplit(tmp2$ANN, "::")), ACC_ORI = rep(tmp2$ACC, each = 2), ACC = unlist(strsplit(tmp2$ACC, ",")), FAMILY = rep(tmp2$FAMILY, each = 2), stringsAsFactors = FALSE)
	
	# hack for AP1
	tmp3$ANN_ORI[tmp3$ANN_ORI == "JUN::FOS"] = c("AP1", "AP1")
	
	tmp4 = data.frame(ID = tmp1$ID, SOURCE = tmp1$SOURCE, ANN_ORI = tmp1$ANN, ANN = tmp1$ANN, ACC_ORI = tmp1$ACC, ACC = tmp1$ACC, FAMILY = tmp1$FAMILY, stringsAsFactors = FALSE)

	
	data.frame(rbind(tmp3, tmp4), stringsAsFactors = FALSE)
}

plotRatio = function(x, main = "", line = 0, col = c("white", "gray")) {
	sel = x$EG_MAPPED != ""
	tot = length(unique(x$ID))
	map = length(unique(x$ID[sel]))
	
	p_map = 100*map/tot
	p_fail = 100*(tot-map)/tot
	
	l_map = paste(paste("MAP (", map, ")", sep = ""), paste(format(p_map,digits = 3), "%", sep = ""), sep = "\n")
	l_fail = paste(paste("FAIL (", tot-map, ")", sep = ""), paste(format(p_fail,digits = 3), "%", sep = ""), sep = "\n")
	
	op = par(mar = rep(2,4))
	pie(c(map, tot - map), main = "", labels = c(l_map, l_fail), col = col)
	title(main, line = line)
	par(op)
}

mapAll2EGuniprobe = function(up) {
	# check organism and choose according to it.
	require(org.Mm.eg.db)
	
	eg = mget(up$ACC_MAPPED, org.Mm.egALIAS2EG, ifnotfound = NA)
	# manage duplicates: collapse:
	eg2 = sapply(eg, function(x) paste(x, collapse = ","))
	eg2 = sapply(eg2, function(x) if(x == "NA") "" else x)
	up$EG = eg2
	up$EG_MAPPED = eg2
	up
}
