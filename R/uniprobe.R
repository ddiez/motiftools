### UNIPROBE API

uniprobeLogin = function(opts = list()) {
	require(RCurl)

	if (is.null(opts$followlocation))
		opts[["followlocation"]] = TRUE # this is neccesary here.

	getCurlHandle(.opts = opts)
}


# get basic information.
getIdsUniprobe = function(curl, filter) {
	require(RCurl)
	require(XML)
	
	url = "http://the_brain.bwh.harvard.edu/uniprobe/browse.php"
	
	urlc = getURLContent(url)
	xml = htmlParse(urlc, asText = TRUE)
	txt = xpathApply(xml, "/html//td", fun = xmlValue)
	txt = txt[1:(grep("Text Search", txt)-1)]
	
	n = length(txt)
	tf.name = unlist(txt[seq(1, n, 8)])
	tf.id = unlist(txt[seq(2, n, 8)])
	tf.family = unlist(txt[seq(3, n, 8)])
	tf.organism = unlist(txt[seq(4, n, 8)])
	tf.description = unlist(txt[seq(5, n, 8)])
	#tf.ref = unlist(txt[seq(6, n, 8)])
	
	txt2 = xpathApply(xml, "/html//td/*", fun = xmlAttrs)
	sel.zip = which(unlist(lapply(txt2, function(x) any(grepl(".zip", x)))))
	sel.link = sel.zip + 1
	tf.zip = unlist(lapply(txt2[sel.zip], function(z) z["href"]), use.names = FALSE)
	tf.link = unlist(lapply(txt2[sel.link], function(z) z["href"]), use.names = FALSE)

	tf.source = sub("/.*", "", sub("downloads/", "", tf.zip))

	#d = data.frame(ID = tf.id, NAME = tf.name, FAMILY = tf.family, ORG = tf.organism, REF = tf.ref, ZIP = tf.zip, LINK = tf.link, SOURCE = tf.source, stringsAsFactors = FALSE)
	d = data.frame(ID = tf.id, NAME = tf.name, FAMILY = tf.family, ORG = tf.organism, LINK = tf.link, SOURCE = tf.source, stringsAsFactors = FALSE)
	
	if (!missing(filter))
		d = d[d$SOURCE %in% filter,]
	d
}

mapSpecies = function(x) {
	map = c(
		"Saccharomyces cerevisiae" = "yeast",
		"Homo sapiens" = "human",
		"Mus musculus" = "mouse",
		"Rattus norvegicus" = "rat",
		"Caenorhabditis elegans" = "worm"
	)
	map[x]
}


# get details.
getMetaDataUniprobe = function(x, curl) {
	d = x
	curl = getCurlHandle()
	base.url = "http://the_brain.bwh.harvard.edu/uniprobe"
	d = sapply((1:nrow(d)), function(k) {
		l = d$LINK[k]
		view.url = paste(base.url, l, sep = "/")
		#print(view.url)
		
		foo.view = getURLContent(view.url, curl = curl)
		xml = htmlParse(foo.view, asText = TRUE)
		txt = xpathApply(xml, "/html//td", fun = xmlValue)

		getField = function(v, t) {
			v = paste("^\n\t  ", v, "\n", sep = "")
			#print(v)
			sel = grep(v, t)
			if(any(sel)) {
				f = paste(v, "|\n|\t", sep = "")
				#print(f)
				gsub(" ", "", gsub(f, "", t[[sel]]))
			} else
				NA
		}
		
		getMultiField = function(v, t) {
			res = lapply(v, function(vv) {
				#print(t)
				getField(vv, t)
			})
			res = unlist(res, use.names = FALSE)
			sel = is.na(res)
			if(all(sel))
				stop(paste("ERROR: [", v, "] NOT FOUND!", sep = ""))
			res[!sel]
		}
				
		print("--------")
		print(d$NAME[k])
		uniprobe = getMultiField("UniPROBE Accession Number", txt)
		print(uniprobe)
		prot = getMultiField(c("Protein", "Gene"), txt)
		print(prot)
		uniprot = getMultiField(c("Swiss-Prot", "Uniprot"), txt)
		print(uniprot)
		refseq = getMultiField(c("RefSeq", "NCBI RefSeq", "NCBI Acc\\."), txt)
		print(refseq)
		#species = getField("Species", txt)
		#species = mapSpecies(species)
		#domain = getField("Domain", txt)
		
		txt = xpathApply(xml, "/html//a", fun = xmlAttrs)
		sel = grepl("\\.pwm|_pwm\\.txt|_pwm_primary\\.txt", txt)
		if(any(sel))
			pwm_link = file.path(base.url, txt[sel])[1]
		else
			pwm_link = NA
		
		print(pwm_link)
		c(ID = uniprobe, SYMBOL = prot, UNIPROT = uniprot, REFSEQ = refseq, LINK = pwm_link)
		#c(ID = uniprobe, SYMBOL = prot, ACC = uniprot, ORG = species, DOMAIN = domain, LINK = pwm_link)
	})
	data.frame(t(d), stringsAsFactors = FALSE)
}

# get PWMs.
getMatrixUniprobe = function(x, file) {
	res = lapply(x$LINK, function(k) {
		tmp = tempfile()
		download.file(k, tmp)
		readPWMuniprobe(tmp)[[1]]
	})
	names(res) = x$ID
	res
}

getDataUniprobe = function(curl, filter) {
	tmp1 = getIdsUniprobe(curl, filter = filter)
	tmp2 = getMetaDataUniprobe(tmp1, curl)
	data.frame(ID = tmp2$ID, SYMBOL = tmp2$SYMBOL, UNIPROT = tmp2$UNIPROT, REFSEQ = tmp2$REFSEQ, FAMILY = tmp1$FAMILY, ORG = tmp1$ORG, SOURCE = tmp1$SOURCE, LINK = tmp2$LINK, stringsAsFactors = FALSE)
}

