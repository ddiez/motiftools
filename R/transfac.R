require(XML)

# matrix: will be mapped to factors.
# factor: will be mapped to gene.
# gene: will be mapped to entrez gene.

transfacLogin = function(username, password) {
	require(RCurl)
	require(RHTMLForms)
	
	# log in into gene-regulation.
	url = "http://www.gene-regulation.com/login"

	opts = list()
#	opts[["followlocation"]] = TRUE # this seems not be necessary here.
#	opts[["cookiejar"]] = "cookies.txt" # either this or "cookiefile" would do the trick.
	opts[["cookiefile"]] = "cookies.txt"

	# get login form and login.
	# maybe it can be done directly without need of RHTMLForms (and hence XML),
	# but I will try that later.
	txt = getURLContent(url)
	doc = htmlParse(txt, asText = TRUE)
	fun = getHTMLFormDescription(doc)
	fun1 = createFunction(fun[[1]])
	curl = getCurlHandle(.opts = opts)
	fun1(password, username, .url = url, .curl = curl)
	curl
}

mapID2IDtransfac = function(ids, curl, type = "matrix") {
	if(missing(curl)) stop("please, first login with transfacLogin() and then provide the returned curl!")
	
	type = match.arg(type, c("matrix", "factor", "gene"))
	
	# query.
	r = sapply(ids, function(x) {
		#print(x)
		params = list()
		params[["AC"]] = x
		url = "http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi/"
		res = postForm(url, .params = params, curl = curl, style = "POST")

		rr = switch(type,
			"matrix" = parseMatrix(res),
			"factor" = parseFactor(res),
			"gene" = parseGene(res)
		)
		#print(rr)
	})
	outPut(r)
}

outPut = function(x) {
	res = unlist(x)
	names(res) = rep(names(x), sapply(x, length))
	res
}

parseMatrix = function(r) {
	h = htmlParse(r)
	res = unlist(xpathApply(h, "/html//a", fun = xmlValue))
	sel <- res == "BF"
	if (any(sel))
		res[which(sel) + 1]
	else ""
}

parseFactor = function(r) {
	h = htmlParse(r)
	res = unlist(xpathApply(h, "/html//a", fun = xmlValue))
	sel <- res == "GE"
	if (any(sel))
		res[which(sel) + 1]
	else ""
}

parseGene = function(r) {
	h = htmlParse(r)
	res = unlist(xpathApply(h, "/html//a", fun = xmlValue))
	if (length(res) == 1) {
		if (res == "TRANSFAC\nProfessional") "NOTFOUND"
	} else {
		res = strsplit(xpathApply(h, "/html", xmlValue)[[1]], "\n")[[1]]
		# TODO: check for the result! maybe no ENTREZGENE entry.
		sel = grepl("ENTREZGENE:", res)
		if (any(sel)) {
			res = res[sel]
			sub(".$", "", sub(".*ENTREZGENE: ", "", res))
		} else ""
	}
}

mapAll2EGtransfac = function(x, curl = curl) {
	res1 = mapID2IDtransfac(x$ID, curl = curl, type = "matrix")
	res2 = mapID2IDtransfac(res1, curl = curl, type = "factor")
	res3 = mapID2IDtransfac(res2, curl = curl, type = "gene")
	
	ids = names(res1)

	src = sapply(names(res1), function(z) x$SOURCE[x$ID == z])
	ann_ori = sapply(names(res1), function(z) x$ANN[x$ID == z])

	d = data.frame(ID = names(res1), SOURCE = src, ANN_ORI = ann_ori, ANN = ann_ori, ACC_ORI = res1, ACC = res1, ACC_MAPPED = res2, EG = res3, EG_MAPPED = res3, stringsAsFactors = FALSE)
	d$EG_MAPPED[d$EG_MAPPED == "NOTFOUND"] = ""
	d
}



