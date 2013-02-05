### UNIPROT API
getUniprotRelease = function() {
	require(RCurl)
	
	u = "http://www.uniprot.org"
	
	h = basicHeaderGatherer()
	foo = getURI(u, headerfunction = h$update)
	h$value()["X-UniProt-Release"]
}
### mapping tool
# acc: vector of accession numbers.

# test1: acc = c("P13368", "P01112")
# test2: acc = c("P13368", "P01112", "FIFJIDJD")
# test3: acc = c("P13368", "P13368", "P01112")

uniprotLogin = function() {
	require(RCurl)
	
	#url = "http://www.uniprot.org/mapping"

	opts = list()
	opts[["followlocation"]] = TRUE # this is neccesary here.

	getCurlHandle(.opts = opts)
}


mapID2IDuniprot = function(acc, curl, from = "ACC", to = "P_ENTREZGENEID") {
	require(RCurl)
	if(missing(curl)) stop("please, first login with uniprotLogin() and then provide the returned curl!")

	u = "http://www.uniprot.org"
	t = "mapping/"

	url = paste(u, t, sep = "/")


	params = list()
	params[["from"]] = from
	params[["to"]] = to
	params[["format"]] = "tab"
	params[["query"]] = paste(acc, collapse = " ")
	
#	opts = list()
#	opts[["followlocation"]] = TRUE
	
	# for testing:
#	furl = paste(paste(u,t,sep = "/"), paste(paste(names(params), unlist(params), sep = "="), collapse = "&"), sep = "?")
	
	# send form:
	foo = postForm(url, .params = params, curl = curl)
	
	# reads data and returns a data.frame.
	parseResult = function(x) {
		r = unlist(strsplit(strsplit(x, "\n")[[1]][-1], "\t"))
	
		nr = length(r)/2
		d = matrix(NA, nrow = nr, ncol = 2)
		for(i in 1:nr) {
			k = i*2-1
			d[i, ] = c(r[k], r[k+1])
		}
		data.frame(FROM = d[,1], TO = d[,2], stringsAsFactors = FALSE)
	}
	
	res  = parseResult(foo)
	
	# TODO:
	# 1. deal with Uniprot that map to several EG (duplicated Uniprot in res)
	# 2. deal with Uniprot that do not map to EG (not in res)
	
	m = matrix("", nrow = length(acc), ncol = 2)
	m[,1] = acc
	apply(res, 1, function(x) {
		m[m[,1] == x[1], 2] <<- x[2]
	})

	data.frame(FROM = m[,1], TO = m[,2], stringsAsFactors = FALSE)
}

# mapAll2EG: several strategies to map UNIPROBE and other ids into Entrez Gene.
# -----------------------------------------------------------------------------
# x: vector of ids to map.
# hack_file: file containing speciall mappings for entries that do not map well.
mapAll2EGuniprot = function(x, curl, hack_file) {
	m1 = mapID2IDuniprot(x$ACC, curl = curl)
	m2 = mapID2IDuniprot(x$ACC, curl = curl, from = "P_REFSEQ_AC", to = "ACC")
	m3 = mapID2IDuniprot(x$ACC, curl = curl, from = "EMBL", to = "ACC")

	acc = x$ACC

	sel = m2$TO != ""
	acc[sel] = m2$TO[sel]
	sel = m3$TO != ""
	acc[sel] = m3$TO[sel]

	x$ACC_MAPPED = acc

	x$EG = mapID2IDuniprot(acc, curl = curl)[,2]
	
	if (!missing(hack_file)) {
		h = read.table(hack_file, fill = TRUE)

		apply(h, 1, function(z) {
			x[x$ANN %in% z[1], ]$ACC_MAPPED <<- z[2]
		})
		x$EG_MAPPED = mapID2IDuniprot(x$ACC_MAPPED, curl = curl)[,2]
	} else {
		x$EG_MAPPED = x$EG
	}
	
	x
}

