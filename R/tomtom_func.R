setOldClass("dendrogram")
setClass("TomTom", representation(
	matrix = "matrix",
	cutoff.type = "character",
	cutoff = "numeric",
	matrix_key = "matrix",
	color_key = "character",
	dendrogram = "dendrogram"
))
setGeneric("getPWMnames", function(object) standardGeneric("getPWMnames"))
setMethod("getPWMnames", "TomTom",
function(object) {
	colnames(object@matrix)
})


cleanMatrix = function(m) {
	#k = 0
	for(i in colnames(m)) {
		for(j in rownames(m)) {
			if(m[i,j] != 0 | m[j,i] != 0) {
				if(m[i,j] != m[j,i]) {
					#message("inconsistency between ", i, " and ", j, " similarities!")
					#k = k + 1
					m[i,j] = 0
					m[j,i] = 0
				}
			}
		}
	}
	m
}


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


clusterMatrix = function(m, dist.method = "pearson", hclust.method = "complete") {
	dcor <- function(x, method = "pearson", use = "everything") {
    	d = 1 - cor(x, method = method, use = use)
    	#d[d < 0] = 0 # remove effect of negative correlation
    	as.dist(d)
	}
	if (dist.method %in% c("pearson", "spearman"))
		d = dcor(m, method = dist.method)
	else
		d = dist(m, method = dist.method)
		
	h = hclust(d, method = hclust.method)
	as.dendrogram(h)
}

cutplot.dendrogram = function(x, h, cluscol, leaflab= "none", horiz=FALSE, lwd=1,
                              ...)
#
# Name: cutplot.dendrogram
# Desc: takes a dendrogram as described in library(mva), cuts it at level h,
#       and plots the dendrogram with the resulting subtrees in different 
#       colors
# Auth: obviously based on plot.dendrogram in library(mva)
#       modifications by Alexander.Ploner@meb.ki.se  211203
#
# Chng: 050204 AP 
#       changed environment(plot.hclust) to environment(as.dendrogram) to
#       make it work with R 1.8.1
#       250304 AP added RainbowPastel() to make it consistent with picketplot
#       030306 AP slightly more elegant access of plotNode
#
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }

    # Not nice, but necessary
    pn  = stats:::plotNode

    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)

    x = cut(x, h)
    plot(x[[1]], leaflab="none", ...)

    x = x[[2]]
    K = length(x)
    if (missing(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE,
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=horiz)
        x1 = x2 + 1
   }


}

##
## plotTomTom()
##
# dedro: (missing); if a dendrogram is provided, a tree and optionally clustering will be plotted.
# h: (missing); sets the height for cutting the tree.
# tree.size: 2; controls the size of the top tree (plot.dendrogram)
# colkey.size: 1; controls the size of the right and bottom color keys (image)
# cluster.box: FALSE; if TRUE, a box is draw around the cluster insted of a colored cluster.
plotTomTom = function(t, h, tree.size = ifelse(missing(zoom), 2, 4), colkey.size = .5, cluster.box = TRUE, col.c, bottom.colkey = TRUE, zoom, col = c("white", "gray70"), plot.labels = ifelse(missing(zoom), FALSE, TRUE), useRaster = TRUE) {
	m = t@matrix

	mat.key = t@matrix_key
	col.key = t@color_key
	
	dendro = t@dendrogram
	
	d.cut = NULL
	# key
	#if (missing(col.key)) col.key = c("white", "black")
	if (! is.null(mat.key) & is.null(col.key)) {
		col.key = rainbow(length(table(mat.key)))
		names(col.key) = paste("Class", 1:length(table(mat.key)))
	}

	if (! missing(dendro)) {
		d.ord = order.dendrogram(dendro)
		m = m[d.ord, d.ord]
		
		if(!is.null(mat.key))
			mat.key = mat.key[d.ord,,drop = FALSE]
		
		if (!missing(h)) {
			d.cut = cut(dendro, h = h)
	
			# random colors.
			cc = grep("gray|grey|light|white", colors(), value = TRUE, invert = TRUE)
			if(missing(col.c))
				col.c = sample(cc, length(d.cut$lower))
				
			if (! missing(zoom)) {
				dendro = d.cut$lower[[zoom]]
				col.c = col.c[zoom]
				z = zoomDendro(m, d.cut, zoom)
				m = m[z$rx, z$ry]
				mat.key = mat.key[z$rx,, drop = FALSE]
				d.cut = dendro
			}
		}
#
	}
	
	if(missing(col.c))
		col.c = "black" # by default, no cluster colors.	
	
	# layout plot.
	if (!bottom.colkey)
		l = layout(matrix(c(0,1,2,3), 2, 2, byrow = TRUE),
			heights = c(tree.size, 10),
			widths = c(colkey.size, 10))
	else
		l = layout(matrix(c(0,1,2,3,0,4), 3, 2, byrow = TRUE),
			heights = c(tree.size, 10, colkey.size),
			widths = c(colkey.size, 10))
	

	# plot tree.
	op = par(mar = c(ifelse(plot.labels, 8.5,.5),.5,.5,.5))
	if(! missing(dendro)) {
		if (missing(h) | !missing(zoom))
			plot(dendro, leaflab = ifelse(plot.labels, "perpendicular", "none"), xaxs = "i", yaxt = "n", edgePar = list(col = col.c))
		else {
			cutplot.dendrogram(dendro, h, xaxs = "i", yaxt = "n", cluscol = col.c)
			abline(h = h, col = "red", lty = "dotted")
		}
		axis(2, line = 0.3, padj = 0.7)
	} else
		plot(0, axes = FALSE)
	
	
	# plot left key.
	if (!is.null(mat.key)) {
		par(mar = c(.5, .5, 0, 0))
		image(1, 1:nrow(m), t(mat.key), useRaster = useRaster, axes = FALSE, col = col.key)
	} else {
		plot(0, axes = FALSE)
	}
	
	# plot main heatmap.
	par(mar = c(.5,.5,0,.5))
	

	if (! missing(dendro) & !cluster.box & !is.null(d.cut) & missing(zoom)) {
		m = assignCluster(m, d.cut)
		col = c("white", col.c)
	}
			
	image(1:ncol(m), 1:nrow(m), m, col = col, useRaster = useRaster, axes = FALSE, ylab = "", xlab = "")
	
	if (! missing(dendro) & !missing(h) & cluster.box & missing(zoom))
		drawClusterBox(m, d.cut, col.c)
	
	if(!is.null(mat.key))
		legend("topleft", legend = names(col.key), fill = col.key, bty = "n")
	
	legend("bottomright", legend = paste(t@cutoff.type, "=", t@cutoff), bty = "n")
	box()
	
	# plot bottom key.
	if (bottom.colkey) {
		if (!missing(mat.key)) {
			par(mar = c(.5,.5,0,.5))
			image(1:nrow(m), 1, mat.key, useRaster = useRaster, axes = FALSE, col = col.key, xlab = "")
		} else {
			plot(0, axes = FALSE)
		}
	}
			
	par(op)
	
	invisible(col.c)
}

zoomDendro = function(m, dc, z) {
	l = labels(dc$lower[[z]])

	rx = range(which(colnames(m) %in% l))
	ry = range(which(rownames(m) %in% l))
	
	list(rx = rx[1]:rx[2], ry = ry[1]:ry[2])
}

assignCluster = function(m, dc) {
	m2 = matrix(0, nrow = nrow(m), ncol = ncol(m), dimnames = list(colnames(m), rownames(m)))
	aCluster = function(x) {
		l = labels(dc$lower[[x]])
		for(i in l) {
			for(j in l) {
				m2[i, j] <<- x
			}
		}
	}
	sapply(1:length(dc$lower), aCluster)
	m2
}

drawClusterBox = function(m, dc, col) {
	dCluster = function(x) {
		l = labels(dc$lower[[x]])
		
		rx = range(which(colnames(m) %in% l))
		ry = range(which(rownames(m) %in% l))
		rect(rx[1], ry[1], rx[2], ry[2], border = col[x], lwd = 2)
	}
	
	sapply(1:length(dc$lower), dCluster)
}


plotBestCut = function(d) {
	H = seq(0, 1.1, 0.1)
	N = rep(NA, length(H))
	names(N) = H
	for(h in H) {
		print(h)
		dc = cut(d, h)
		N[as.character(h)] = length(dc$lower)
	}
	plot(H, N, las = 1, xlab = "Height", ylab = "# Clusters", pch = 21, bg = "gray", main = "Distribution of # clusters\nvs. height in 'cut'")
}

getGenesFromCluster = function(query, cl, ml, what = "EG") {
	what = match.arg(what, c("EG", "ENS"))
	m = unlist(cl[sapply(cl, function(c) any(c %in% query))])
	g = ml@annotation[m, what]
	g = g[g != ""]
	unlist(strsplit(g, ","))
}

getMotifsFromCluster = function(query, cl) {
	unlist(cl[sapply(cl, function(c) any(c %in% query))])
}

