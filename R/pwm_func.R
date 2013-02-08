readPWM = function(file, format = "uniprobe") {
	switch(format, uniprobe = readPWMuniprobe(file), jaspar = readPWMjaspar(file), 
		transfac = readPWMtransfac(file))
}


readPWMtransfac = function(file) {
	foo <- read.table(file = file, fill = TRUE, na.string = "", as.is = TRUE)
	sel.id = which(foo[, 1] == "AC")

	foo[sel.id, 2]
}


readPWMtransfac = function(file) {
	table <- read.table(file = file, fill = TRUE, na.string = "")

	table.vector <- cbind(as.vector(table[[1]]), as.vector(table[[2]]), as.vector(table[[3]]), 
		as.vector(table[[4]]), as.vector(table[[5]]))
	XX <- which(table.vector == "XX")
	DE <- which(table.vector == "AC")
	names <- table.vector[DE, 2]
	listPWM <- list()
	for (i in seq(length(DE))) {
		listPWM[[i]] <- matrix(as.numeric(table.vector[(DE[i] + 1):(XX[i] - 1), 
			2:5]), nrow = 4, byrow = T, dimnames = list(c("A", "C", "G", "T")))
		colnames(listPWM[[i]]) <- 1:(length(listPWM[[i]])/4)
		names(listPWM)[i] <- table.vector[DE[i], 2]
	}
	listPWM
}

readPWMuniprobe = function(file) {
	foo = scan(file, what = "", sep = "\n")
	sel.a = grep("A:", foo)
	ids = foo[sel.a - 1]
	a = strsplit(foo[sel.a], "\t")
	c = strsplit(foo[sel.a+1], "\t")
	g = strsplit(foo[sel.a+2], "\t")
	t = strsplit(foo[sel.a+3], "\t")
	listPWM = lapply(1:length(ids), function(k) {
		aa = a[[k]][-1]
		cc = c[[k]][-1]
		gg = g[[k]][-1]
		tt = t[[k]][-1]
		matrix(as.numeric(c(aa,cc,gg,tt)), nrow = 4, byrow = TRUE, dimnames = list(c("A", "C", "G", "T"), 1:length(aa)))
	})
	names(listPWM) = ids
	listPWM
}

readPWMjaspar = function(file, return.pfm = FALSE) {
    foo = read.table(file, sep = "\t", as.is = TRUE)
    ids = foo[, 1]
    path = dirname(file)

    readMatrix = function(x) {
        m = scan(file.path(path, paste(x, ".pfm", sep = "")), quiet = TRUE)
        m = matrix(m, nrow = 4, byrow = TRUE, dimnames = list(c("A", "C", "G", "T"), 1:(length(m)/4)))
        if(return.pfm)
        	m
        else {
			s = colSums(m)
			t(t(m)/s)
		}
    }
    listPWM = lapply(ids, readMatrix)
    names(listPWM) = ids
    listPWM
}

#readPWMjaspar = function(file) {
	#foo = read.table(file, sep = "\t", as.is = TRUE)
	#ids = foo[, 1]
	#path = dirname(file)

	#readMatrix = function(x) {
		#m = read.table(file.path(path, paste(x, ".pfm", sep = "")))
		#colnames(m) = 1:ncol(m)
		#rownames(m) = c("A", "C", "G", "T")
		#s = colSums(m)
		#m/s
	#}
	#pwm = sapply(ids, readMatrix)
	#names(pwm) = ids
	#pwm
#}

writePWMuniprobe = function(listPWM, file, ann = NULL) {
	file.create(file)
	foo = lapply(names(listPWM), function(pwm) {
		if(!is.null(ann)) s = paste(pwm, ann[pwm], sep = "\t") else s = pwm
		write(s, file = file, append = TRUE)
		m = listPWM[[pwm]]
		sapply(rownames(m), function(l) {
			l = paste(paste(l, ":", sep = ""), paste(m[l, ], collapse = "\t"), 
				sep = "\t")
			write(l, file, append = TRUE)
		})
	})
}

# this can be used with HOMER.
writePWMtable = function (listPWM, file, ann = NULL) {
  file.create(file)
  foo = lapply(names(listPWM), function(pwm) {
    s = paste(">", pwm, sep = "")
    if (!is.null(ann)) 
      s = paste(s, paste(pwm, ann[pwm], sep = "_"), 0, sep = "\t")
    else
      s = paste(s, pwm, 0, sep = "\t")
    write(s, file = file, append = TRUE)
    m = listPWM[[pwm]]
    write(t(m), file, append = TRUE, sep = "\t")
  })
}
