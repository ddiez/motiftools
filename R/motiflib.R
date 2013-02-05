setClass("MotifLib", representation = list(
	library = "list",
	annotation = "data.frame"
))

setMethod("[", "MotifLib",
function(x, i, j, drop = TRUE) {
	x@library = x@library[i]
	x@annotation = x@annotation[i,]
	x
})

setMethod("show", "MotifLib",
function(object) {
	message("MotifLib object with: ", length(object@library), " motifs.")
})

setMethod("names", "MotifLib",
function(x) {
		names(x@library)
})
