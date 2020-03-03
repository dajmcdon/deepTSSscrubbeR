
#' Retrieve DNA Shape
#'
#' Use DNAshapeR to get shape information on signal surrounding TSSs
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom reticulate array_reshape
#' @importFrom DNAshapeR getShape encodeSeqShape
#' @importFrom Biostrings writeXStringSet DNAStringSet
#'
#' @param deep_obj deep_tss object
#' @param shape_features Shape features to consider
#'
#' @rdname encode_shape-function
#'
#' @export

encode_shape <- function(
	deep_obj,
	shape_features = c(
		"MGW", "ProT", "Roll", "HelT", "Rise",
		"Shift", "Slide", "Tilt", "Buckle", "Opening",
		"Shear", "Stagger", "Stretch", "EP"
	)
) {

	## Grab shape sequences from object
	shape_parameters <- deep_obj@ranges$shape

	## Write temporary fasta file required for DNAShapeR
	shape_parameters$shape_seq %>%
		DNAStringSet %>%
		writeXStringSet("temp.fa")

	## Get encoded DNA shape features.
	shape <- getShape("temp.fa", shapeType = shape_features) %>%
		encodeSeqShape("temp.fa", ., featureNames = paste0("1-", shape_features)) %>%
		array_reshape(dim = c(nrow(.), 1, ncol(.), 1))

	deep_obj@encoded$shape <- shape
	return(deep_obj)
}
