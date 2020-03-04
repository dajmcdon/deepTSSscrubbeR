
#' Get Softclipped Bases
#'
#' Extract softclipped bases from bam
#'
#' @include deep_tss.R
#'
#' @param deep_obj deep_tss object
#'
#' @import tibble
#' @importFrom stringr str_extract str_replace str_sub
#'
#' @rdname get_softclip-function
#'
#' @export

get_softclip <- function(deep_obj) {

	## Get soft-clipped bases.
	soft <- as.data.table(deep_obj@experiment)

	soft[, softclipped := ifelse(
			strand == "+", str_extract(cigar.first, "^\\d+S"),
			str_extract(cigar.first, "\\dS$")
		)
	][,
		softclipped := as.numeric(str_replace(softclipped, "S", "")), NULL
	][,
		softclipped := ifelse(is.na(softclipped), 0, softclipped)
	][,
		soft_bases := ifelse(softclipped == 0, NA, str_sub(seq_firstinread, end = softclipped))
	][,
		c("softclipped", "seq_firstinread") := NULL
	][,
		soft_group := .GRP,
		by = .(seqnames, start, end, strand, soft_bases)
	]

	soft_unique <- soft[, .(seqnames, start, end, strand, soft_bases, soft_group)]
	soft_unique <- unique(soft_unique)
	
	## Add soft-clipped index back to original data.
	soft[, soft_bases := NULL]
	deep_obj@experiment <- as.data.frame(soft)

	## Add soft-clipped ranges to deep tss object.
	soft_unique <- makeGRangesFromDataFrame(soft_unique, keep.extra.columns = TRUE)
	deep_obj@ranges$soft <- soft_unique

	return(deep_obj)
}
