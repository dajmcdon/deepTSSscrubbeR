
#' Expand Ranges
#'
#' Expand TSS GRanges for downstream analysis
#'
#' @importFrom GenomicRanges GRanges
#'
#' @param deep_obj deep_tss object
#' @param sequence_expansion Number of bases to expand on each side for surrounding sequence analysis
#' @param signal_expansion Number of bases to expand on each side for surrounding signal analysis
#' @param shape_expansion Number of bases to expand on each side for surrounding DNA shape analysis
#'
#' @return Expanded GRanges objects
#'
#' @export

expand_ranges <- function(
	deep_obj,
	sequence_expansion = 10,
	signal_expansion = 10,
	shape_expansion = 10
) {

	## Pull out data.
	original_ranges <- as.data.table(deep_obj@experiment)[,
		.(seqnames, start, end, strand, tss)
	]
	original_ranges <- unique(original_ranges)

	## Expand sequence ranges.
	sequence_expanded <- copy(original_ranges)
	sequence_expanded[,
		c("start", "end") := list(
			start - sequence_expansion,
			end + sequence_expansion
		)
	]
	sequence_expanded[, sequence_index := seq_len(.N)]

	sequence_expanded <- makeGRangesFromDataFrame(
		sequence_expanded, keep.extra.columns = TRUE
	)

	## Expand surrounding signal ranges.
	signal_expanded <- copy(original_ranges)
	signal_expanded[,
		c("start", "end") := list(
			start - signal_expansion,
			end + signal_expansion
		)
	]
	signal_expanded[, signal_index := seq_len(.N)]

	signal_expanded <- makeGRangesFromDataFrame(
		signal_expanded, keep.extra.columns = TRUE
	)

	## Expand ranges to get DA shape.
	shape_expanded <- copy(original_ranges)
	shape_expanded[,
		c("start", "end") := list(
			start - shape_expansion,
			end + shape_expansion
		)
	]
	shape_expanded[, shape_index := seq_len(.N)]

	shape_expanded <- makeGRangesFromDataFrame(
		shape_expanded, keep.extra.columns = TRUE
	)

	## Add expanded ranges back to deep tss object.
	deep_obj@ranges$sequence <- sequence_expanded
	deep_obj@ranges$signal <- signal_expanded
	deep_obj@ranges$shape <- shape_expanded

	deep_obj@settings$sequence_expansion <- sequence_expansion
	deep_obj@settings$signal_expansion <- signal_expansion
	deep_obj@settings$shape_expansion <- shape_expansion

	return(deep_obj)
}
