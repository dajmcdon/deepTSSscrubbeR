
#' Expand Ranges
#'
#' Expand TSS GRanges for downstream analysis
#'
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges strand makeGRangesFromDataFrame
#' @importFrom dplyr filter
#' @importFrom purrr map
#'
#' @param deep_obj deep_tss object
#' @param set Either 'building' for the train and test set, or 'all'
#' @param sequence_expansion Number of bases to expand on each side for surrounding sequence analysis
#' @param signal_expansion Number of bases to expand on each side for surrounding signal analysis
#'
#' @return Expanded GRanges objects
#'
#' @export

expand_ranges <- function(
	deep_obj,
	set = "building",
	sequence_expansion = 10,
	signal_expansion = 10
) {
	sliced <- deep_obj@experiment %>%
		filter(index %in% c("train", "test")) %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)

	sequence_expanded <- expand_range(sliced, sequence_expansion)
	signal_expanded <- expand_range(sliced, signal_expansion)

	deep_obj@ranges <- list(
		"sequence" = sequence_expanded,
		"signal" = signal_expanded
	)
	deep_obj@settings$sequence_expansion <- sequence_expansion
	deep_obj@settings$signal_expansion <- signal_expansion

	return(deep_obj)
}

#' Expand Ranges Function
#'
#' Backend function to expand granges
#'
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges strand resize sort width score
#'
#' @param grange GRange object to expand
#' @param expand_size Size to expand ranges on each side of TSS
#'
#' @return Expanded GRanges object

expand_range <- function(grange, expand_size) {
	tss_pos <- grange[strand(grange) == "+"] %>%
                resize(width = expand_size, fix = "start") %>%
                resize(width = width(.) + expand_size, fix = "end")

        tss_neg <- grange[strand(grange) == "-"] %>%
                resize(width = expand_size, fix = "end") %>%
                resize(width = width(.) + expand_size, fix = "start")

        tss_ranges <- c(tss_pos, tss_neg) %>% sort

        return(tss_ranges)
}
