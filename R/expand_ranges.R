
#' Expand Ranges
#'
#' Expand TSS GRanges for downstream analysis
#'
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr filter
#' @importFrom purrr map
#'
#' @param deep_obj deep_tss object
#' @param sequence_expansion Number of bases to expand on each side for surrounding sequence analysis
#' @param signal_expansion Number of bases to expand on each side for surrounding signal analysis
#'
#' @return Expanded GRanges objects
#'
#' @export

expand_ranges <- function(
	deep_obj,
	sequence_expansion = 10,
	signal_expansion = 10
) {

	sequence_expanded <- deep_obj@ranges$subset %>%
		expand_range(., sequence_expansion)

	signal_expanded <- deep_obj@ranges$subset %>%
		expand_range(., signal_expansion)

	deep_obj@ranges$sequence <- sequence_expanded
	deep_obj@ranges$signal <- signal_expanded

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
#' @importFrom purrr map
#'
#' @param grange GRange object to expand
#' @param expand_size Size to expand ranges on each side of TSS
#'
#' @return Expanded GRanges object

expand_range <- function(grange, expand_size) {

	new_ranges <- map(grange, function(x) {
		tss_ranges <- x %>%
                	resize(width = expand_size + 1, fix = "end") %>%
                	resize(width = (expand_size * 2) + 1, fix = "start")

		return(tss_ranges)
	})

        return(new_ranges)
}
