
#' Mark Status
#'
#' Mark status as spurious based on read threshold
#'
#' @importFrom GenomicRanges GRanges score
#' @importFrom purrr map
#'
#' @param experiment deep_tss object
#' @param threshold Reads less than or equal to threshold will be marked as spurious
#'
#' @return GRanges object with status column indicating spurious (1) or not (0)
#'
#' @export

mark_status <- function(experiment, threshold = 1) {
	granges <- experiment@experiment %>%
		map(function(x) {
			x$status <- ifelse(score(x) <= threshold, 1, 0)
			return(x)
		})
	
	experiment@status <- granges
	return(experiment)
}
