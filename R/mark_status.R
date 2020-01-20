
#' Mark Status
#'
#' Mark status as spurious based on read threshold
#'
#' @importFrom dplyr mutate case_when
#'
#' @param deep_obj deep_tss object
#' @param lower Reads less than or equal to 'lower' will be marked as spurious
#' @param upper Reads greater than or equal to 'upper' will be marked as good
#'
#' @return GRanges object with status column indicating spurious (1) or not (0)
#'
#' @export

mark_status <- function(deep_obj, lower = 2, upper = 5) {
	status_marked <- deep_obj@experiment %>%
		mutate("status" = case_when(
			score <= lower ~ 0,
			score >= upper ~ 1
		))

	deep_obj@experiment <- status_marked
	return(deep_obj)
}
