
#' Mark Status
#'
#' Mark status as spurious based on read threshold
#'
#' @importFrom dplyr case_when
#'
#' @param deep_obj deep_tss object
#' @param lower Reads less than or equal to 'lower' will be marked as spurious
#' @param upper Reads greater than or equal to 'upper' will be marked as good
#' @param use_annotation If TRUE will ignore upper and lower bounds and instead
#' pick a sample of TSSs within and outside of annotated promoters
#'
#' @return GRanges object with status column indicating spurious (1) or not (0)
#'
#' @export

mark_status <- function(deep_obj, lower = 2, upper = 10, use_annotation = FALSE) {

	if (!use_annotation) {
		status_marked <- as.data.table(deep_obj@experiment)
		status_marked[, status := case_when(
			score <= lower ~ 0,
			score >= upper ~ 1
		)]

		deep_obj@settings$lower <- lower
		deep_obj@settings$upper <- upper
		deep_obj@settings$use_annotation <- use_annotation
	} else {
		status_marked <- as.data.table(deep_obj@experiment)
		status_marked[, status := case_when(
			simple_annotations == "exon" ~ 0,
			simple_annotations == "promoter" & score >= upper ~ 1
		)]
		deep_obj@settings$use_annotation <- use_annotation
	}

	deep_obj@experiment <- as.data.frame(status_marked)
	return(deep_obj)
}
