setClass(
	"deep_tss_object",
	representation(
		experiment = "list",
		indices = "list",
	),
	prototype(
		experiment = list(),
		indices = list()
	)
)

#' deepTSSscrubbeR constructor function
#'
#' @param TSSs Named list of TSS granges
#'
#' @importFrom GenomicRanges GRanges
#'
#' @export

deep_tss <- function(TSSs) {
	deep_tss_object <- new(
		"deep_tss_object",
		TSSs = TSSs
	)

	return(deep_tss_obejct)
}
