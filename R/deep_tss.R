setClass(
	"deep_tss_object",
	representation(
		experiment = "list",
		status = "list",
		ranges = "list",
		split = "list"
	),
	prototype(
		experiment = list(),
		status = list(),
		ranges = list(),
		split = list()
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
	deep_tss_object <- new("deep_tss_object", experiment = TSSs)
	return(deep_tss_object)
}
