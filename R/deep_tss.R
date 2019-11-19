setClass(
	"deep_tss_object",
	representation(
		experiment = "list",
		settings = "list",
		status = "list",
		ranges = "list",
		split = "list",
		status_encoded = "list"
	),
	prototype(
		experiment = list(),
		settings = list(),
		status = list(),
		ranges = list(),
		split = list(),
		status_encoded = list()
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
