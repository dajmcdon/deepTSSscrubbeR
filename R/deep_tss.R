setClass(
	"deep_tss_object",
	representation(
		experiment = "list",
		indices = "list",
		ranges = "list"
	),
	prototype(
		experiment = list(),
		indices = list(),
		ranges = list()
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
