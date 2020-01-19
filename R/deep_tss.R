setClass(
	"deep_tss",
	representation(
		experiment = "list",
		settings = "list",
		status = "list",
		ranges = "list",
		split_index = "list",
		split = "list",
		status_encoded = "list",
		sequences = "list",
		encoded_sequences = "list"
	),
	prototype(
		experiment = list(),
		settings = list(),
		status = list(),
		ranges = list(),
		split_index = list(),
		split = list(),
		status_encoded = list(),
		sequences = list(),
		encoded_sequences = list()
	)
)

#' deepTSSscrubbeR constructor function
#'
#' @importFrom Rsamtools scanBam scanBamWhat scanBamFlag ScanBamParam
#'
#' @param bam Bam file from five-prime mapping data
#'
#' @export

deep_tss <- function(bam) {
	
	params = ScanBamParam(
		what = c("qname", "flag", "rname", "strand", "pos", "cigar", "seq"),
		reverseComplement = TRUE,
		flag = scanBamFlag(isFirstMateRead = TRUE)
	)

	bam <- scanBam(bam, param = params)[[1]]
	bam$seq <- as.character(bam$seq)

	bam <- bam %>%
		stack %>% unstack %>%
		as_tibble(.name_repair = "unique")

	deep_tss_object <- new("deep_tss", experiment = bam)
	return(deep_tss_object)
}
