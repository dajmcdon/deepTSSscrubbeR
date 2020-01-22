setClass(
	"deep_tss",
	representation(
		experiment = "data.frame",
		settings = "list",
		ranges = "list",
		encoded = "list",
		model = "list",
		results = "list"
	),
	prototype(
		experiment = data.frame(),
		settings = list(),
		ranges = list(),
		encoded = list(),
		model = list(),
		results = list()
	)
)

#' deepTSSscrubbeR constructor function
#'
#' @importFrom Rsamtools scanBam scanBamWhat scanBamFlag ScanBamParam
#' @import tibble
#' @importFrom purrr pluck
#' @importFrom dplyr mutate left_join mutate_if select
#' @importFrom GenomicAlignments readGAlignmentPairs
#'
#' @param bam Bam file from five-prime mapping data
#'
#' @export

deep_tss <- function(bam) {
	
	## Get bam pair info.
	bampe <- readGAlignmentPairs(bam, use.names = TRUE) %>%
		as.data.frame %>%
		as_tibble(rownames = "qname")

	## Get bam first in read seq.
	params <- ScanBamParam(
		what = c("qname", "flag", "seq"),
		reverseComplement = TRUE,
		flag = scanBamFlag(isFirstMateRead = TRUE)
	)

	bamseq <- scanBam(bam, param = params)[[1]]
	bamseq$seq <- as.character(bamseq$seq)

	bamseq <- tibble(
		qname = pluck(bamseq, "qname"),
		flag_firstinread = pluck(bamseq, "flag"),
		seq_firstinread = pluck(bamseq, "seq")
	)
	
	## Combine paired bam and bam seq info.
	combined <- left_join(bampe, bamseq, by = "qname") %>%
		mutate(
			"start" = ifelse(flag_firstinread == "99", start.first, end.first),
			"end" = start,
			"strand" = strand.first,
			"seqnames" = seqnames.first,
			"tss" = start
		) %>%
		mutate_if(is.integer, as.double) %>%
		mutate_if(is.factor, as.character) %>%
		select(
			seqnames, start, end, strand, qname,
			cigar.first, flag_firstinread, seq_firstinread
		) %>%
		add_count(seqnames, strand, start, name = "score")

	deep_tss_object <- new("deep_tss", experiment = combined)
	return(deep_tss_object)
}
