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
#' @import methods
#' @importFrom Rsamtools scanBam scanBamWhat scanBamFlag ScanBamParam
#' @import data.table
#' @importFrom purrr pluck
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges
#'
#' @param bam Bam file from five-prime mapping data
#'
#' @export

deep_tss <- function(bam) {
	
	## Get bam pair info.
	bampe <- readGAlignmentPairs(bam, use.names = TRUE) %>%
		as.data.frame %>%
		as.data.table(keep.rownames = "qname")
	setkey(bampe, qname)

	## Get bam first in read seq.
	params <- ScanBamParam(
		what = c("qname", "flag", "seq"),
		reverseComplement = TRUE,
		flag = scanBamFlag(isFirstMateRead = TRUE)
	)

	bamseq <- scanBam(bam, param = params)[[1]]
	bamseq$seq <- as.character(bamseq$seq)

	bamseq <- data.table(
		qname = pluck(bamseq, "qname"),
		flag_firstinread = pluck(bamseq, "flag"),
		seq_firstinread = pluck(bamseq, "seq"),
		key = "qname"
	)
	
	## Combine paired bam and bam seq info.
	combined <- merge(bampe, bamseq)

	combined <- combined[,
		.(start = ifelse(flag_firstinread == "99", start.first, end.first),
		strand = strand.first, seqnames = seqnames.first,
		cigar.first, flag_firstinread, seq_firstinread, qname)
	]

	combined[, c("end", "tss") := start]

	combined <- makeGRangesFromDataFrame(combined, keep.extra.columns = TRUE)
	combined <- sort(combined)
	combined <- as.data.table(combined)

	combined[,
		c("score", "tss_group") := list(.N, .GRP),
		by = .(seqnames, start, end, strand)
	][,
		rowid := seq_len(.N)
	]

	deep_tss_object <- new("deep_tss", experiment = as.data.frame(combined))
	return(deep_tss_object)
}
