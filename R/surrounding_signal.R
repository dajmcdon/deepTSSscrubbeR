
#' Get Surrounding Signal
#'
#' Get TSS signal around selected TSS
#'
#' @import tibble
#' @importFrom dplyr select contains count rename left_join bind_rows vars mutate_at
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom IRanges findOverlapPairs
#' @importFrom purrr map map2
#'
#' @param deep_obj deep tss object
#'
#' @rdname get_signal-function
#'
#' @export

get_signal <- function(deep_obj) {

	## Grab all TSS positions and scores.
	all_ranges <- as.data.table(deep_obj@experiment)
	all_ranges <- all_ranges[, .(seqnames, start, end, strand, score)]
	all_ranges <- unique(all_ranges)
	all_ranges <- makeGRangesFromDataFrame(all_ranges, keep.extra.columns = TRUE)

	## Find surrounding TSSs.
	overlaps <- deep_obj@ranges$signal %>%
		findOverlapPairs(., all_ranges) %>%
		as.data.table

	overlaps <- overlaps[,
		.(first.X.seqnames, first.X.start, second.X.start, first.X.tss,
		second.X.score, first.X.end, first.X.strand, first.signal_index)
	]

	overlaps[,
		position := ifelse(
			first.X.strand == "+",
			second.X.start - first.X.start,
			first.X.end - second.X.start
		)
	][,
		second.X.start := NULL
	]

	setnames(
		overlaps,
		old = c(
			"first.X.seqnames", "first.X.start", "first.X.end",
			"first.X.strand", "second.X.score", "first.X.tss",
			"first.signal_index"
		),
		new = c("seqnames", "start", "end", "strand", "score", "tss", "signal_index")
	)

	## Create matrix of surrounding signal.
	signal_length <- (deep_obj@settings$signal_expansion * 2) + 1
	dummy <- data.table(
		seqnames = "__dummy__",
		position = seq(0, signal_length - 1, 1)
	)
	positions <- bind_rows(dummy, overlaps)

	positions <- dcast(
		overlaps, seqnames + start + end + strand + tss + signal_index ~ position,
		fill = 0, value.var = "score"
	)[
		seqnames != "__dummy__"
	]

	## Merge into original data.
	original_data <- as.data.table(deep_obj@experiment)
	pos <- positions[, .(seqnames, tss, strand, signal_index)]
	original_data <- merge(original_data, pos, on = c("seqnames", "tss", "strand"))

	## Return surrounding signal to deep tss object.
	positions <- makeGRangesFromDataFrame(positions, keep.extra.columns = TRUE)
	deep_obj@ranges$signal <- positions
	deep_obj@experiment <- as.data.frame(original_data)

	return(deep_obj)
}
