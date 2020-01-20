
#' Get Surrounding Signal
#'
#' Get TSS signal around selected TSS
#'
#' @import tibble
#' @importFrom dplyr select contains count group_by ungroup rename left_join
#' @importFrom tidyr complete pivot_wider
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom IRanges findOverlapPairs
#'
#' @rdname get_signal-function
#'
#' @export

get_signal <- function(deep_obj) {
	overlaps <- deep_obj@ranges$signal %>%
		findOverlapPairs(deep_obj@ranges$all) %>%
		as.data.frame %>%
		as_tibble(.name_repair = "unique") %>%
		select(
			first.X.qname,
			first.X.seqnames,
			first.X.start, second.X.start,
			first.X.end,
			first.X.strand
		) %>%
		mutate("position" = ifelse(
			first.X.strand == "+",
			second.X.start - first.X.start,
			first.X.end - second.X.start
		)) %>%
		select(-second.X.start)

	signal_length <- deep_obj@settings$signal_expansion * 2

	positions <- overlaps %>%
		count(
			first.X.qname, first.X.seqnames, first.X.start,
			first.X.end, first.X.strand, position
		) %>%
		group_by(
			first.X.qname, first.X.seqnames, first.X.start,
			first.X.end, first.X.strand
		) %>%
		complete(position = seq(0, signal_length - 1, 1), fill = list(n = 0)) %>%
		ungroup %>%
		select(first.X.qname, position, n) %>%
		rename("qname" = first.X.qname) %>%
		pivot_wider(names_from = position, values_from = n)

	overlap <- deep_obj@ranges$signal %>%
		as_tibble(.name_repair = "unique") %>%
		left_join(positions, by = "qname") %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)

	deep_obj@ranges$signal <- overlap
	return(deep_obj)
}
