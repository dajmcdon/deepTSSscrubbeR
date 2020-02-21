
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
	overlaps <- map(deep_obj@ranges$signal,
		~ findOverlapPairs(., deep_obj@ranges$original$all) %>%
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
	)

	signal_length <- (deep_obj@settings$signal_expansion * 2) + 1
	dummy <- tibble(first.X.qname = "__dummy__", position = seq(0, signal_length - 1, 1))

	positions <- map(overlaps,
		~ count(.,
			first.X.qname, first.X.seqnames, first.X.start,
			first.X.end, first.X.strand, position
		) %>%
		bind_rows(dummy, .) %>%
		pivot_wider(names_from = position, values_from = n) %>%
		rename("qname" = first.X.qname) %>%
		select(qname, matches("^\\d+$")) %>%
		mutate_at(vars(-qname), ~ replace_na(., 0)) %>%
		filter(qname != "__dummy__")
	)

	overlap <- map2(deep_obj@ranges$signal, positions,
		~ as_tibble(.x, .name_repair = "unique") %>%
			left_join(.y, by = "qname") %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE)
	)

	deep_obj@ranges$signal <- overlap
	return(deep_obj)
}
