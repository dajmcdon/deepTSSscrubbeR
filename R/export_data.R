
#' Export TSS Bedgraphs
#'
#' Export bedgraphs
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame strand
#' @importFrom dplyr mutate distinct group_by ungroup summarize filter add_count
#' @importFrom rtracklayer export
#'
#' @param deep_obj tss_obj
#' @param n_reads filter TSSs with less than this number of reads
#' @param cutoff filter reads with less than this probability of being a real TSS
#'
#' @rdname export_bedgrpahs-function
#'
#' @export

export_bedgraphs <- function(deep_obj, n_reads = 1, cutoff = 0.9) {
	tss <- deep_obj@results$all %>%
		as_tibble(.name_repair = "unique") %>%
		mutate(start = tss, end = tss)

	unfiltered <- tss %>%
		distinct(seqnames, start, end, strand, score) %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)
	
	filtered <- tss %>%
		filter(score >= n_reads & probs >= cutoff) %>%
		add_count(seqnames, start, end, strand, name = "score") %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)

	probs <- tss %>%
		group_by(seqnames, start, end, strand) %>%
		summarize(score = mean(probs)) %>%
		ungroup %>%
		makeGRangesFromDataFrame(keep.extra.columns = TRUE)

	unfiltered_pos <- unfiltered[strand(unfiltered) == "+"]
	unfiltered_neg <- unfiltered[strand(unfiltered) == "-"]

	export(unfiltered_pos, "unfiltered_pos.bedgraph")
	export(unfiltered_neg, "unfiltered_neg.bedgraph")

	filtered_pos <- filtered[strand(filtered) == "+"]
	filtered_neg <- filtered[strand(filtered) == "-"]

	export(filtered_pos, "filtered_pos.bedgraph")
	export(filtered_neg, "filtered_neg.bedgraph")

	probs_pos <- probs[strand(probs) == "+"]
	probs_neg <- probs[strand(probs) == "-"]

	export(probs_pos, "probs_pos.bedgraph")
	export(probs_neg, "probs_neg.bedgraph")
}
