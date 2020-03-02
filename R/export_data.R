
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

#' Export Raw
#'
#' Export the raw values that eventially get encoded.
#'
#' @importFrom stringr str_c
#'
#' @param deep_obj tss_obj
#'
#' @rdname export_raw-function
#' @export

export_raw <- function(deep_obj, sequence_file, signal_file) {
	
	## Pull out the sequences surrounding the TSS and soft-clipped bases.
	tss_sequences <- deep_obj@ranges$sequence$all %>%
		as_tibble(.name_repair = "unique") %>%
		select(qname, seqnames, start, end, strand, tss, score, soft_bases, seqs) %>%
		rename(tss_position_score = score, surrounding_seqs = seqs) %>%
		mutate(start = tss, end = tss) %>%
		select(-tss)

	## Pull out the sequences surrounding the TSSs that will be used as input for DNAShapeR.
	tss_shape_sequences <- deep_obj@ranges$shape$all %>%
		as_tibble(.name_repair = "unique") %>%
		select(qname, seqs) %>%
		rename(shape_seqs = seqs)

	## merge sequence based data.
	merged_sequences <- left_join(tss_sequences, tss_shape_sequences, by = "qname")

	## Mark status.
	merged_sequences <- merged_sequences %>%
		mutate(status = case_when(
			tss_position_score >= deep_obj@settings$upper ~ 1,
			tss_position_score <= deep_obj@settings$lower ~ 0,
			TRUE ~ 2
		))

	## Pull out the signal surrounding the TSS.
	signal_indices <- deep_obj@settings$signal_expansion
	signal_indices <- str_c("X", as.character(seq(0, signal_indices * 2, 1)))

	surrounding_signal <- deep_obj@ranges$signal$all %>%
		as_tibble(.name_repair = "unique") %>%
		select(qname, all_of(signal_indices))

	## Export sequence data.
	write.table(
		merged_sequences, sequence_file, sep = "\t",
		col.names = TRUE, row.names = FALSE, quote = FALSE
	)

	# Export signal data.
	write.table(
		surrounding_signal, signal_file, sep = "\t",
		col.names = TRUE, row.names = FALSE, quote = FALSE
	)
}

#' Export Encoded
#'
#' Export encoded data
#'
#' @importFrom reticulate import
#'
#' @param deep_obj tss_obj
#'
#' @rdname export_encoded-function
#' @export

export_encoded <- function(deep_obj) {
	
	## Import numpy.
	np <- import("numpy")

	## Export encoded genomic data.
	np$save("genomic.npy", deep_obj@encoded$genomic$all)

	## Export shape data.
	np$save("shape.npy", deep_obj@encoded$shape$all)

	## Export softclipped data.
	np$save("softclipped.npy", deep_obj@encoded$softclipped$all)

	## Export signal data.
	np$save("signal.npy", deep_obj@encoded$signal$all)
}
