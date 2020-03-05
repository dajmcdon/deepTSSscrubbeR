
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
#'
#' @rdname export_bedgrpahs-function
#'
#' @export

export_bedgraphs <- function(deep_obj) {

	tss_all <- as.data.table(deep_obj@results$all)
	tss_all[, c("start", "end") := tss]

	if (deep_obj@settings$model_type == "score") {
		tss_score <- tss_all[, .(seqnames, start, end, strand, predicted_log2_score)]
		tss_score <- tss_score[,
			.(score = 2^mean(predicted_log2_score)),
			by = .(seqnames, start, end, strand)
		]

		tss_score <- makeGRangesFromDataFrame(tss_score, keep.extra.columns = TRUE)
		
		tss_score_pos <- tss_score[strand(tss_score) == "+"]
		tss_score_neg <- tss_score[strand(tss_score) == "-"]

		export(tss_score_pos, "predicted_score_pos.bedgraph")
		export(tss_score_neg, "predicted_score_neg.bedgraph")
	}


	tss <- tss_all[, .(seqnames, start, end, strand, score)]
	tss <- unique(tss)	
	tss <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)
	
	pos <- tss[strand(tss) == "+"]
	neg <- tss[strand(tss) == "-"]

	export(pos, "pos.bedgraph")
	export(neg, "neg.bedgraph")

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
