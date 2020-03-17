
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
	} else {
		tss_status <- tss_all[, .(seqnames, start, end, strand, status_prob)]
		tss_status <- tss_status[,
			.(score = mean(status_prob)),
			by = .(seqnames, start, end, strand)
		]

		tss_status <- makeGRangesFromDataFrame(tss_status, keep.extra.columns = TRUE)

		tss_status_pos <- tss_status[strand(tss_status) == "+"]
		tss_status_neg <- tss_status[strand(tss_status) == "-"]

		export(tss_status_pos, "predicted_status_pos.bedgraph")
		export(tss_status_neg, "predicted_status_neg.bedgraph")
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
#' @param raw_report Raw report directory and file name
#'
#' @rdname export_raw-function
#' @export

export_raw <- function(deep_obj, raw_report = "raw_report.tsv") {
	
	# Prepare data for export.
	raw_df <- as.data.table(deep_obj@experiment)[,
		.(rowid, qname, seqnames, start, end, strand, tss, score,
		tss_group, soft_group, gene_id, transcript_id,
		tss_distance, simple_annotations) 
	]

	# Export signal data.
	write.table(
		raw_df, raw_report, sep = "\t", col.names = TRUE,
		row.names = FALSE, quote = FALSE
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
	np$save("genomic.npy", deep_obj@encoded$genomic)

	## Export shape data.
	np$save("shape.npy", deep_obj@encoded$shape)

	## Export softclipped data.
	np$save("softclipped.npy", deep_obj@encoded$soft)

	## Export signal data.
	np$save("signal.npy", deep_obj@encoded$signal)
}

#' Export Report
#'
#' Export a table report for predictions.
#'
#' @param deep_obj deep tss object
#'
#' @rdname export_report-function
#' @export

export_report <- function(deep_obj) {
	
	## Pull out useful info.
	results <- as.data.table(deep_obj@results$all)
	results[, c("cigar.first", "flag_firstinread", "tss", "status", "index") := NULL]

	## Merge soft-clipped base info.
	soft <- as.data.table(deep_obj@ranges$soft)[, .(soft_group, soft_bases)]
	setkey(soft, soft_group)
	setkey(results, soft_group)
	results <- merge(results, soft, all.x = TRUE)

	## Merge sequence info.
	sequences <- as.data.table(deep_obj@ranges$sequence)[, .(tss_group, genomic_seq)]
	setkey(results, tss_group)
	results <- merge(results, sequences, all.x = TRUE)

	## Clean up results.
	results <- results[,
		.(seqnames, start, end, strand, width, score, gene_id,
		transcript_id, tss_distance, annotation, simple_annotations,
		rowid, soft_group, soft_bases, tss_group, genomic_seq, status_prob)
	]
	results <- distinct(results, soft_group, .keep_all = TRUE)

	## Export the report table.
	write.table(
		results, "prediction_report.tsv", sep = "\t",
		col.names = TRUE, row.names = FALSE, quote = FALSE
	)		
}
