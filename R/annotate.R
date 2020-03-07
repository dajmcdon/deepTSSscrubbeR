
#' Annotate the TSSs
#'
#' Annotate TSSs to the nearest transcript.
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom stringr str_detect
#' @importFrom dplyr case_when
#' @importFrom GenomicFeatures makeTxDbFromGFF
#'
#' @param deep_obj deep tss object
#' @param annotation genome annotation
#' @param upstream Bases upstream to consider promoter
#' @param downstream Bases downstream to consider promoter
#'
#' @rdname tss_annotate-function
#' @export

tss_annotate <- function(deep_obj, annotation, upstream = 1000, downstream = 1000) {

	## Turn annotation data into TxDb.
	genome_annotation <- makeTxDbFromGFF(annotation)

	## Prepare data for annotation.
	original_data <- as.data.table(deep_obj@experiment)

	ranges <- original_data[,
		.(seqnames, start, end, strand)
	]
	ranges <- unique(ranges)
	ranges <- makeGRangesFromDataFrame(ranges)

	## Annotate ranges.
	annotated_ranges <- annotatePeak(
		ranges, tssRegion = c(-upstream, downstream),
		TxDb = genome_annotation, sameStrand = TRUE,
		level = "transcript"
	)

	annotated_ranges <- as.data.table(annotated_ranges)[,
		.(seqnames, start, end, strand, geneId, transcriptId, annotation, distanceToTSS)
	]

	setnames(
		annotated_ranges,
		old = c("geneId", "transcriptId", "distanceToTSS"),
		new = c("gene_id", "transcript_id", "tss_distance")
	)

	annotated_ranges[,
		simple_annotations := case_when(
			annotation == "Promoter" ~ "promoter",
			str_detect(annotation, pattern="(Exon|UTR)") ~ "exon",
			str_detect(annotation, pattern="Intron") ~ "intron",
			str_detect(annotation, pattern="Downstream") ~ "downstream",
			annotation == "Distal Intergenic" ~ "intergenic"
		)
	]

	## Merge annotation data back into original data.
	keys <- c("seqnames", "start", "end", "strand")
	setkeyv(annotated_ranges, keys)
	setkeyv(original_data, keys)

	merged <- merge(original_data, annotated_ranges, all.x = TRUE)

	## Put data back into deep tss object.
	deep_obj@experiment <- as.data.frame(merged)
	return(deep_obj)

}
