#' Retrieve Sequences
#'
#' Retrieve sequences for expanded ranges
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools FaFile
#' @importFrom Biostrings getSeq
#' @importFrom purrr set_names
#'
#' @param deep_obj dep_tss object
#' @param assembly Path to genome assembly
#'
#' @rdname get_sequences-function
#'
#' @export

get_sequences <- function(deep_obj, assembly) {
        
	## Load the genome assembly.
	genome <- FaFile(assembly)

	## Retrieve the sequences used for genomic encoding.
	genomic_ranges <- deep_obj@ranges$sequence
        genomic_seqs <- genomic_ranges %>%
		getSeq(genome, .) %>%
		as.character %>%
		set_names(NULL)

	genomic_ranges$genomic_seq <- genomic_seqs

	## Merge genomic sequence index back into original data.
	original_data <- as.data.table(deep_obj@experiment)
	genomic <- as.data.table(genomic_ranges)[,
		.(seqnames, tss, strand, sequence_index)
	]

	original_data <- merge(original_data, genomic, on = c("seqnames", "tss", "strand"))

	## Retrieve the sequences used for shape encoding.
	shape_ranges <- deep_obj@ranges$shape
	shape_seqs <- genomic_ranges %>%
		getSeq(genome, .) %>%
		as.character %>%
		set_names(NULL)

	shape_ranges$shape_seq <- shape_seqs

	## Merge shape index back into original data.
	shape <- as.data.table(shape_ranges)[,
		.(seqnames, tss, strand, shape_index)
	]

	original_data <- merge(original_data, shape, on = c("seqnames", "tss", "strand"))

	## Return the deep tss object.
	deep_obj@experiment <- original_data
	deep_obj@ranges$sequence <- genomic_ranges
	deep_obj@ranges$shape <- shape_ranges

        return(deep_obj)
}
