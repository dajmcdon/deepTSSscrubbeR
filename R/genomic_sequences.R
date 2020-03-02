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
        genomic_seqs <- deep_obj@ranges$sequence %>%
		getSeq(genome, .) %>%
		as.character %>%
		set_names(NULL)
	deep_obj@ranges$sequence$genomic_seq <- genomic_seqs

	## Retrieve the sequences used for shape encoding.
	shape_seqs <- deep_obj@ranges$shape %>%
		getSeq(genome, .) %>%
		as.character %>%
		set_names(NULL)
	deep_obj@ranges$shape$shape_seq <- shape_seqs

	## Return the deep tss object.
        return(deep_obj)
}
