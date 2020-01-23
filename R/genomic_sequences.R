#' Retrieve Sequences
#'
#' Retrieve sequences for expanded ranges
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools FaFile
#' @importFrom Biostrings getSeq
#' @importFrom purrr set_names map map2
#'
#' @param deep_obj dep_tss object
#' @param assembly Path to genome assembly
#'
#' @rdname get_sequences-function
#'
#' @export

get_sequences <- function(deep_obj, assembly) {
        genome <- FaFile(assembly)

        genomic_seqs <- deep_obj@ranges$sequence %>%
                map(
                        ~ getSeq(genome, .) %>%
                        as.character %>%
                        set_names(NULL)
                )

	deep_obj@ranges$sequence$train$seqs <- genomic_seqs$train
	deep_obj@ranges$sequence$test$seqs <- genomic_seqs$test
	deep_obj@ranges$sequence$all$seqs <- genomic_seqs$all

	shape_seqs <- deep_obj@ranges$shape %>%
		map(
			~ getSeq(genome, .) %>%
			as.character %>%
			set_names(NULL)
		)

	deep_obj@ranges$shape$train$seqs <- shape_seqs$train
	deep_obj@ranges$shape$test$seqs <- shape_seqs$test
	deep_obj@ranges$shape$all$seqs <- shape_seqs$all

        return(deep_obj)
}
