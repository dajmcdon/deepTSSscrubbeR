
#' Retrieve Sequences
#'
#' Retrieve sequences for expanded GRanges
#'
#' @importFrom purrr map
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools FaFile getSeq
#'
#' @param experiment deep_tss object
#' @param assembly File path to genome assembly
#'
#' @return Sequences stored in @sequences
#'
#' @export

retrieve_sequences <- function(experiment, assembly) {
	genome_assembly <- FaFile(assembly)

	seq_data <- experiment@split[[experiment@settings$reference_sample]] %>%
		map(~ getSeq(genome_assembly, .))

	experiment@sequences <- seq_data
	experiment@settings$assembly <- assembly

	return(experiment)
}
