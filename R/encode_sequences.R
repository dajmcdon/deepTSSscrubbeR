
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

	seq_data <- experiment@split$sequence %>%
		map(~ getSeq(genome_assembly, .))

	experiment@sequences <- seq_data
	experiment@settings$assembly <- assembly

	return(experiment)
}

#' One Hot Encode a Sequence
#'
#' Internal function to one hot encode a sequence
#'
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom reticulate array_reshape
#'
#' @param sequence String representation of bases
#'
#' @return one hot encoded sequence

one_hot_seqs <- function(sequence, length) {
	split_seq <- sequence %>%
		str_split(pattern = "") %>%
		unlist

	one_hot_seq <- c("A", "T", "G", "C") %>%
		map(~ {split_seq == .x} %>% as.numeric) %>%
		do.call(rbind, .) %>%
		array_reshape(dim = c(4, length))

	return(one_hot_seq)
}

#' One Hot Encode Sequences
#'
#' One hot encode the sequences from the expanded ranges
#'
#' @importFrom purrr map
#' @importFrom reticulate array_reshape
#'
#' @param experiment deep_tss object
#'
#' @return one hot encoded tensor of sequence stored in @encoded_sequences
#'
#' @export

encode_sequences <- function(experiment) {
	encoded_sequences <- experiment@sequences %>%
		map(
			~ as.character(.x) %>% as.character %>%
				map(~ one_hot_seqs(.))
		)
}
