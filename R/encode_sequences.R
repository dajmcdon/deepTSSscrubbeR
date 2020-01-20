#' Retrieve Sequences
#'
#' Retrieve sequences for expanded ranges
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools FaFile
#' @importFrom Biostrings getSeq
#'
#' @param deep_obj dep_tss object
#' @param assembly Path to genome assembly
#'
#' @rdname get_sequences-function
#'
#' @export

get_sequences <- function(deep_obj, assembly) {
        genome <- FaFile(assembly)

        seqs <- deep_obj@ranges$sequence %>%
                getSeq(genome, .) %>%
                as.character %>%
                set_names(NULL)

        deep_obj@ranges$sequence$seq <- seqs
        return(deep_obj)
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

one_hot_seqs <- function(sequence) {
	one_hot_seq <- sequence %>%
		str_split(pattern = "", simplify = TRUE) %>%
		factor(levels = c("A", "T", "G", "C", "N")) %>%
		data.frame("seq" = .) %>%
		model.matrix(~ 0 + seq, data = .) %>%
		t %>%
		array_reshape(dim = c(5, ncol(.)))

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

encode_sequences <- function(deep_obj) {
	encoded_sequences <- deep_obj@sequences %>%
		map()
			
}
