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
#' @rdname one_hot_encode-function
#'
#' @return one hot encoded sequence

one_hot_seqs <- function(sequence) {
	one_hot_seq <- sequence %>%
		str_split(pattern = "", simplify = TRUE) %>%
		factor(levels = c("A", "T", "G", "C", "N", "U")) %>%
		data.frame("seq" = .) %>%
		model.matrix(~ 0 + seq, data = .) %>%
		t

	return(one_hot_seq)
}

#' One Hot Encode Genomic Sequences
#'
#' One hot encode the genomic sequences from the expanded ranges
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom purrr pmap
#' @importFrom reticulate array_reshape
#'
#' @param deep_obj deep_tss object
#'
#' @return one hot encoded tensor of sequence stored in @encoded_sequences
#'
#' @rdname encode_genomic-function
#'
#' @export

encode_genomic <- function(deep_obj) {
	sequence_length <- deep_obj@settings$sequence_expansion * 2

	encoded_sequences <- deep_obj@ranges$sequence %>%
		as_tibble(.name_repair = "unique") %>%
		pmap(function(...) {
			args <- list(...)
			onehot <- one_hot_seqs(args$seq)
			return(onehot)
		}) %>%
		array_reshape(dim = c(length(.), 6, sequence_length, 1))

	deep_obj@encoded$genomic <- encoded_sequences
	return(deep_obj)
}

#' One Hoe Encode Soft-clipped Bases
#'
#' One hot encode the soft-clipped sequences from the bam file
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom purrr pmap
#' @importFrom stringr str_pad
#' @importFrom reticulate array_reshape
#'
#' @param deep_obj object
#'
#' @rdname encode_soft-function
#'
#' @export

encode_soft <- function(deep_obj) {
	encoded_soft <- deep_obj@ranges$sequence %>%
		as_tibble(.name_repair = "unique") %>%
		pmap(function(...) {
			args <- list(...)
			if (is.na(args$soft_bases)) {
				soft <- "UUU"
			} else {
				soft <- str_pad(args$soft_bases, 3, "right", "U")
			}
			onehot <- one_hot_seqs(soft)
			return(onehot)
		}) %>%
		array_reshape(dim = c(length(.), 6, 3, 1))

	deep_obj@encoded$softclipped <- encoded_soft
	return(deep_obj)
}

#' Encode TSS Status
#'
#' Encode the 'true' status of the TSSs
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom reticulate array_reshape
#'
#' @param deep_obj deep_tss object
#'
#' @rdname encode_status-function
#'
#' @export

encode_status <- function(deep_obj) {
	status <- deep_obj@ranges$sequence %>%
		as_tibble(.name_repair = "unique") %>%
		pull(status) %>%
		array_reshape(length(.))

	deep_obj@encoded$status <- status
	return(deep_obj)
}

#' Encode Surrounding Signal
#'
#' Encode signal around TSS
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr select matches
#' @importFrom reticulate array_reshape
#'
#' @rdname encode_signal-function
#'
#' @export

encode_signal <- function(deep_obj) {
	signal <- deep_obj@ranges$signal %>%
		as_tibble(.name_repair = "unique") %>%
		select(matches("X\\d+$")) %>%
		mutate_all(~ ifelse(. > 0, 1, 0)) %>%
		as.matrix %>%
		array_reshape(dim = c(nrow(.), 1, ncol(.)))

	deep_obj@encoded$signal <- signal
	return(deep_obj)
}
