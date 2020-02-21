
#' One Hot Encode a Sequence
#'
#' Internal function to one hot encode a sequence
#'
#' @importFrom stats model.matrix
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
#' @importFrom caret dummyVars contr.ltfr
#' @importFrom stats predict
#' @importFrom GenomicRanges GRanges
#' @importFrom purrr pmap map
#' @importFrom dplyr mutate_all
#' @importFrom stringr str_pad
#' @importFrom tidyr separate
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
	sequence_length <- (deep_obj@settings$sequence_expansion * 2) + 1

	encoded_sequences <- map(
		deep_obj@ranges$sequence,
		~ as_tibble(., .name_repair = "unique") %>%
			select(seqs) %>%
			separate(seqs, 1:sequence_length, into = sprintf("P%s", 1:sequence_length)) %>%
			mutate_all(~ factor(., levels = c("A", "T", "G", "C", "N")))
	)

	onehot <- map(
		encoded_sequences,
		~ dummyVars(' ~ .', data = .x) %>%
			predict(newdata = .x) %>%
			as_tibble(.name_repair = "unique") %>%
			as.matrix %>%
			array_reshape(dim = c(nrow(.), sequence_length, 5, 1))
	)

	deep_obj@encoded$genomic <- onehot
	return(deep_obj)
}

#' One Hot Encode Soft-clipped Bases
#'
#' One hot encode the soft-clipped sequences from the bam file
#'
#' @import tibble
#' @importFrom caret dummyVars contr.ltfr
#' @importFrom stats predict
#' @importFrom GenomicRanges GRanges
#' @importFrom purrr pmap map
#' @importFrom stringr str_pad
#' @importFrom reticulate array_reshape
#'
#' @param deep_obj object
#'
#' @rdname encode_soft-function
#'
#' @export

encode_soft <- function(deep_obj) {
	encoded_soft <- map(
		deep_obj@ranges$sequence,
		~ as_tibble(., .name_repair = "unique") %>%
			select(soft_bases) %>%
			mutate(
				soft_bases = replace_na(soft_bases, "U"),
				soft_bases = str_pad(soft_bases, 3, "right", "U")
			) %>%
			separate(soft_bases, 1:3, into = sprintf("P%s", 1:3)) %>%
			mutate_all(~ factor(., levels = c("A", "T", "G", "C", "N", "U")))
	)

	onehot <- map(
                encoded_soft,
                ~ dummyVars(' ~ .', data = .x) %>%
                        predict(newdata = .x) %>%
                        as_tibble(.name_repair = "unique") %>%
                        as.matrix %>%
                        array_reshape(dim = c(nrow(.), 3, 6, 1))
        )


	deep_obj@encoded$softclipped <- onehot
	return(deep_obj)

}

#' Encode TSS Status
#'
#' Encode the 'true' status of the TSSs
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom magrittr %>% extract
#' @importFrom reticulate array_reshape
#' @importFrom purrr map
#'
#' @param deep_obj deep_tss object
#'
#' @rdname encode_status-function
#'
#' @export

encode_status <- function(deep_obj) {
	status <- deep_obj@ranges$sequence %>%
		extract(c("train", "test")) %>%
		map(
			~ as_tibble(., .name_repair = "unique") %>%
			pull(status) %>%
			array_reshape(length(.))
		)	

	deep_obj@encoded$status <- status
	return(deep_obj)
}

#' Encode Surrounding Signal
#'
#' Encode signal around TSS
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr select matches mutate_all
#' @importFrom reticulate array_reshape
#'
#' @param deep_obj deep tss object
#'
#' @rdname encode_signal-function
#'
#' @export

encode_signal <- function(deep_obj) {
	signal <- map(
		deep_obj@ranges$signal,
		~ as_tibble(., .name_repair = "unique") %>%
		select(matches("X\\d+$")) %>%
		mutate_all(~ ifelse(. > 0, 1, 0)) %>%
		as.matrix %>%
		array_reshape(dim = c(nrow(.), 1, ncol(.), 1))
	)

	deep_obj@encoded$signal <- signal
	return(deep_obj)
}
