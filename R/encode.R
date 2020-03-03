
#' One Hot Encode a Sequence (Deprecated)
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
#' @export

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
#' @import caret
#' @importFrom stats predict
#' @importFrom GenomicRanges GRanges
#' @importFrom stringr str_pad str_split
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

	encoded_sequences <- as.data.table(deep_obj@ranges$sequence)
	encoded_sequences <- encoded_sequences[, str_split(genomic_seq, "", simplify = TRUE)]
	encoded_sequences <- as.data.table(encoded_sequences)
	encoded_sequences <- encoded_sequences[, lapply(.SD, factor, levels = c("A", "T", "G", "C", "N"))]

	onehot <- encoded_sequences %>%
		dummyVars(' ~ .', data = .) %>%
		predict(newdata = encoded_sequences) %>%
		as.matrix %>%
		array_reshape(dim = c(nrow(.), sequence_length, 5, 1))

	deep_obj@encoded$genomic <- onehot
	return(deep_obj)
}

#' One Hot Encode Soft-clipped Bases
#'
#' One hot encode the soft-clipped sequences from the bam file
#'
#' @import tibble
#' @import caret
#' @importFrom stats predict
#' @importFrom stringr str_pad str_split
#' @importFrom reticulate array_reshape
#' @importFrom tidyr replace_na
#'
#' @param deep_obj object
#'
#' @rdname encode_soft-function
#'
#' @export

encode_soft <- function(deep_obj) {
	
	encoded_soft <- as.data.table(deep_obj@ranges$soft)

	encoded_soft[,
		soft_bases := replace_na(soft_bases, "U")
	][,
		soft_bases := str_pad(soft_bases, 3, "right", "U")
	]

	encoded_soft <- encoded_soft[, str_split(soft_bases, "", simplify = TRUE)]
	encoded_soft <- as.data.table(encoded_soft)
	encoded_soft <- encoded_soft[, lapply(.SD, factor, levels = c("A", "T", "G", "C", "N", "U"))]

	onehot <- encoded_soft %>%
		dummyVars(' ~ .', data = .) %>%
		predict(newdata = encoded_soft) %>%
		as.matrix %>%
		array_reshape(dim = c(nrow(.), 3, 6, 1))

	deep_obj@encoded$soft <- onehot
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

	## Prepare surrounding signal matrix.
	signal <- as.data.table(deep_obj@ranges$signal)
	signal[, c(
			"seqnames", "start", "end", "strand",
			"tss", "signal_index", "width"
		) := NULL
	]

	## Convert to binary representation of signal.
	signal <- as.matrix(signal)
	signal <- (signal >= 1) * 1
	
	## Convert to array and save to deep tss object.
	signal <- array_reshape(signal, dim = c(nrow(signal), 1, ncol(signal), 1))

	deep_obj@encoded$signal <- signal
	return(deep_obj)
}
