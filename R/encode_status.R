
#' Encode TSS Status
#'
#' Encode TSS status into proper tensor
#'
#' @importFrom reticulate array_reshape
#' @importFrom purrr map
#' @importFrom GenomicRanges GRanges
#'
#' @param experiment deep_tss object
#'
#' @return status encoded into tensor in slot @status_encoded
#'
#' @export

encode_status <- function(experiment) {

	ref_sample <- experiment@split_index %>%
		map(
			~ experiment@status[[experiment@settings$reference_sample]][.x] %>%
				sort %>%
				.$status %>%
				array_reshape(dim = length(.))
		)

	experiment@status_encoded <- ref_sample

	return(experiment)
}
