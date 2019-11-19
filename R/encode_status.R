
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
	ref_sample <- experiment@split[[reference_sample]] %>%
		map(function(x) {
			status_tensor <- array_reshape(x$status, dim = length(x$status))
			return(status_tensor)
		})

	experiment@status_encoded <- ref_sample
	experiment@settings$reference_sample <- reference_sample

	return(experiment)
}
