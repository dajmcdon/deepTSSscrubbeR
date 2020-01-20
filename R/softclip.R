
#' Get Softclipped Bases
#'
#' Extract softclipped bases from bam
#'
#' @import tibble
#' @importFrom dplyr pull bind_rows
#' @importFrom stringr str_extract_all
#'
#' @rdname softclipped-function
#'
#' @export

softclipped <- function(deep_obj) {
	soft_bases <- deep_obj@experiment %>%
		pull(cigar) %>%
		str_extract_all("(\\d+\\w)") %>%
		map_df(~ matrix(., nrow = 1, byrow = TRUE) %>% as_tibble)

	return(soft_bases)
}
