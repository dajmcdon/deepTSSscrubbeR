
#' Get Softclipped Bases
#'
#' Extract softclipped bases from bam
#'
#' @include deep_tss.R
#'
#' @param deep_obj deep_tss object
#'
#' @import tibble
#' @importFrom dplyr pull mutate case_when mutate_if select add_count
#' @importFrom stringr str_extract_all str_match str_sub
#' @importFrom tidyr replace_na
#' @importFrom purrr map2_chr
#'
#' @rdname softclipped-function
#'
#' @export

softclipped <- function(deep_obj) {
	
	neg_strand <- deep_obj@experiment %>%
		pull(strand) %>%
		{. == "-"}

	soft_bases <- deep_obj@experiment %>%
		pull(cigar) %>%
		str_extract_all("(\\d+\\w)") %>%
		map2_chr(neg_strand, function(x, y) {
			if (y) x <- rev(x)
			return(x[1])
		})

	soft_bases <- deep_obj@experiment %>%
		add_column("fiveprime_cigar" = soft_bases) %>%
		mutate("fiveprime_soft" = str_match(fiveprime_cigar, "^(\\d+)S")[,2] %>% as.double) %>%
		replace_na(list(fiveprime_soft = 0)) %>%
		mutate_if(is.integer, as.double) %>%
		mutate_if(is.factor, as.character) %>%
		mutate(
			"soft_bases" = ifelse(fiveprime_soft == 0, NA, str_sub(seq, end = fiveprime_soft)),
			"start" = case_when(
				fiveprime_soft > 0 & strand == "+" ~ (bam_start + fiveprime_soft),
				fiveprime_soft > 0 & strand == "-" ~ (bam_start - fiveprime_soft),
				fiveprime_soft == 0 ~ bam_start
			),
			"end" = start
		) %>%
		select(-seq) %>%
		add_count(seqnames, strand, start, name = "score")
		
	
	deep_obj@experiment <- soft_bases
	return(deep_obj)
}
