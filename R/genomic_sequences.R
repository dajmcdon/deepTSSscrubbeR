#' Retrieve Sequences
#'
#' Retrieve sequences for expanded ranges
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsamtools FaFile
#' @importFrom Biostrings getSeq
#' @importFrom purrr set_names map map2
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
                map(
                        ~ getSeq(genome, .) %>%
                        as.character %>%
                        set_names(NULL)
                ) %>%
                map2(deep_obj@ranges$sequence, function(x, y) {
                        y$seq <- x
                        return(y)
                })

        deep_obj@ranges$sequence$seq <- seqs
        return(deep_obj)
}
