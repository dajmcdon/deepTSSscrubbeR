
#' Train/Test Split
#'
#' Split data into training and test set
#'
#' @import tibble
#' @importFrom dplyr mutate group_by ungroup sample_n filter
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#'
#' @param deep_obj deep_tss object
#' @param train_split Total number of TSSs to use on training set, split evenly between marked status
#' @param test_split Total number of TSSs to use on test set, split evenly between marked status
#'
#' @return Training and Test sets saved into slot @split
#'
#' @export

split_data <- function(deep_obj, train_split, test_split) {
	
	## Train and test will have equal split between potential spurious and likely good TSSs.
	train_size <- ceiling(train_split / 2)
	test_size <- ceiling(test_split / 2)

	## Select randomly one of same TSS with same softclip status.
	select_TSSs <- as_tibble(deep_obj@experiment, name_repair = "unique") %>%
		filter(!is.na(status)) %>%
		group_by(tss_group) %>%
		sample_n(1) %>%
		ungroup

	## Mark TSSs for training.
	train_samples <- select_TSSs %>%
		group_by(status) %>%
		sample_n(train_size) %>%
		ungroup %>%
		mutate(index = "train") %>%
		as.data.table

	## Mark TSSs for testing.
	test_samples <- select_TSSs %>%
		filter(!tss_group %in% train_samples[, tss_group]) %>%
		group_by(status) %>%
		sample_n(test_size) %>%
		ungroup %>%
		mutate(index = "test") %>%
		as.data.table

	## Add index back to data.
	merged <- rbind(train_samples, test_samples)
	merged <- merged[, .(rowid, index)]
	setkey(merged, rowid)
	original <- as.data.table(deep_obj@experiment, key = "rowid")
	merged <- merge(original, merged, all.x = TRUE)

	## Return data back to deep_tss object
	deep_obj@experiment <- as.data.frame(merged)
	return(deep_obj)
}
