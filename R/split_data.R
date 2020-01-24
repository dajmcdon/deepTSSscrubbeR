
#' Train/Test Split
#'
#' Split data into training and test set
#'
#' @importFrom dplyr mutate group_by ungroup sample_n bind_rows left_join filter select
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
	select_TSSs <- deep_obj@experiment %>%
		group_by(tss_group) %>%
		sample_n(1) %>%
		ungroup %>%
		filter(!is.na(status))

	## Mark TSSs for training.
	train_samples <- select_TSSs %>%
		group_by(status) %>%
		sample_n(train_size) %>%
		ungroup %>%
		mutate(index = "train")

	## Mark TSSs for testing.
	test_samples <- select_TSSs %>%
		filter(!tss_group %in% pull(train_samples, tss_group)) %>%
		group_by(status) %>%
		sample_n(test_size) %>%
		ungroup %>%
		mutate(index = "test")

	## Add index back to data.
	merged <- bind_rows(train_samples, test_samples) %>%
		select(qname, index) %>%
		left_join(deep_obj@experiment, ., by = "qname")
	
	gr <- makeGRangesFromDataFrame(merged, keep.extra.columns = TRUE)
	gr_train <- gr[!is.na(gr$index) & gr$index == "train"]
	gr_test <- gr[!is.na(gr$index) & gr$index == "test"]

	deep_obj@ranges$original$all <- gr
	deep_obj@ranges$original$train <- gr_train
	deep_obj@ranges$original$test <- gr_test
	deep_obj@experiment <- merged

	return(deep_obj)
}
