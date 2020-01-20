
#' Train/Test Split
#'
#' Split data into training and test set
#'
#' @importFrom dplyr mutate
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
	
	select_sample <- deep_obj@experiment

	train_index <- c(
		sample(which(select_sample$status == 1), size = ceiling(train_split / 2)),
		sample(which(select_sample$status == 0), size = ceiling(train_split / 2))
	)

	test_index <- c(
		sample(which(select_sample$status == 1 & !(select_sample$status %in% train_index)), size = ceiling(test_split / 2)),
		sample(which(select_sample$status == 0 & !(select_sample$status %in% train_index)), size = ceiling(test_split / 2))
	)
	
	select_sample <- mutate(select_sample, "index" = NA)
	select_sample[train_index, "index"] <- "train"
	select_sample[test_index, "index"] <- "test"
	
	gr <- makeGRangesFromDataFrame(select_sample, keep.extra.columns = TRUE)

	deep_obj@ranges$all <- gr
	deep_obj@experiment <- select_sample
	return(deep_obj)
}
