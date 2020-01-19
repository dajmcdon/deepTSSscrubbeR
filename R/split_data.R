
#' Train/Test Split
#'
#' Split data into training and test set
#'
#' @importFrom reticulate array_reshape
#' @importFrom purrr map
#' @importFrom GenomicRanges sort
#'
#' @param experiment deep_tss object
#' @param reference_sample Sample to analyze
#' @param train_split Total number of TSSs to use on training set, split evenly between marked status
#' @param test_split Total number of TSSs to use on test set, split evenly between marked status
#'
#' @return Training and Test sets saved into slot @split
#'
#' @export

split_data <- function(experiment, reference_sample, train_split, test_split) {
	select_sample <- experiment@status[[reference_sample]]
	
	train_index <- c(
		sample(which(select_sample$status == 1), size = ceiling(train_split / 2)),
		sample(which(select_sample$status == 0), size = ceiling(train_split / 2))
	)

	test_index <- c(
		sample(which(select_sample$status == 1 & !(x$status %in% train_index)), size = ceiling(test_split / 2)),
		sample(which(select_sample$status == 0 & !(x$status %in% train_index)), size = ceiling(test_split / 2))
	)
	
	split_index <- list("train_split" = train_index, "test_split" = test_index)
		
	split_list <- list(
		"sequence" = list(
			"train_split" = sort(experiment@ranges$sequence[[reference_sample]][train_index]),
			"test_split" = sort(experiment@ranges$sequence[[reference_sample]][test_index])
		),
		"signal" = list(
			"train_split" = sort(experiment@ranges$signal[[reference_sample]][train_index]),
			"test_split" = sort(experiment@ranges$signal[[reference_sample]][test_index])
		)
	)

	experiment@split_index <- split_index
	experiment@split <- split_list
	experiment@settings$reference_sample <- reference_sample
	experiment@settings$train_split <- train_split
	experiment@settings$test_split <- test_split

	return(experiment)
}
