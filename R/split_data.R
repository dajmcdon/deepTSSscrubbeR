
#' Train/Test Split
#'
#' Split data into training and test set
#'
#' @importFrom reticulate array_reshape
#' @importFrom purrr map
#'
#' @param experiment deep_tss object
#' @param train_split Total number of TSSs to use on training set, split evenly between marked status
#' @param test_split Total number of TSSs to use on test set, split evenly between marked status
#'
#' @return Training and Test sets saved into slot @split
#'
#' @export

split_data <- function(experiment, train_split, test_split) {
	split <- experiment@status %>%
		map(function(x) {
			train_index <- c(
				sample(which(x$status == 1), size = ceiling(train_split / 2)),
				sample(which(x$status == 0), size = ceiling(train_split / 2))
			)

			test_index <- c(
				sample(which(x$status == 1 & !(x$status %in% train_index)), size = ceiling(test_split / 2)),
				sample(which(x$status == 0 & !(x$status %in% train_index)), size = ceiling(test_split / 2))
			)
			
			split_list <- list(
				"train_split" = x[train_index],
				"test_split" = x[test_index]
			)

			return(split_list)
		})

	experiment@split <- split
	return(experiment)
}
