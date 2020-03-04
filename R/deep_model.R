
#' Build TSS Model
#'
#' Build model of TSS Confidence
#'
#' @import keras
#'
#' @param deep_obj depp_tss object
#' @param optimizer Optimizer for compiling model
#' @param metric Metric to optimize model
#'
#' @rdname tss_model
#'
#' @export

tss_model <- function(deep_obj, optimizer = "adam", metric = "accuracy") {

	## Specify input sizes.

	# Genomic sequence input.
	genomic_length <- (deep_obj@settings$sequence_expansion * 2) + 1
	genomic_input <- layer_input(shape = c(genomic_length, 5, 1), name = "genomicinput")

	# Surrounding signal input.
	signal_length <- (deep_obj@settings$signal_expansion * 2) + 1
	signal_input <- layer_input(shape = c(1, signal_length, 1), name = "signalinput")

	# Softclipped input.
	soft_input <- layer_input(shape = c(3, 6, 1), name = "softinput")

	# DNA shape input.
	shape_length <- dim(deep_obj@encoded$shape)[3]
	shape_input <- layer_input(shape = c(1, shape_length, 1), name = "shapeinput")

	## Specify input layers.
	
	# Genomic sequence layer.
	genomic_layer <- genomic_input %>%
		layer_conv_2d(filters = 32, kernel_size = c(3, 2), activation = "relu") %>%
		layer_dropout(0.5) %>%
		layer_conv_2d(filters = 64, kernel_size = c(3, 2), activation = "relu", padding = "same") %>%
		layer_dropout(0.5) %>%
		layer_flatten()

	# Softclipped base layer.
	soft_layer <- soft_input %>%
		layer_conv_2d(filters = 32, kernel_size = c(1, 2), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_conv_2d(filters = 64, kernel_size = c(2, 2), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_flatten()

	# Surrounding signal layer.
	signal_layer <- signal_input %>%
		layer_conv_2d(filters = 32, kernel_size = c(1, 2), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_conv_2d(filters = 64, kernel_size = c(1, 4), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_flatten()

	# DNA shape layer.
	shape_layer <- shape_input %>%
		layer_conv_2d(filters = 32, kernel_size = c(1, 3), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_conv_2d(filters = 64, kernel_size = c(1, 5), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_flatten()

	## Concatenate input layers.
	concatenated_inputs <- layer_concatenate(list(
		genomic_layer,
		soft_layer,
		signal_layer,
		shape_layer
	))

	## Answer layer.
	answer <- concatenated_inputs %>%
		layer_dense(units = 512, activation = "relu") %>%
		layer_dense(units = 512, activation = "relu") %>%
		layer_dense(units = 256, activation = "relu") %>%
		layer_dense(units = 256, activation = "relu") %>%
		layer_dense(units = 1, activation = "sigmoid")

	## Final model.
	model <- keras_model(
		list(genomic_input, soft_input, signal_input, shape_input),
		answer
	)

	## Compile model.
	model %>% compile(
		optimizer = optimizer,
		loss = "binary_crossentropy",
		metrics = metric
	)

	deep_obj@model$model <- model
	return(deep_obj)
}

#' Train TSS Model
#'
#' Train the deepTSSscrubbeR model
#'
#' @import keras
#'
#' @param deep_obj tss_obj
#' @param epochs epochs
#' @param batch_size batch size
#' @param validation_split validation split
#'
#' @rdname tss_train-function
#'
#' @export

tss_train <- function(deep_obj, epochs = 25, batch_size = 150, validation_split = 0.25) {
	
	## Get the training data.
	indices <- as.data.table(deep_obj@status_indices$train)
	tss_groups <- indices[, tss_group]
	soft_groups <- indices[, soft_group]

	sequence_input <- deep_obj@encoded$genomic[tss_groups, , , , drop = FALSE]
	soft_input <- deep_obj@encoded$soft[soft_groups, , , , drop = FALSE]
	signal_input <- deep_obj@encoded$signal[tss_groups, , , , drop = FALSE]
	shape_input <- deep_obj@encoded$shape[tss_groups, , , , drop = FALSE]

	## Pull out the previously generated model.
	deep_model <- deep_obj@model$model
	
	## Train model.
	history <- deep_model %>%
		fit(
			list(
				genomicinput = sequence_input,
				softinput = soft_input,
				signalinput = signal_input,
				shapeinput = shape_input
			),
			deep_obj@encoded$status$train,
			epochs = epochs,
			batch_size = batch_size,
			validation_split = validation_split
		)

	## Add trained model to deep tss object.
	deep_obj@model$train_history <- history
	deep_obj@model$trained_model <- deep_model

	return(deep_obj)
}

#' Evaluate TSS Model
#'
#' Evaluate the trained deepTSSScrubbeR model
#'
#' @import keras
#'
#' @param deep_obj tss_obj with trained model
#'
#' @rdname tss_evaluate-function
#'
#' @export



tss_evaluate <- function(deep_obj) {

	## Get the test data.
        indices <- as.data.table(deep_obj@status_indices$test)
        tss_groups <- indices[, tss_group]
        soft_groups <- indices[, soft_group]

        sequence_input <- deep_obj@encoded$genomic[tss_groups, , , , drop = FALSE]
        soft_input <- deep_obj@encoded$soft[soft_groups, , , , drop = FALSE]
        signal_input <- deep_obj@encoded$signal[tss_groups, , , , drop = FALSE]
        shape_input <- deep_obj@encoded$shape[tss_groups, , , , drop = FALSE]

	accuracy_results <- deep_obj@model$trained_model %>%
		evaluate(
			list(
				genomicinput = sequence_input,
				softinput = soft_input,
				signalinput = signal_input,
				shapeinput = shape_input
			),
			deep_obj@encoded$status$test
		)

	deep_obj@model$evaluation <- accuracy_results
	return(deep_obj)
}

#' Predict TSS artifact probabilities
#'
#' Predict the probability of being an artifact based on the deepTSSscrubbeR model
#'
#' @import keras
#' @import tibble
#' @importFrom GenomicRanges GRanges
#'
#' @param deep_obj deep_tss object with trained model
#'
#' @rdname tss_predict-function
#'
#' @export

tss_predict <- function(deep_obj) {
	predicted_probs <- deep_obj@model$trained_model %>%
		predict(list(
			genomicinput = deep_obj@encoded$genomic$all,
			softinput = deep_obj@encoded$softclipped$all,
			signalinput = deep_obj@encoded$signal$all,
			shapeinput = deep_obj@encoded$shape$all
		))

	predictions <- deep_obj@ranges$sequence$all
	predictions$probs <- predicted_probs

	deep_obj@results$all <- predictions
	return(deep_obj)				
}
