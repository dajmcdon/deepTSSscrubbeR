
#' Build TSS Model
#'
#' Build model of TSS Confidence
#'
#' @import keras
#'
#' @param deep_obj deep_tss object
#' @param model_type Whether to model 'status' or 'score'
#' @param optimizer Optimizer for compiling model
#' @param metric Metric to optimize model.
#' Defaults are "accuracy" for classification, and "mean_absolute_error"
#' for score regression
#'
#' @rdname tss_model
#'
#' @export

tss_model <- function(deep_obj, model_type = "status", optimizer = "adam", metric = "default") {

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
	if (model_type == "status") {

		answer <- concatenated_inputs %>%
			layer_dense(units = 512, activation = "relu") %>%
			layer_dense(units = 512, activation = "relu") %>%
			layer_dense(units = 256, activation = "relu") %>%
			layer_dense(units = 256, activation = "relu") %>%
			layer_dense(units = 1, activation = "sigmoid")

	} else if (model_type == "score") {

		answer <- concatenated_inputs %>%
			layer_dense(units = 512, activation = "relu") %>%
			layer_dense(units = 512, activation = "relu") %>%
			layer_dense(units = 256, activation = "relu") %>%
			layer_dense(units = 256, activation = "relu") %>%
			layer_dense(units = 1)

	}

	## Final model.
	model <- keras_model(
		list(genomic_input, soft_input, signal_input, shape_input),
		answer
	)

	## Compile model.
	if (model_type == "status") {
		if (metric == "default") metric <- "accuracy"

		model %>% compile(
			optimizer = optimizer,
			loss = "binary_crossentropy",
			metrics = metric
		)
	} else if (model_type == "score") {
		if (metric == "default") metric <- "mean_absolute_error"

		model %>% compile(
			optimizer = optimizer,
			loss = "mean_squared_error",
			metrics = metric
		)
	}

	deep_obj@settings$model_type <- model_type
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
#' @param subset Optional vector of read qnames to train on
#' @param weights Weights as list to pass to fit function
#'
#' @rdname tss_train-function
#'
#' @export

tss_train <- function(
		deep_obj, epochs = 25, batch_size = 150, validation_split = 0.25,
		subset = NA, weights = NA
) {
	
	## Get the training data.
	if (!is.na(subset)) {
		table_subset <- as.data.table(deep_obj@experiment)[qname %in% subset]
		tss_groups <- table_subset[, tss_group]
		soft_groups <- table_subset[, soft_group]
	} else {
		indices <- as.data.table(deep_obj@status_indices$train)
		tss_groups <- indices[, tss_group]
		soft_groups <- indices[, soft_group]
	}

	sequence_input <- deep_obj@encoded$genomic[tss_groups, , , , drop = FALSE]
	soft_input <- deep_obj@encoded$soft[soft_groups, , , , drop = FALSE]
	signal_input <- deep_obj@encoded$signal[tss_groups, , , , drop = FALSE]
	shape_input <- deep_obj@encoded$shape[tss_groups, , , , drop = FALSE]

	## Pull out the previously generated model.
	deep_model <- deep_obj@model$model
	
	## Get original scores or status.
	if (!is.na(subset)) {

		if (deep_obj@settings$model_type == "status") {
			original_values <- table_subset[, status]
			original_values <- array_reshape(original_values, dim = length(original_values))
		} else {
			original_values <- log2(table_subset[, score])
			original_values <- array_reshape(original_values, dim = length(original_values))
		}

	} else {

		if (deep_obj@settings$model_type == "status") {
			original_values <- deep_obj@encoded$status$train
		} else {
			original_values <- deep_obj@encoded$score$train
		}
	
	}

	## Train model.
	if (is.na(weights)) {

		history <- deep_model %>%
			fit(
				list(
					genomicinput = sequence_input,
					softinput = soft_input,
					signalinput = signal_input,
					shapeinput = shape_input
				),
				original_values,
				epochs = epochs,
				batch_size = batch_size,
				validation_split = validation_split
			)

	} else {

		history <- deep_model %>%
			fit(
				list(
					genomicinput = sequence_input,
					softinput = soft_input,
					signalinput = signal_input,
					shapeinput = shape_input
				),
				original_values,
				epochs = epochs,
				batch_size = batch_size,
				validation_split = validation_split,
				class_weights = weights
			)

	}

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

        if (deep_obj@settings$model_type == "status") {
                original_values = deep_obj@encoded$status$test
        } else {
                original_values = deep_obj@encoded$score$test
        }


	accuracy_results <- deep_obj@model$trained_model %>%
		evaluate(
			list(
				genomicinput = sequence_input,
				softinput = soft_input,
				signalinput = signal_input,
				shapeinput = shape_input
			),
			original_values
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

	## Get appropriate data.
	tss_groups <- deep_obj@experiment$tss_group
	soft_groups <- deep_obj@experiment$soft_group

        sequence_input <- deep_obj@encoded$genomic[tss_groups, , , , drop = FALSE]
        soft_input <- deep_obj@encoded$soft[soft_groups, , , , drop = FALSE]
        signal_input <- deep_obj@encoded$signal[tss_groups, , , , drop = FALSE]
        shape_input <- deep_obj@encoded$shape[tss_groups, , , , drop = FALSE]

	## predict.
	predicted <- deep_obj@model$trained_model %>%
		predict(list(
			genomicinput = sequence_input,
			softinput = soft_input,
			signalinput = signal_input,
			shapeinput = shape_input
		))

	results <- as.data.table(deep_obj@experiment)

	if (deep_obj@settings$model_type == "score") {
		results[, c("log2_score", "predicted_log2_score") := list(log2(score), predicted[, 1])]
	} else {
		results[, status_prob := predicted[, 1]]
	}

	deep_obj@results$all <- as.data.frame(results)
	return(deep_obj)				
}
