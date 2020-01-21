
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
	genomic_length <- deep_obj@settings$sequence_expansion * 2
	genomic_input <- layer_input(shape = c(6, genomic_length, 1), name = "genomicinput")

	# Surrounding signal input.
	signal_length <- deep_obj@settings$signal_expansion * 2
	signal_input <- layer_input(shape = c(1, signal_length, 1), name = "signalinput")

	# Softclipped input.
	soft_input <- layer_input(shape = c(6, 3, 1), name = "softinput")

	## Specify input layers.
	
	# Genomic sequence layer.
	genomic_layer <- genomic_input %>%
		layer_conv_2d(filter = 32, kernel_size = c(2, 3), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_conv_2d(filter = 64, kernel_size = c(2, 3), activation = "relu", padding = "same") %>%
		layer_dropout(0.25) %>%
		layer_flatten()

	# Softclipped base layer.
	soft_layer <- soft_input %>%
		layer_conv_2d(filter = 32, kernel_size = c(2, 1), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_conv_2d(filter = 64, kernel_size = c(2, 2), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_flatten()

	# Surrounding signal layer.
	signal_layer <- signal_input %>%
		layer_conv_2d(filter = 32, kernel_size = c(1, 2), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_conv_2d(filter = 64, kernel_size = c(1, 4), activation = "relu") %>%
		layer_dropout(0.25) %>%
		layer_flatten()

	## Concatenate input layers.
	concatenated_inputs <- layer_concatenate(list(
		genomic_layer,
		soft_layer,
		signal_layer
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
		list(genomic_input, soft_input, signal_input),
		answer
	)

	## Compile model.
	model %>% compile(
		optimizer = optimizer,
		loss = "binary_crossentropy",
		metrics = metrics
	)
}
