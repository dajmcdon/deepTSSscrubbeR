# deepTSSscrubbeR

Provide confidence levels to TSSs as not being artifacts.

## About

Template Switching Reverse Transcription (TSRT) based 5' mapping technologies such as RAMPAGE, STRIPE-seq, and nanoCAGE
have gained favor because of their low input requirments and/or relative simplicity.
However, certain artifacts tend to be somewhat common in the results of these methods.
These artifacts can include strand invasion (invasion of the template switching primer into the extending cDNA),
template switching at uncapped ends (although with less probability than capped ends),
and the addition of an extra cytidine to the 5' end of the cDNA (common to even orthogonal techniques).

This software aims to reduce the prevalence of artifactual TSSs by using Convolution Neural Networks (CNNs) to model 
characteristics of spurious TSSs, and then comparing these characteristics to all TSSs.
The CNN takes into consideration a few factors when assessing the validity of a TSS.
First, in most 5' mapping technologies, there is an almost universal prevalence of the addition 
of non-templated cytosines just adjacent to the TSS potentially due to the cap itself acting as a template for reverse transcription.
Second, most organisms tend to have a sequence bias at TSSs, such as pyrimidine purine bias in humans and yeast.
Third, most TSSs tend to cluster into Transcription Start Regions (TSRs), so sparse singeltons tend to be less common at true TSSs.
Finally, true promoters tend to have certain DNA structural properties.

The software uses a bam file as input, and should still include the bases that were softclipped adjacent to the TSS by the alignment software.

## Instlaling Software

It is recommended to use the provided singularity container with all required software installed.
Singularity are containers similar to docker containers that allow compatability and reproducibilty for software and workflows.
You must first install the [singularity software](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) 
onto your machine to use containers.



## Quick Start

Pull the singularity container, shell into it, activate the conda environment,
then launch R.
```
singularity pull library://rpolicastro/default/deep_tss_scrubber:0.7.0
singularity shell -eCB "$(pwd)" -H "$(pwd)" deep_tss_scrubber_0.7.0.sif

. /opt/conda/etc/profile.d/conda.sh
conda activate r-reticulate; R
```
You can now start using the deepTSSscrubbeR library.

```
library("deepTSSscrubbeR")
library("magrittr")
## 'dummyVars' seems to be broken unless you load the caret library
library("caret")

bam <- system.file("extdata", "S288C.bam", package = "deepTSSscrubbeR")
assembly_fasta <- system.file("extdata", "yeast_assembly.fasta", package = "deepTSSscrubbeR")
annotation_gtf <- system.file("extdata", "yeast_annotation.gtf", package = "deepTSSscrubbeR")

tss_obj <- deep_tss(bam) %>%
	tss_annotate(annotation = annotation_gtf, upstream = 250, downstream = 1) %>%
	get_softclip %>%
	mark_status(use_annotation = TRUE) %>%
	split_data(train_split = 7500, test_split = 2500) %>%
	expand_ranges %>%
	get_sequences(assembly_fasta) %>%
	get_signal

export_raw(tss_obj, "sequence_data.tsv", "signal_data.tsv")

tss_encoded <- tss_obj %>%
	encode_genomic %>%
	encode_shape %>%
	encode_soft %>%
	encode_signal %>%
	encode_status

export_encoded(tss_encoded)

# Classifier

deep_model <- tss_encoded %>%
	tss_model(metric = c("AUC", "accuracy")) %>%
	tss_train(epochs = 10, batch_size = 250) %>%
	tss_evaluate %>%
	tss_predict

export_bedgraphs(deep_model)
export_report(deep_model)

# Regression

deep_model <- tss_encoded %>%
	tss_model(model_type = "score") %>%
	tss_train(epochs = 20) %>%
	tss_evaluate %>%
	tss_predict

export_bedgraphs(deep_model)

# Training on selected TSSs.

qnames <- system.file("extdata", "stable_set.tsv", package = "deepTSSscrubbeR") %>%
	read.delim(header = TRUE, stringsAsFactors = FALSE, sep = "\t")
qnames <- qnames[, 1]

deep_model <- tss_encoded %>%
	tss_model(metric = c("AUC", "accuracy")) %>%
	tss_train(
		validation_split = 0.10,
		epochs = 20,
		batch_size = 10000,
		subset = qnames,
		weights = list("0" = 1, "1" = 101232 / 6931)
	) %>%
	tss_predict

```

## Detailed

**WORK IN PROGRESS**

Read in the TSSs

```
library("reticulate")
use_condaenv("keras")

bam <- system.file("extdata", "S288C.bam", package = "deepTSSscrubbeR")
```

Make a deepTSSscrubbeR object

```
tss_obj <- deep_tss(bam)
```

Analyze 5' ends for soft clipped bases

```
tss_obj <- get_softclip(tss_obj)
```

Mark likely spurious TSSs based on number of reads

```
tss_obj <- mark_status(tss_obj, lower = 2, upper = 5)
```

split data into training and test sets

```
tss_obj <- split_data(tss_obj, train_split = 1000, test_split = 1000)
```

Expand ranges for downstream analysis

```
tss_obj <- expand_ranges(tss_obj,
	sequence_expansion = 10,
	signal_expansion = 10,
	shape_expansion = 10
)
```

Get genomic sequences for ranges

```
assembly_fasta <- system.file("extdata", "yeast_assembly.fasta", package = "deepTSSscrubbeR") 

tss_obj <- get_sequences(tss_obj, assembly_fasta)
```

Get surrounding signal

```
tss_obj <- get_signal(tss_obj)
```

Encode the genomic sequences

```
tss_obj <- encode_genomic(tss_obj)
```

Encode soft-clipped bases

```
tss_obj <- encode_soft(tss_obj)
```

Encode TSS status

```
tss_obj <- encode_status(tss_obj)
```

Encode signal around TSS

```
tss_obj <- encode_signal(tss_obj)
```

.
.
.
