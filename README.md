Read in the TSSs

```
TSSs <- system.file("extdata", "TSSs.RDS", package = "deepTSSscrubbeR")
TSSs <- readRDS(TSSs)
```

Make a deepTSSscrubbeR object

```
tss_obj <- deep_tss(TSSs)
```

Mark likely spurious TSSs based on number of reads

```
tss_obj <- mark_status(tss_obj, threshold = 2)
```

Expand ranges for downstream analysis

```
tss_obj <- expand_ranges(tss_obj, sequence_expansion = 10, signal_expansion = 15)
```

Split data into training and test sets

```
tss_obj <- split_data(tss_obj, train_split = 1000, test_split = 1000)

Set reference set to analyze, and encode the TSS status of that set

```
tss_obj <- encode_status(tss_obj, reference_sample = "set_1")
```

Get sequences using expanded GRanges, and then one-hot encode them

```
genome_assembly <- system.file("extdata", "yeast_assembly.fasta", package = "deepTSSscrubbeR")
tss_obj <- retrieve_sequences(tss_obj, assembly = genome_assembly)
```
