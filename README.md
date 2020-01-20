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
tss_obj <- softclipped(tss_obj)
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
tss_obj <- expand_ranges(tss_obj, sequence_expansion = 10, signal_expansion = 10)
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
