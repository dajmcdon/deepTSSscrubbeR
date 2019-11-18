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
tss_obj <- expand_ranges(tss_obj, 10, 15)
```
