Read in the TSSs

```
TSSs <- system.file("extdata", "TSSs.RDS", package = "deepTSSscrubbeR")
TSSs <- readRDS(TSSs)
```

Make a deepTSSscrubbeR object

```
tss_obj <- deep_tss(TSSs)
```

Expand ranges for downstream analysis

```
tss_obj <- expand_ranges(tss_obj, 10, 15)
```
