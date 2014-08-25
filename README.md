3dgenome
========

Scripts to recreate figures in paper: _Integrative modelling reveals the principles of multi-scale chromatin boundary formation in human nuclear organisation_ (2014, in preparation).

Instructions on how to run, dependencies etc. will follow...

# Dependencies

## blmR

These scripts rely on an ad-hoc R package I've written called "blmR". To install, use the `devtools` package:

```r
devtools::install_github("blmR", "blmoore")
```

## External programs

* `bedtools`
* `bigWigAverageOverBed`^
* `ICE`^

^ To generate new input files, otherwise not required
