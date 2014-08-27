3dgenome
========

Scripts to recreate figures in paper: _Integrative modelling reveals the principles of multi-scale chromatin boundary formation in human nuclear organisation_ (2014, in preparation).

Full instructions on how to run, input files etc. will follow...

# Dependencies

## blmR

These scripts rely on an _ad hoc_ R package I've written called "blmR". To install, use the `devtools` package:

```r
devtools::install_github("blmR", "blmoore")
```

## R packages

Lots of commonly-installed R packages are also used, including but not limited to: 

### CRAN

* `caret`
* `dplyr` 
* `ggplot2`
* `gridExtra`
* `plotrix`
* `randomForest`
* `RColorBrewer`
* `reshape2`
* `RHmm`

### Bioconductor

* `BSgenome` (and UCSC hg19)
* `GenomicRanges`

## External programs

* `bedtools`
* `bigWigAverageOverBed`<sup>*</sup>
* `ICE`<sup>*</sup>

<sup>*</sup> To generate new input files, otherwise not required

## Raw data

Serialised rds data files are provided under [`data/rds`](data/rds). These were built from ENCODE uniformly processed .bigWig files which are available to download [from ENCODE](http://encodedcc.sdsc.edu/ftp/modENCODE_VS_ENCODE/Regulation/Human/signal/foldChange/) (list of filenames used [here](data/ENCODE_filelist.txt)). These files (used in [Boyle _et al._ (2014)](http://www.nature.com/nature/journal/v512/n7515/full/nature13668.html)) come from the July 2012 data freeze and use MACSv2 to convert aligned ChIP-seq reads into a measure of "signal" relative to input chromatin background.

Some of the scripts rely on other large files not included in this repository (in order to reanalyse data from scratch), but reasonable intermediates are provided where possible and stored under [`data/`](data/). For notes on which scripts are effected, see the [how to run](https://github.com/blmoore/3dgenome#how-to-run) section of this README.

## How to run


## sessionInfo()

Below is an output of sessionInfo() for troubleshooting purposes, some loaded packages may not be required and likewise, some required packages may not be loaded.

```r
R version 3.1.1 (2014-07-10)
Platform: x86_64-apple-darwin13.1.0 (64-bit)

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Rmpi_0.6-5           snow_0.3-13          pvclust_1.2-2        naturalsort_0.1.2    GenomicRanges_1.16.4 GenomeInfoDb_1.0.2  
 [7] IRanges_1.22.10      BiocGenerics_0.10.0  dplyr_0.2            data.table_1.9.2     RHmm_2.0.3           nlme_3.1-117        
[13] plyr_1.8.1           blmR_0.1             RColorBrewer_1.0-5   infotheo_1.2.0       ROCR_1.0-5           gplots_2.14.1       
[19] calibrate_1.7.2      randomForest_4.6-10  caret_6.0-30         lattice_0.20-29      gridExtra_0.9.1      ggplot2_1.0.0       
[25] fGarch_3010.82       fBasics_3010.86      MASS_7.3-33          timeSeries_3010.97   timeDate_3010.98     corrgram_1.5        
[31] seriation_1.0-13     reshape2_1.4         R.matlab_3.0.1       plotrix_3.5-7        scales_0.2.4        

loaded via a namespace (and not attached):
 [1] assertthat_0.1      bitops_1.0-6        BradleyTerry2_1.0-5 brglm_0.5-9         car_2.0-20          caTools_1.17       
 [7] cluster_1.15.2      codetools_0.2-8     colorspace_1.2-4    devtools_1.5        digest_0.6.4        evaluate_0.5.5     
[13] foreach_1.4.2       gclus_1.3.1         gdata_2.13.3        gtable_0.1.2        gtools_3.4.1        httr_0.4           
[19] iterators_1.0.7     KernSmooth_2.23-12  labeling_0.2        lme4_1.1-7          magrittr_1.0.1      Matrix_1.1-4       
[25] memoise_0.2.1       minqa_1.2.3         munsell_0.4.2       nloptr_1.0.4        nnet_7.3-8          proto_0.3-10       
[31] R.methodsS3_1.6.1   R.oo_1.18.0         R.utils_1.32.4      Rcpp_0.11.2         RCurl_1.95-4.3      roxygen2_4.0.1.99  
[37] splines_3.1.1       stabledist_0.6-6    stats4_3.1.1        stringr_0.6.2       tools_3.1.1         TSP_1.0-9          
[43] whisker_0.3-2       XVector_0.4.0   
```
