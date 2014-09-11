3dgenome
========

Scripts to recreate figures in paper: _Integrative modelling reveals the principles of multi-scale chromatin boundary formation in human nuclear organisation_ (2014, in preparation).

### Table of contents

* [Dependencies](https://github.com/blmoore/3dgenome#dependencies)
* [Raw data](https://github.com/blmoore/3dgenome#raw-data)
* [How to run](https://github.com/blmoore/3dgenome#how-to-run)
* [sessionInfo](https://github.com/blmoore/3dgenome#sessioninfo)

# Dependencies

Scripts will run on OS X and other Unix-based systems. External dependencies should be installed somewhere on your `$PATH`.

## blmR

These scripts rely on an _ad hoc_ R package I've written called "blmR". To install, use the `devtools` package:

```r
devtools::install_github("blmR", "blmoore")
```

## R packages

Lots of commonly-installed R packages are also used, including but not limited to: 

### CRAN

* `caret`
* `corrgram`
* `dplyr` 
* `ggplot2` (w/ dependencies `reshape2`, `scales`, `RColorBrewer` etc.)
* `gridExtra`
* `naturalsort`
* `plotrix`
* `pvclust`
* `randomForest`
* `RHmm`
* `snow` (and `Rmpi`)

### Bioconductor

* `BSgenome` (and UCSC hg19)
* `GenomicRanges`

## External programs

* `bedtools`
* `bigWigAverageOverBed`<sup>*</sup>
* `ICE` / [hiclib](http://mirnylab.bitbucket.org/hiclib)<sup>*</sup>

<sup>*</sup> To generate new input files, otherwise not required

# Raw data

Serialised rds data files are provided under [`data/rds`](data/rds). These were built from ENCODE uniformly processed .bigWig files which are available to download [from ENCODE](http://encodedcc.sdsc.edu/ftp/modENCODE_VS_ENCODE/Regulation/Human/signal/foldChange/) (list of filenames used [here](data/text/ENCODE_filelist.txt)). These files (used in [Boyle _et al._ (2014)](http://www.nature.com/nature/journal/v512/n7515/full/nature13668.html)) come from the July 2012 data freeze and use MACSv2 to convert aligned ChIP-seq reads into a measure of "signal" relative to input chromatin background.

Some of the scripts rely on other large files not included in this repository (in order to reanalyse data from scratch), but reasonable intermediates are provided where possible and stored under [`data/`](data/). For notes on which scripts are effected, see the [how to run](https://github.com/blmoore/3dgenome#how-to-run) section of this README.

# How to run

First clone the repository (warning: ~250 MB total) and enter the directory:

```bash
git clone https://github.com/blmoore/3dgenome
cd 3dgenome
```

Now scripts can be run with Rscript, i.e.:

```bash
> Rscript R/0_buildDatFiles.R 

# output follows...

Use all variables is set to:  FALSE 

Reading in files... 

[1] "HaibTfbsH1hescAtf3V0416102AlnRep0"
[1] "SydhTfbsH1hescCebpbIggrabAlnRep0"
[1] "BroadHistoneH1hescChd1a301218aStdAlnRep0"
[1] "SydhTfbsH1hescChd2IggrabAlnRep0"
[1] "OpenChromChipH1hescCmycAlnRep1"
[1] ...


Building Random Forest...

     |      Out-of-bag   |
Tree |      MSE  %Var(y) |
  25 | 0.005201    35.61 |
  50 | 0.004946    33.87 |
  75 | 0.004833    33.09 |
 100 | 0.004842    33.15 |
 ...
```

Scripts are numbered in the order of which they should be run. You may wish to delete existing PDFs under `figures/` to convince yourself they're being regenerated (if not they'll be overwritten).

#### 0_buildDatFiles.R

Builds a tidy `data.frame` in the format: {**Y**, **x**_1, **x**_2, ..., **x**_n} where **Y** is a vector of values to be predicted (compartment eigenvector profiles in this instance) and each **x** is a numeric vector input feature. Here each cell represents an averaged feature (fold-change relative to input chromatin) over 1 Mb.

This script generates two sets of `.rds` data files under [`data/`](data/rds):

* `cellType_XVars.rds` - the aforementioned `data.frame` containing input features and empirical independent variable data.
* `cellType_XVars_RFmod.rds` - the Random Forest regression model for the respective `data.frame`.

#### 1_fig1CompartmentProfiles.R

Plots genome-wide compartment profiles for each cell type under study. The resulting figures make up part of **figure 1** ([`figures/f1_chr2.pdf`](figures/f1_chr2.pdf)) and supplementary figure 1 ([`figures/suppl/s1_GenomewideWigglePlots.pdf`](figures/suppl/s1_GenomewideWigglePlots.pdf)).

#### 2_fig1TADsDiagram.R

Generates detailed view of a ~22 Mb region of chromosome 2, indicating conservation of higher order genome organisation at multiple levels. Plots [`figures/f1_ZoomedRegion.pdf`](figures/f1_ZoomedRegion.pdf) which makes up the second half of **figure 1** in the manuscript.

#### 2b_suppl2boundaryDistances.R

Compares boundaries for compartments and TADs across cell types and tests the significance of observed relationships. Plots supplementary figure 2 in two parts: [`2a`](figures/suppl/s2a_boundsEcdf.pdf) and [`2b`](figures/suppl/s2b_compartmentCorrgram.pdf).

#### 3_fig2ModelResults.R

Plots Random Forest modelling results per cell type as three separate plots ([`f2_gmRes.pdf`](figures/f2_gmRes.pdf), [`f2_h1Res.pdf`](figures/f2_h1Res.pdf) and [`f2_k5Res.pdf`](figures/f2_gmRes.pdf)). Also calculates variable importance per model and plots a summary of the top ten ([`f2b_varImpPerModel.pdf`](figures/f2b_varImpPerModel.pdf)). Combining these plots gives **figure 2** in the manuscript.

#### 4_fig3CrossApplication.R

Performs cross-application of cell type specific models and generates a summary plot ([`f3b_crossApplyBars.pdf`](figures/f3b_crossApplyBars.pdf)). Also plots a reciprocal example of cross application between two of the cell types ([`f3ai`](figures/f3ai_gmXapplyK5.pdf) and [`f3aii`](figures/f3aii_k5XapplyGm.pdf)). Together these are **figure 3** in the manuscript.

#### 5_fig4aStratifyByVariability.R

This script has a non-neglible runtime (~10 mins on a modern processor) but could easily be parallelised / optimised. Here we're splitting the genome into an equal number of bins based on how variable chromatin structure is across the cell types under study. Then a model (as in script 3) is built per split, and the results compared ([`f4a`](figures/f4a_stratByVar.pdf)).

#### 6_fig4bFlippedBoxplots.R

**Requires external files**. To run this script you need to download the ENCODE predicted chromatin states from [here](http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/segmentations/jan2011/Combined_7_state/) (i.e. gzip archives for the three cell types being used). They should be unzipped and placed under `data/bedfiles/`, e.g.

```bash
# download bedfiles from ENCODE ftp (35 MB total)
wget -P data/bedfiles/ -i data/text/chromStateFiles.txt
# un-gzip archives
gunzip data/bedfiles/*.gz
```

Also requires R package `RHmm`, which at the time of writing has been removed from CRAN. Download the newest available version from [the archive](http://cran.r-project.org/src/contrib/Archive/RHmm/) and install with e.g. `R CMD INSTALL <RHmm_2.0.0.tar.gz>`.

Finally, this script also has the external dependancy **bedtools** which should be on your `$PATH`. Should be easy enough to install through your linux central repo (e.g. `apt-get install bedtools`) or via your OS X pkg manager of choice (`brew install bedtools`). bedtools may throw an "out of memory" exception if running on a machine with < approx. 4 GB RAM.

The output of this script includes [`f4b`](figures/f4b_enhancerEnrichFlipped.pdf) as well as supplementary figures [`s6`](figures/suppl/s6_transcribedFlipped.pdf) and [`s7`](figures/suppl/s7_allBeans.pdf). Combine f4b with f4a (script 5) and you'll get **figure 4**.

#### 7_fig5boundaryEnrichments.R

Reproducing figure 5 requires a lot of leg work, I am happy to send you intermediate files if you wish to run this script (ping me [on twitter](https://twitter.com/benjaminlmoore)). Comments within the Rscript explain what's needed but the essential steps are: 

1. Call compartment and TAD boundaries (or use those provided in [`data/bedfiles`](data/bedfiles))
2. Generate a series of equally spaced bins around each one (Python script: [binAroundBed.py](py/binAroundBed.py))
3. Use `bigWigAverageOverBed` to average all input features (~300 .bigWig files, maybe 200 GB of raw data) — ideally use a cluster for this

This gives you a series of text files with intervals and averaged signal per bin, with which you can then test for enrichments per boundary and plot **figure 5** ([f5a](figures/f5a_boundaryEnrichmentProfiles.pdf), [f5b](figures/f5b_boundaryEnrichmentBubble.pdf)) as well as supplementary figures ([s9](figures/suppl/s9_tadBoundaries.pdf) and [s10](figures/suppl/s10_compartmentBoundaries.pdf)).

#### 8_suppl10FeatureHeatmaps.R

This generates a single figure in a resource-intensive, uninteresting way (probably a good one to skip). As written uses the `pvclust` R package to find "significant" clusters in the input features through multiscale bootstrap (if a cluster is stable through undersampled and oversample data, they reason, it's significantly more clustery than would be expected by chance --- give or take, read the docs). Parallelisation is achieved with `snow` and `Rmpi` running on an openmpi backend (4 threads), but by all means run it single-threaded if you so desire. In the final(ish) paper, these heatmaps became supplementary figure 4 (not 10 as the title suggests, TO FIX).

#### 9_additionalFigures.R

Finally this script will plot the remaining supplementary figures: [`s3`](figures/suppl/s3_VarImpDifferences.pdf) and [`s5`](figures/suppl/s5_featureBoxplots.pdf). The only figure not regenerated through these scripts is supplementary figure 8 (MS numbering), which is [UCSC browser](http://genome.ucsc.edu/cgi-bin/hgTracks) screenshots, laboriously generated by hand (well, ish).

## sessionInfo()

Below is an output of sessionInfo() for troubleshooting purposes, some loaded packages may not be required and likewise, some required packages may not be loaded. An exception caused by attached packages is likely due to version issues.

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
[25] fGarch_3010.82       fBasics_3010.86      MASS_7.3-33          timeSeries_3010.97   timeDate_3010.98     corrgram_1.6        
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
