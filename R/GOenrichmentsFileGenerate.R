## Just generates genomic ranges in format for use with ensembl BioMart, 
## i.e. to get gene IDs within these regions for enrichment testing (which
## is done by the DAVID online tool: http://david.abcc.ncifcrf.gov/)

# suppress sci notation
options(scipen=999)

## All regions:
capture.output(
  cat(paste0(gsub("chr", "", paste0(gsub("-", ":", rownames(h.dat)), ":", 
       as.numeric(gsub(".*-", "", rownames(h.dat))) + 1e6)), collapse=",")),
  file=file("~/hvl/ice/bedfiles/allRegionsEns.txt", "w"))

g.open <- read.table("~/hvl/ice/bedfiles/g.open.bed")
capture.output(
  cat(paste0(gsub("chr", "", g.open[,1]), ":", g.open[,2], ":", g.open[,3], collapse=",")),
  file=file("~/hvl/ice/bedfiles/gopen_biomartRegions.txt", "w"))

h.open <- read.table("~/hvl/ice/bedfiles/h.open.bed")
capture.output(
  cat(paste0(gsub("chr", "", h.open[,1]), ":", h.open[,2], ":", h.open[,3], collapse=",")),
  file=file("~/hvl/ice/bedfiles/hopen_biomartRegions.txt", "w"))

k.open <- read.table("~/hvl/ice/bedfiles/k.open.bed")
capture.output(
  cat(paste0(gsub("chr", "", k.open[,1]), ":", k.open[,2], ":", k.open[,3], collapse=",")),
  file=file("~/hvl/ice/bedfiles/kopen_biomartRegions.txt", "w"))
