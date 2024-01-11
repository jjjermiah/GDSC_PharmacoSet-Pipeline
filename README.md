# This is a working repository for the GDSC Pipeline


# current pipeline status:
- no summarized experiments have been created yet.
  - rules and scripts have been created to outline the process but actual SummarizedExperiment objects have not been created yet.
- treatment data is yet to be downloaded.
  - idea is to build the MAE first as that is the same between GDSC1 and GDSC2
  - then build the GDSC1 and GDSC2 treatmentResponseExperiment in two new repos
- TODO:: annotate treatments and samples
- TODO:: extract `GRanges` object into its own rule to be used across rules
- TODO:: remove dependency of metadata file in each preprocessing rule, save for SE or build PSet
- TODO:: annotate metadata for each `Experiment` object using config details
- TODO:: molecularProfiles Create SummarizedExperiments
  - rnaseq : DONE
  - cnv : need pre-processing 
  - fusion : need pre-processing
  - mutation : pre-processing DONE
  - microarray : need pre-processing
- TODO:: create conda environments + docker images for each rule and use them in the pipeline
  

``` bash
snakemake \
  --snakefile workflow/Snakefile \
  --rulegraph | dot -Tsvg > resources/rulegraph.svg
```

### The following dag shows the pipeline steps, though the steps are not implemented for all. see TODOs above.
![pipeline status](resources/rulegraph.svg)
