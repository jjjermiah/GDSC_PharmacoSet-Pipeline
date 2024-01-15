#' RULE: preprocessMetadata
#' AUTHOR: Jermiah Joseph
#' CREATED: 01-08-2024
#' This script takes in the following files:
#' - INPUT$sampleMetadata
#' - INPUT$treatmentMetadata
#' - INPUT$cellModelPassportsMetadata
#' - INPUT$cellmodelpassportsGeneAnnotation
#' - INPUT$gencode_annotation_file
#' and outputs the following files:
#' - OUTPUT$metadata
#' 
#' PACKAGE DEPENDENCIES:
#' - data.table
#' - GenomicRanges
#' - log4r
#' - qs
#' - readxl
#' - rtracklayer
#' 

## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting preprocessCNV.R\n")
    
# 0.2 Read in the sample metadata
# -------------------------------
info("Reading in sample metadata")
sample <- readxl::read_excel(INPUT$sampleMetadata, sheet = 1, col_names = TRUE, na = "NA")
sample <- data.table::as.data.table(sample)
info(paste0("Number of samples: ", nrow(sample)))

# 0.3 Read in the treatment metadata
# ----------------------------------
info("Reading in treatment metadata")
treatment <- data.table::fread(INPUT$treatmentMetadata, header = TRUE)
info(paste0("Number of treatments: ", nrow(treatment)))

# 0.4 Read in the Cell Model Passport metadata
# --------------------------------------------
info("Reading in Cell Model Passport metadata")
cellmodels <- data.table::fread(INPUT$cellModelPassportsMetadata, header = TRUE)

info("Subsetting Cell Model Passport metadata to only samples in GDSC sample metadata")
gdsc_cmp <- cellmodels[model_name %in% sample$`Sample Name`]
cols <- c(
    'model_id', 'sample_id', 'model_name', 'tissue', 'cancer_type', 
    'cancer_type_ncit_id', 'sample_site', 'COSMIC_ID', 'BROAD_ID', 'CCLE_ID',
    'gender', 'ethnicity', 'age_at_sampling', 'clinical_staging'
)

gdsc_cmp <- gdsc_cmp[, ..cols]

sampleMetadata <- merge(sample[,.(`Sample Name`)], gdsc_cmp, by.x = "Sample Name", by.y = "model_name", all.x = TRUE)
# rename Sample Name to sample_id
data.table::setnames(sampleMetadata, "Sample Name", "sampleid")

info(capture.output(str(sampleMetadata)))

# 3. Build Genomic Annotation
# ---------------------------
geneAnnot <- data.table::fread(INPUT$cellmodelpassportsGeneAnnotation, header = TRUE)

gencodeAnnot <- rtracklayer::import(INPUT$gencode_annotation_file)

dsGencode_dt <- data.table::as.data.table(gencodeAnnot)

# remove version numbers from gene_id, transcript_id, havana_transcript, exon_id
# i.e ENST00000456328.2  -> ENST00000456328
dsGencode_dt[, 
    c("gene_id", "transcript_id", "havana_transcript", "exon_id") 
        := lapply(.SD, function(x) gsub("\\.\\d+$", "", x)), 
    .SDcols = c("gene_id", "transcript_id", "havana_transcript", "exon_id")]

# remove columns that are not needed
cols <- c(
    "seqnames", "start", "end", "strand", 
    "gene_id", "gene_name", 
    "gene_type", "source", "tag")

genes_dt <- merge(
    geneAnnot,
    dsGencode_dt[type == "gene", ..cols],
    by.x = "ensembl_gene_id",
    by.y = "gene_id",
    all.x=TRUE,      #
    sort=FALSE
)
data.table::setnames(genes_dt, c("gene_name", "gene_id"), c("symbol", "CMP_gene_id"))

GRanges <- GenomicRanges::makeGRangesFromDataFrame(
    df = genes_dt, keep.extra.columns=TRUE, na.rm = T,
    start.field = "start", end.field = "end", seqnames.field = "seqnames")

GRanges <- GRanges[!duplicated(GRanges$symbol), ]

outputFiles <- list(
    "sample" = sampleMetadata,
    "treatment" = treatment,
    "geneAnnot" = genes_dt,
    "GRanges" = GRanges
)

dir.create(dirname(OUTPUT$metadata), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$metadata)
