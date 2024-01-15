
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(CoreGx))
suppressPackageStartupMessages(library(GenomicRanges))


# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting build_treatmentResponseExperiment.R\n")

# 0.2 read in treatmentResponse data
# ----------------------------------
info(paste0("Reading in: ", INPUT$rawdata))
rawdata <- data.table::fread(INPUT$rawdata)

info(paste0("Reading in: ", INPUT$processed))
data <- readxl::read_xlsx(INPUT$processed)
data <- data.table::as.data.table(data)

# 0.3 read in metadata 
# --------------------
info(paste0("Reading in: ", INPUT$metadata))
metadata <- qs::qread(INPUT$metadata)
treatment <- metadata$treatment

# 1.0 Subsetting Data
# -------------------

rawdata[, TAG.1 := sapply(strsplit(TAG, split = "-"), function(x) {
  return(x[[1]])
})]

rawdata

## Taking mean if more than 5 measurements, otherwise median

rawdata[, Viability := {
  control <- .SD[TAG == ..control.column, ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
  background <- .SD[TAG.1 == "B", ifelse(length(INTENSITY) > 5, mean(INTENSITY), median(INTENSITY))]
  Viability <- (INTENSITY - background) / (control - background) * 100
  list(Viability = Viability)
}, .(MASTER_CELL_ID, BARCODE)]

# NEED TO DO TRE DATA MAPPING! 