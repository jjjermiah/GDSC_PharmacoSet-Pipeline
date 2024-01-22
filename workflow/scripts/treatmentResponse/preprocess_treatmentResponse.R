
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    # save.image()
}

suppressPackageStartupMessages(library(data.table, quietly = TRUE))

# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(log4r::file_appender(LOGFILE, append = TRUE)))

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting preprocess_treatmentResponse.R\n")

# 0.2 read in treatmentResponse data
# ----------------------------------
info(paste0("Reading in: ", INPUT$rawdata))
rawdata <- data.table::fread(INPUT$rawdata)
info(paste0(capture.output(str(rawdata)), collapse = "\n"))
info(paste0("Number of unique cell lines: ", uniqueN(rawdata$CELL_LINE_NAME)))
info(paste0("Number of unique drugs: ", uniqueN(rawdata$DRUG_ID), "\n\n"))

info(paste0("Reading in: ", INPUT$processed))
data <- readxl::read_xlsx(INPUT$processed)
data <- data.table::as.data.table(data)
info(paste0(capture.output(str(data)), collapse = "\n"))
info(paste0("Number of unique cell lines: ", uniqueN(data$CELL_LINE_NAME)))
info(paste0("Number of unique drugs: ", uniqueN(data$DRUG_ID), "\n\n"))

# 0.3 read in metadata 
# --------------------
info(paste0("Reading in: ", INPUT$metadata))
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample
treatment <- metadata$treatment

# 1.0 Subsetting Data
# -------------------
dim_before_subset_rawdata <- dim(rawdata)
dim_before_subset_data <- dim(data)
info(paste0("Subsetting rawdata & data to only include samples in metadata"))
rawdata_s <- rawdata[CELL_LINE_NAME %in% sample$sampleid]

info(paste0("Subsetting rawdata to only include treatments in metadata"))
rawdata_s_t <- rawdata_s[
    (!is.na(DRUG_ID) & (DRUG_ID %in% treatment$DRUG_ID)) | 
        grepl("^(L|R|A|N|B)", TAG)]

if(!all(unique(data$DRUG_ID) %in% treatment$DRUG_ID)){
    info(paste0("ERROR: Not all drugs in data are in metadata"))
    # data <- data[DRUG_ID %in% treatment$DRUG_ID]
}

subset_rawdata <- 
    merge(rawdata_s_t, treatment, by.x = "DRUG_ID", by.y = "DRUG_ID", all.x = T) 

if(!all(unique(subset_rawdata[!is.na(DRUG_NAME),DRUG_NAME]) %in% treatment$DRUG_NAME)){
    info(paste0("ERROR: Not all drugs in subset_rawdata are in metadata"))
    # subset_rawdata <- subset_rawdata[DRUG_NAME %in% treatment$DRUG_NAME]
}



dim_after_subset_rawdata <- dim(subset_rawdata)

info(paste0("Number of rows removed rawdata: ", dim_before_subset_rawdata[1] - dim_after_subset_rawdata[1]))

###############################################################################
# THESE STEPS ARE GDSC SPECIFIC. 
# -----------------------------------------------------------------------------
# ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/GDSC_Raw_Data_Description.pdf
# for more information

# need to install gdscIC50 package from github


while (file.exists(".lock")) {
    info("Directory is locked. Waiting 60 seconds and trying again...")
    Sys.sleep(60)
}

if (!require("gdscIC50", quietly = TRUE))
        devtools::install_github("cancerrxgene/gdscIC50", build_vignettes=FALSE)

library(gdscIC50)

# drop all rows where TAG starts with L and DRUG_ID is NA 
dt <- subset_rawdata[!(grepl("^L", TAG) & is.na(DRUG_ID)) | !TAG == 'FAIL']


dt <- gdscIC50::removeFailedDrugs(subset_rawdata)

dt <- gdscIC50::removeMissingDrugs(subset_rawdata)

# Need to identify the negative control tag for normalizeData
# Should be either NC-0 or NC-1 

neg_control_TAGS <- dt[grepl("^NC", TAG), unique(TAG)]
if(length(neg_control_TAGS) > 1){
    info(paste0("More than one negative control found: ", paste0(neg_control_TAGS, collapse = ", ")))
    info(paste0("Setting all negative controls to: ", neg_control_TAGS[1]))
    dt[TAG %in% neg_control_TAGS, TAG := neg_control_TAGS[1]]
    neg_control_tag <- neg_control_TAGS[1]
}else if(length(neg_control_TAGS) == 0){
    stop("No negative control found")
}else{
    info(paste0("Negative control found: ", neg_control_TAGS))
    neg_control_tag <- neg_control_TAGS[1]
}



# NOTE: the GDSC documentation describes "B" as BLANK (no drug, no cells, just media)
# but also PC-1 as positive control, No titration of this positive control in the 
# drug set.

info(paste0("Normalizing data with negative control: ", neg_control_tag))
normData <- suppressWarnings(
    gdscIC50::normalizeData(dt, trim=T, neg_control=neg_control_tag, pos_control="B"))


# normalizeData drops the DRUG_NAME column
info(paste0("Merging DRUG_NAME back into normData"))
normData <- merge(normData, treatment[,.(DRUG_ID, DRUG_NAME)], by.x = "DRUG_ID_lib",  by.y = "DRUG_ID", all.x = T)

# # FOR NOW, REMOVE ALL ROWS WHERE DRUG_NAME HAS A PUNCTUATION OR SPACE
normData <- unique(normData[!grepl("[[:punct:]]", DRUG_NAME)])
# # TODO:: FIXME:: WHEN MAPPING TREATMENTS, CLEAN NAMES!!

info(paste0("Saving normData to: ", OUTPUT$preprocessed))
outputfiles <- list(
    rawdata = rawdata,
    fitted_data = data,
    normData = normData,
    metadata = metadata
)

if(!dir.exists(dirname(OUTPUT$preprocessed))){
    dir.create(dirname(OUTPUT$preprocessed))
}

qs::qsave(outputfiles, OUTPUT$preprocessed, nthreads = THREADS)
