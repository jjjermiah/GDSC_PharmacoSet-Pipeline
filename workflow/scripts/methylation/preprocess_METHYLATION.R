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
info("Starting preprocess_METHYLATION.R\n")


# 0.2 Load input data
# -------------------
# meth_cell_data is a zip file 
zipDir <- dirname(INPUT$meth_cell_data)
unzipDir <- file.path(zipDir, "meth_cell_data")
info(paste0("Unzipping ", INPUT$meth_cell_data, " to ", unzipDir))
unzip(INPUT$meth_cell_data, exdir = unzipDir)

# should be one file and a __MACOSX folder (ignore)
files <- list.files(unzipDir)
info(sprintf("Files in %s :\n%s", unzipDir, paste(files, collapse = "\n")))

dataFile <- file.path(unzipDir, list.files(unzipDir, pattern = "*.txt"))
dt <- data.table::fread(dataFile, sep = "\t", header = TRUE) 

# summarize dt for logs 
info(paste0("dt has ", nrow(dt), " rows and ", ncol(dt), " columns"))
info(paste(c("Summary: \n", capture.output(dt[1:5,1:5])), collapse = "\n"))

# 0.3 Load metadata and sample mapping file
# -----------------------------------------
metadata <- qs::qread(INPUT$metadata)
sample <- metadata$sample
rowRanges <- metadata$GRanges
info(paste(capture.output(str(metadata, 2)),  collapse = "\n"))

map_dt <- data.table::as.data.table(readxl::read_excel(INPUT$sample_mapping_metadata))

#' the column names of dt are some combination of two IDs:
#' "8221924127_R04C02" "8221924165_R06C02" "8359018054_R03C01"
#' it looks like a combination of the "Sentrix_ID" and "Sentrix_Position" from
#' the map_dt
map_dt[, MethylationDataID := paste0(Sentrix_ID, "_", Sentrix_Position)]
stopifnot(all(colnames(dt)[-1] %in% map_dt$MethylationDataID))
col_match <- match(map_dt$MethylationDataID, colnames(dt)[-1])
# replace the column names with the "Sample_Name" from the map_dt
colnames(dt)[-1] <- map_dt$Sample_Name[col_match]

# remove duplicate columns
dt <- dt[, !duplicated(colnames(dt)), with = FALSE]



# # Now that the column names are the sample names, we subset for only columns
# # that are in the sampleid column of sample
common_cols <- intersect(colnames(dt)[-1], sample$sampleid)

# # subset dt to only the common columns
dt <- dt[, c("V1", common_cols), with = FALSE]



# # Convert the "V1" Column into 3 columns: "Chromosome", "Start", "End"
# # -------------------------------------------------------------------
# # the "V1" column is a string of the form "chr1:100000-200000"
# # we want to split this into 3 columns: "Chromosome", "Start", "End"

# dt[, c("Chromosome", "Start", "End") := data.table::tstrsplit(V1, "[:-]")]
# dt[, .(Chromosome, Start, End, V1, `COR-L303`)]

# meth_granges <- GRanges(
#     seqnames = dt$Chromosome, 
#     ranges = IRanges(as.numeric(dt$Start), as.numeric(dt$End)))


# meth_granges
# library(GenomicFeatures)
# path <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/GDSC/metadata/human/GRCh38_v44/annotation.gtf"

# txdb <- GenomicFeatures::makeTxDbFromGFF(path, format = "gtf")

# annotated_granges <- annotatePeakInBatch(meth_granges, txdb, output="nearest")

# meth_granges <- meth_granges[order(meth_granges@seqnames, meth_granges@ranges)]
# meth_granges
# # find overlaps meth_granges and rowRanges
# # ------------------------------------
# # meth_granges is the GRanges object created from the methylation data
# # rowRanges is the GRanges object created from the cell model passport data 

# overlaps_ <- findOverlaps(meth_granges, rowRanges)
# overlaps_

# # create GRanges object from overlaps
# # -----------------------------------
# # the overlaps_ object is a S4 object of class "Hits"
# # we can use the hits to subset the meth_granges object and get the additional 
# # metadata from the rowRanges object

# new_granges <- rowRanges[queryHits(overlaps_)]
# # order by range column
# new_granges <- new_granges[order(new_granges@seqnames, new_granges@ranges)]

# new_granges

# # add the additional metadata from rowRanges to the methylation data
# # ------------------------------------------------------------------
# # get all possible columns from rowRanges 

# rowRanges_cols <- names(rowRanges)





# # subset dt to only the common columns
# dt[, c("V1", common_cols), with = FALSE]


# annotation <- read.csv("GPL13534_HumanMethylation450_15017482_v.1.1.csv", skip = 7)
# annotation <- data.table::as.data.table(annotation)

# probe_mapping <- annotation[, c("IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name")]
# merge(dt, probe_mapping, by.x = "V1", by.y = "IlmnID", all.x = TRUE)
library(data.table)
gpl <- GEOquery::getGEO('GPL13534')
gpl <- data.table::as.data.table(gpl@dataTable@table)

# remove rows where UCSC_RefGene_Name is empty
gpl <- gpl[!UCSC_RefGene_Name == "",]

# for each row of gpl, create an additional row if there are multiple values in
# the UCSC_RefGene_Name column
gpl[, symbol := tstrsplit(UCSC_RefGene_Name, ";"), by = 1:nrow(gpl)]
gpl <- unique(gpl[symbol %in% rowRanges$symbol,])
gpl <- gpl[UCSC_CpG_Islands_Name %in% dt$V1,]

gpl_s <- gpl[, .(symbol, Name, UCSC_CpG_Islands_Name)]

dt_s <- dt[V1 %in% gpl_s$UCSC_CpG_Islands_Name,]

merged_dt <- merge(dt_s, gpl_s, by.x = "V1", by.y = "UCSC_CpG_Islands_Name", all.x = TRUE)

merged_dt <- merged_dt[!duplicated(Name),]
merged_dt <- merged_dt[!duplicated(symbol), !c("Name", "V1")]

# order both rows and columns alphabetically by symbol
merged_dt <- merged_dt[order(symbol), order(colnames(merged_dt)), with = FALSE]

# convert to dataframe with rownames as symbol and colnames as sampleid and drop symbol column

merged_df <- as.data.frame(merged_dt[, !c("symbol")], row.names = merged_dt$symbol)
rownames(merged_df) <- merged_dt$symbol
merged_df[1:10, 1:10]

mtx <- as.matrix(merged_df)
mtx[1:10, 1:10]

rownames(mtx) %in% rowRanges$symbol
colnames(mtx) %in% sample$sampleid

assay_mtx <- unique(mtx)

gr <- rowRanges[rowRanges$symbol %in% rownames(mtx),]
sample_ <- sample[sample$sampleid %in% colnames(mtx),][!duplicated(sampleid),]
rse <- SummarizedExperiment::SummarizedExperiment(
    assays = list("meth" = mtx),
    rowRanges = gr,
    colData = sample_
)
rse

# "chr1:10003165-10003585" %in% annot$UCSC_CpG_Islands_Name
# annot[annot$UCSC_CpG_Islands_Name == "chr1:10003165-10003585", ]



# # remove extra first column for row number
# mtx_dt <- data.table::fread("GSE68379_Matrix.processed.txt")[,-1]

# mtx_dt[,1:10]


# # only keep "Row.names" column and any column that ends in ".Beta"
# mtx_dt_subset <- mtx_dt[, c("Row.names", grep(".Beta", colnames(mtx_dt), value = TRUE)), with = FALSE]

# # remove the .Beta from the column names
# colnames(mtx_dt_subset)[-1] <- gsub("_AVG.Beta", "", colnames(mtx_dt_subset)[-1])

# # only keep columns that are in sample$sampleid 
# common_cols <- intersect(colnames(mtx_dt_subset)[-1], sample$sampleid)
# mtx_dt_subset <- mtx_dt_subset[, c("Row.names", common_cols), with = FALSE]

# stopifnot(all(mtx_dt_subset$Row.names %in% annot$Name))

# probe.annotation <- rnb.get.annotation("probes450")

# # combine probe.annotation into one GRanges object
# probe.annotation <- do.call(c, probe.annotation)
# probe.annotation # GRangesList object
# probe.annotation$Name <- rownames(probe.annotation)
# rbindlist(lapply(probe.annotation, data.table::as.data.table, keep.rownames=T), use.names = TRUE)

# gene.annotation <- rnb.get.annotation("genes")
# attr(gene.annotation, "version")
# gene.annotation
