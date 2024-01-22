
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}
pak::pkg_install("bhklab/CoreGx")

suppressPackageStartupMessages(library(data.table, quietly = TRUE))
suppressPackageStartupMessages(library(PharmacoGx))
suppressPackageStartupMessages(library(GenomicRanges))

# 0.1 Setup Logger
# ----------------
# create a logger from the LOGFILE path in append mode
logger <- log4r::logger(
    appenders = list(
      log4r::file_appender(LOGFILE, append = TRUE),
      log4r::console_appender()
    )
)

# make a function to easily log messages to the logger
info <- function(msg) log4r::info(logger, msg)
info("Starting fit_treatmentResponseExperiment.R\n")


# 0.2 Load Data
# -------------
info(paste0("Reading in: ", INPUT$preprocessed))
preprocessed <- qs::qread(INPUT$preprocessed, nthreads = THREADS)
info(paste0(capture.output(str(preprocessed)), collapse = "\n"))

info(paste0("Reading in: ", INPUT$tre))
tre <- qs::qread(INPUT$tre, nthreads = THREADS)


info(paste0(capture.output(tre), collapse = "\n"))

drugs <- rowData(tre)[!grepl("[[:punct:]]|[[:space:]]",DRUG_NAME), unique(DRUG_NAME)]
samples <- colData(tre)[, unique(CELL_LINE_NAME)]

# # FOR NOW, REMOVE ALL ROWS WHERE DRUG_NAME HAS A PUNCTUATION OR SPACE
# subset_procdata <- unique(subset_procdata[!grepl("[[:punct:]]|[[:space:]]", DRUG_NAME),])
# # TODO:: FIXME:: WHEN MAPPING TREATMENTS, CLEAN NAMES!!

tre_s <- tre[.(DRUG_NAME %in% drugs), .(CELL_LINE_NAME %in% samples),]

info("Endoaggregating tre with recomputed profiles")
tre_fit <- tre_s |>
    CoreGx::endoaggregate(
        {  # the entire code block is evaluated for each group in our group by
            # 1. fit a log logistic curve over the dose range
            fit <- PharmacoGx::logLogisticRegression(CONC, Viability,
                viability_as_pct=FALSE)
            # 2. compute curve summary metrics
            ic50 <- PharmacoGx::computeIC50(CONC, Hill_fit=fit)
            aac <- PharmacoGx::computeAUC(CONC, Hill_fit=fit)
            # 3. assemble the results into a list, each item will become a
            #   column in the target assay.
            list(
                HS=fit[["HS"]],
                E_inf = fit[["E_inf"]],
                EC50 = fit[["EC50"]],
                Rsq=as.numeric(unlist(attributes(fit))),
                aac_recomputed=aac,
                ic50_recomputed=ic50
            )
        },
        assay="raw",
        target="profiles_recomputed",
        enlist=FALSE,  # this option enables the use of a code block for aggregation
        by=c("DRUG_NAME", "CELL_LINE_NAME"),
        nthread=THREADS  # parallelize over multiple cores to speed up the computation
)
info(paste0(capture.output(tre_fit), collapse = "\n"))

# 0.3 Save Data
# -------------
qs::qsave(tre_fit, OUTPUT$tre, nthreads = THREADS)
