## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    LOGFILE <- snakemake@log[[1]]
    save.image()
}

se_list <- unlist(lapply(INPUT$summarizedExperiments, qs::qread, nthreads = THREADS))


sampleid <- unique(unlist(lapply(se_list, colnames)))
sample <- data.frame(sampleid, row.names = sampleid)

se_list <- lapply(se_list, function(se){
    se[, colnames(se) %in% sampleid]
})

colData <- data.frame(
    sampleid = sampleid,
    batchid = rep(NA, length(sampleid)),
    row.names = sampleid
)

ExpList <- MultiAssayExperiment::ExperimentList(se_list)

sampleMapList <- lapply(se_list, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})
names(sampleMapList) <- names(ExpList)
MultiAssayExperiment::listToMap(sampleMapList)

mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList)
)
mae
