
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    # save.image()
}

sample <- readxl::read_excel(INPUT$sampleMetadata, sheet = 1, col_names = TRUE, na = "NA")
sample <- data.table::as.data.table(sample)

treatment <- data.table::fread(INPUT$treatmentMetadata, header = TRUE)

cellmodels <- data.table::fread(INPUT$cellModelPassportsMetadata, header = TRUE)

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

# 3. Build Genomic Annotation
# ---------------------------
geneAnnot <- data.table::fread(INPUT$cellmodelpassportsGeneAnnotation, header = TRUE)
path <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/GDSC/metadata/human/GRCh38_v44/annotation.gtf"
dsGencode <- rtracklayer::import(path)
dsGencode_dt <- data.table::as.data.table(dsGencode)

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


outputFiles <- list(
    "sample" = sampleMetadata,
    "treatment" = treatment,
    "geneAnnot" = genes_dt,
    "GRanges" = GRanges
)

dir.create(dirname(OUTPUT$metadata), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputFiles, file = OUTPUT$metadata)
