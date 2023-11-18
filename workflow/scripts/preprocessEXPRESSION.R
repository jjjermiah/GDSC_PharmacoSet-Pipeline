# FILEPATH: /home/bioinf/bhklab/jermiah/psets/pset-pipelines/GDSC/scripts/processExpressionData.R
# PACKAGE DEPENDENCIES:
#  - affy
#  - readr

# Set the directory paths for the raw data and metadata
rawdata_Dir <- "/home/bioinf/bhklab/jermiah/psets/pset-pipelines/GDSC/rawdata"
metadata_Dir <- "/home/bioinf/bhklab/jermiah/psets/pset-pipelines/GDSC/metadata"

## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
       
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}

celDirectory <- INPUT$CELfiles
sample_annotations <- INPUT$CEL_metadata

# ----------------------------- EXPRESSION DATA ----------------------------- ##
# Get a list of all the .cel files in the rawdata expression directory
fileList <- list.files(celDirectory, pattern = ".cel", full.names = TRUE)

eSet <- affy::rma(affy::read.affybatch(fileList))   # Read in the .cel files using the affy package and perform RMA normalization
mtx <- affy::exprs(eSet)                            # Extract the expression matrix from the normalized data

colnames(mtx) <- sub(pattern = ".cel", "",          # Remove the ".cel" extension from the column names of the expression matrix
                    colnames(mtx), fixed = TRUE) 


## ----------------------------- SAMPLE ANNOTATIONS ----------------------------- ##
# NOTE: These sample annotations are specific to the affymetrix data
# The rawdata files use Assay IDs specific to this technology and have provided annotations 
# to map these IDs back to the sample name 
sample_annotations <- data.table::fread(sample_annotations, header=TRUE, sep="\t")

# Get the common column names between the matrix and sample annotations
common_colnames <- intersect(colnames(mtx), sample_annotations$`Assay Name`)

# Subset the matrix and sample annotations to only include the common columns
mtx <- mtx[, colnames(mtx) %in% common_colnames]
sample_annotations <- sample_annotations[sample_annotations$`Assay Name` %in% common_colnames, ]
sample_annotations <- sample_annotations[match(colnames(mtx), sample_annotations$`Assay Name`), ]

# Check if the column names of the matrix and sample annotations match
all.equal(colnames(mtx), sample_annotations$`Assay Name`)

## ----------------------------- GENE ANNOTATIONS ----------------------------- ##
library(hgu219.db)
k <- keys(hgu219.db,keytype="PROBEID")

columns <- c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "UNIPROT", "PMID") 

annotations <- BiocParallel::bplapply(columns, function(col){
        data <- data.table::as.data.table(select(hgu219.db, keys=k, columns=c(col), keytype="PROBEID"))
        # remove rows with dupllicate "PROBEID" values
        data <- data[!duplicated(data$PROBEID), ]
        data
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = THREADS)
)

# merge all data.tables on "PROBEID"
gene_annotations <- Reduce(function(x, y) merge(x, y, by = "PROBEID", all = TRUE), annotations)

# # Get common gene names
common_genes <- intersect(rownames(mtx), gene_annotations$PROBEID)

# # Subset and match both matrices
mtx <- mtx[rownames(mtx) %in% common_genes, ]
gene_annotations <- gene_annotations[gene_annotations$PROBEID %in% common_genes, ]
gene_annotations <- gene_annotations[match(rownames(mtx), gene_annotations$PROBEID), ]
all.equal(rownames(mtx), gene_annotations$PROBEID)  # Check if matching worked


mtx_collapsed <- WGCNA::collapseRows(datET = mtx, rowGroup = gene_annotations$SYMBOL, rowID = rownames(mtx))$datETcollapsed
colnames(mtx_collapsed) <- sample_annotations$`Characteristics[cell line]`

expr <- data.table::as.data.table(mtx_collapsed, keep.rownames = "gene.symbol")

# subset gene_annotations on SYMBOL with the values in expr$gene.symbol
# only keep the rows in gene_annotations that have a match in expr

expr_annot <- gene_annotations[gene_annotations$SYMBOL %in% expr$gene.symbol, .(SYMBOL, GENENAME, ENSEMBL, ENTREZID, UNIPROT, PMID)]


# ----------------------------- OUTPUT ----------------------------- ##


# combined expr_annot and sample_annotations and mtx into a list of 3 items
outputList <- list(expr_annot, sample_annotations, expr)
names(outputList) <- c("expr_annot", "sample_annotations", "expr")


# make output directory if it doesnt exist
dir.create(dirname(OUTPUT[['preprocessedEXPRESSION']]), recursive = TRUE, showWarnings = FALSE)
qs::qsave(outputList, OUTPUT[['preprocessedEXPRESSION']], nthreads = THREADS)
