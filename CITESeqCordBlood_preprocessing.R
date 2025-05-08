### ------------------------------------------------------------------------ ###
### CITEseq Cord Blood - PREPROCESSING
### ------------------------------------------------------------------------ ###

rm(list=ls())

# ---------------------------------------------------------------------------- #
# PACKAGES
# ---------------------------------------------------------------------------- #

library(cellbaseR)
library(EnsDb.Hsapiens.v79)
library(MultiAssayExperiment)
library(SingleCellMultiModal)
library(SingleCellExperiment)

# ---------------------------------------------------------------------------- #
# DATA
# ---------------------------------------------------------------------------- #

sce <- CITEseq(DataType="cord_blood", modes="*", dry.run=FALSE, version="1.0.0",
               DataClass="SingleCellExperiment")
# class: SingleCellExperiment
# dim: 36280 7858
# metadata(0):
# assays(1): counts
# rownames(36280): ERCC_ERCC-00104 HUMAN_A1BG ... MOUSE_n-R5s25
#   MOUSE_n-R5s31
# rowData names(0):
# colnames(7858): TACAGTGTCTCGGACG GTTTCTACATCATCCC ... TTGCCGTGTAGATTAG
#   GGCGTGTAGTGTACTC
# colData names(6): adt.discard mito.discard ... celltype markers
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(1): scADT

# ---------------------------------------------------------------------------- #
# human genes
# ---------------------------------------------------------------------------- #

# gene names start with "ERCC", "HUMAN" or "MOUSE"
unique(sub("_.*$", "", rownames(sce)))

# keep only the genes starting with "HUMAN"
genes2keep <- sub("_.*$", "", rownames(sce)) == "HUMAN"

# we obtain a smaller dataset with 20400 genes (from 36280)
sce <- sce[genes2keep, ]
# class: SingleCellExperiment
# dim: 20400 7858
# metadata(0):
#   assays(1): counts
# rownames(20400): HUMAN_A1BG HUMAN_A1BG-AS1 ... HUMAN_hsa-mir-7515
# HUMAN_hsa-mir-8072
# rowData names(0):
#   colnames(7858): TACAGTGTCTCGGACG GTTTCTACATCATCCC ... TTGCCGTGTAGATTAG
# GGCGTGTAGTGTACTC
# colData names(6): adt.discard mito.discard ... celltype markers
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(1): scADT

# remove "HUMAN_" from genes name
rownames(sce) <- sub("HUMAN_", "", rownames(sce))

# ---------------------------------------------------------------------------- #
# gene annotation (ensembl IDs)
# ---------------------------------------------------------------------------- #

geneMap <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(sce),
                             keytype = "SYMBOL", columns = c("SYMBOL","GENEID"),
                             multiVals = "first")
# keep the first row for each gene
geneMap <- geneMap[!duplicated(geneMap[,1]), ]

# unfortunately we are not able to map all the genes
nrow(sce)       # 20400
nrow(geneMap)   # 19518
# we just drop those genes
genes2keep <- rownames(sce) %in% geneMap[,1]

# we obtain a smaller dataset with 19518 genes (from 20400)
sce <- sce[genes2keep, ]
# class: SingleCellExperiment
# dim: 19518 7858
# metadata(0):
#   assays(1): counts
# rownames(19518): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
# rowData names(0):
#   colnames(7858): TACAGTGTCTCGGACG GTTTCTACATCATCCC ...
# TTGCCGTGTAGATTAG GGCGTGTAGTGTACTC
# colData names(6): adt.discard mito.discard ... celltype markers
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(1): scADT

# finally we rename the genes
rownames(sce) <- geneMap[geneMap[,1]==rownames(sce), 2]

# ---------------------------------------------------------------------------- #
# pathways selection
# ---------------------------------------------------------------------------- #

pathways <- readRDS("data/CITESeqCordBlood/pathways.RDS")

# 20 immune pathways
immune_pathways = c("Complement and coagulation cascades",
                    "Platelet activation",
                    "Neutrophil extracellular trap formation",
                    "Toll-like receptor signaling pathway",
                    "NOD-like receptor signaling pathway",
                    "RIG-I-like receptor signaling pathway",
                    "Cytosolic DNA-sensing pathway",
                    "C-type lectin receptor signaling pathway",
                    "Natural killer cell mediated cytotoxicity",
                    "Antigen processing and presentation",
                    "T cell receptor signaling pathway",
                    "Th1 and Th2 cell differentiation",
                    "Th17 cell differentiation",
                    "IL-17 signaling pathway",
                    "B cell receptor signaling pathway",
                    "Fc epsilon RI signaling pathway",
                    "Fc gamma R-mediated phagocytosis",
                    "Leukocyte transendothelial migration",
                    "Intestinal immune network for IgA production",
                    "Chemokine signaling pathway")

W_b <- matrix(0, nrow(sce), length(pathways))
rownames(W_b) <- rownames(sce)
colnames(W_b) <- names(pathways)
W_b <- as.data.frame(W_b)
for(i in 1:length(pathways)){
  W_b[ rownames(W_b) %in% pathways[[i]], i] <- 1
  W_b[,i] <- as.factor(W_b[,i])
  col_name <- gsub(" ", "_", names(pathways)[i])
  col_name <- gsub("/", "_", col_name)
  col_name <- gsub("-", "_", col_name)
  colnames(W_b)[i] <- col_name
}
# number of genes in each pathway
genes_num <- apply(W_b, 2, function(x){table(x)[2]})

# let's see how many distinct genes are in the immune system pathways
sum(!is.na(apply(W_b[,names(pathways) %in% immune_pathways], 1, function(x) table(x)[2])))
# only 750 distinct genes in the immune system pathways

# Lets see the distinct genes by keeping always the 20 pathways of the immune system
sum_distinct <- sum(!is.na(apply(W_b[,names(pathways) %in% immune_pathways], 1, function(x) table(x)[2])))
removed_path <- which(names(pathways) %in% immune_pathways)
pathways2keep <- removed_path
genes2rm <- unique(unlist(apply(W_b[,removed_path], 2, function(x) which(x==1))))
remaining <- W_b[-genes2rm, ]
cum_sum <- sum_distinct

# then add the others
while((sum_distinct < 5000) & (length(pathways2keep) < 50)){
  genes_num_tmp <- apply(remaining, 2, function(x){table(x)[2]})
  removed_path <- which.max(genes_num_tmp)
  pathways2keep <- c(pathways2keep, removed_path)
  genes2rm <- which(remaining[,removed_path]==1)
  remaining <- remaining[-genes2rm, ]
  sum_distinct <- sum_distinct + max(genes_num_tmp, na.rm = TRUE)
  cum_sum <- c(cum_sum, sum_distinct)
}
cum_sum
# elbow after 4 further pathways (2613)
plot(cum_sum, type = "l")

# we keep the 20 immune pathways and the remaining 4 largest pathways
W_b <- W_b[, pathways2keep[1:24]]
genes2keep <- !(is.na(apply(W_b, 1, function(x){table(x)[2]})))

sce <- sce[genes2keep, ]
# class: SingleCellExperiment
# dim: 2613 7858
# metadata(0):
#   assays(1): counts
# rownames(2613): ENSG00000128274 ENSG00000081760 ...
#       ENSG00000109576 ENSG00000149313
# rowData names(0):
#   colnames(7858): TACAGTGTCTCGGACG GTTTCTACATCATCCC ...
# TTGCCGTGTAGATTAG GGCGTGTAGTGTACTC
# colData names(6): adt.discard mito.discard ... celltype markers
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(1): scADT

# ---------------------------------------------------------------------------- #
# gene annotation (metacovariates)
# ---------------------------------------------------------------------------- #

cb <- CellBaseR()
# metacovariates of interest: biotype, status, name, chromosome, start, end, strand, source
get_info <- function(mat){
  obj <- data.frame(row.names = mat$id)
  obj$id <- mat$id
  obj$name <- mat$name
  obj$biotype <- mat$biotype
  obj$status <- mat$status
  obj$chromosome <- mat$chromosome
  obj$start <- mat$start
  obj$end <- mat$end
  obj$strand <- mat$strand
  obj$source <- mat$source
  return(obj)
}
# set genes2keep[i] to FALSE if metadata of gene i cannot be found
genes2keep <- rep(TRUE, nrow(sce))
# allocate memory for new gene
newGene <- matrix(NA, nrow = 1, ncol = 9)
# matrix of full info for metacovariates
w_info <- data.frame()
for(i in 1:length(rownames(sce))) {
  cat(i, rownames(sce)[i], "\n")
  newGene <- NULL
  try(newGene <- get_info(getGene(object = cb, ids = rownames(sce)[i], resource = "info")))
  if (is.null(newGene)) {
    genes2keep[i] <- FALSE
    next
  }
  w_info <- dplyr::bind_rows(w_info, newGene)
}
# keep only genes with metadata
 sum(genes2keep)
# we obtain a slightly smaller dataset with 2587 genes (from 2613)
sce <- sce[genes2keep, ]
# class: SingleCellExperiment
# dim: 2587 7858
# metadata(0):
#   assays(1): counts
# rownames(3440): ENSG00000128274 ENSG00000081760 ...
#       ENSG00000109576 ENSG00000149313
# rowData names(9): id name ... strand source
# colnames(7858): TACAGTGTCTCGGACG GTTTCTACATCATCCC ...
# TTGCCGTGTAGATTAG GGCGTGTAGTGTACTC
# colData names(6): adt.discard mito.discard ... celltype markers
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(1): scADT

# ---------------------------------------------------------------------------- #
# technical metacovariates
# ---------------------------------------------------------------------------- #

gc_content <- readRDS("data/CITESeqCordBlood/gccontent.RDS")

W_t <- matrix(0, nrow(sce), 2)
rownames(W_t) <- rownames(sce)

W_t[,1] <- w_info$end- w_info$start
gc_content_clean <- gc_content[names(gc_content) %in% rownames(W_t)]
W_t[names(gc_content_clean), 2] <- gc_content_clean
colnames(W_t) <- c("length", "gc_content")

# ---------------------------------------------------------------------------- #
# save metacovariates
# ---------------------------------------------------------------------------- #

W_b <- W_b[rownames(W_t),]
rowData(sce) <- cbind(W_t, W_b)

# ---------------------------------------------------------------------------- #
# covariates
# ---------------------------------------------------------------------------- #

dim(colData(sce))       # 7858 x 6
colnames(colData(sce))

colData(sce)$adt.discard <- NULL  # FALSE only
colData(sce)$mito.discard <- NULL # FALSE only
colData(sce)$discard <- NULL      # FALSE only
colData(sce)$species <- NULL      # HUMAN only

colData(sce)$celltype <- as.factor(colData(sce)$celltype)
colData(sce)$markers <- as.factor(colData(sce)$markers)

colData(sce)$num.genes <- colSums(as.matrix(counts(sce))>0)
colData(sce)$tot.counts <- colSums(as.matrix(counts(sce)))

dim(colData(sce))       # 7858 x 4
colnames(colData(sce))  # "celltype" "markers" "num.genes", "tot.counts"

# ---------------------------------------------------------------------------- #
# data saving
# ---------------------------------------------------------------------------- #

saveRDS(sce, "data/CITESeqCordBlood/CITEseqCordBlood_an.RDS")

# ---------------------------------------------------------------------------- #
