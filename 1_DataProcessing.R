#### Load libraries ####

library(SummarizedExperiment)
library(TCGAbiolinks)
require(dplyr)
require(NOISeq)
library(DESeq2)
library(biomaRt)
library(EDASeq)

#### Prepare functions ####

filter_TCGA <- function(x) {dataFilt <- TCGAanalyze_Filtering(tabDF = x,
                                                              method = "quantile",
                                                              qnt.cut = 0.25)
threshold <- round(dim(x)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
x <- x[rownames(x) %in% rownames(dataFilt), ]
print(dim(x))
return(x)
}

get_annot <- function(x,y) {
  inter <- intersect(rownames(x), y$Ensembl_ID_Version)
  length(inter)
  annot1 <- y[y$Ensembl_ID_Version  %in% inter,]
  print(dim(annot1))
  annot1 <- annot1[!duplicated(annot1$Ensembl_ID_Version),]
  print(dim(annot1))
  x <- x[rownames(x) %in% annot1$Ensembl_ID_Version,]
  print(dim(x))
  x <- x[!duplicated(rownames(x)),] #
  print(dim(x))
  print(head(rownames(x)))
  print(head(annot1$Ensembl_ID_Version))
  annot1 <- annot1[match(rownames(x), annot1$Ensembl_ID_Version), ]
  print(dim(annot1))
  print(dim(x))
  print(head(rownames(x)))
  print(head(annot1$Ensembl_ID_Version))
  return(list(x,annot1))
}

norm <- function(x, y, z) {
  ln.data <- withinLaneNormalization(x, y$Length, which = "full")
  gcn.data <- withinLaneNormalization(ln.data , y$GC, which = "full")
  norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
  noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(z$Group))
  mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
  rnas2 <- exprs(mydata2corr1)
  return(rnas2)
}

#### Get annotation file ####
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")

features <- c("ensembl_gene_id", "chromosome_name", 
              "start_position", "end_position", "hgnc_symbol",	
              "percentage_gene_gc_content", "gene_biotype", "ensembl_gene_id_version", "hgnc_id")
chrs <- c(1:22, "X", "Y")

annot <- getBM(attributes = features,
               filters = "chromosome_name",
               values = chrs, 
               mart = ensembl)

colnames(annot)<-c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type", "Ensembl_ID_Version", "HGNC_ID")
annot$Length <- abs(annot$End - annot$Start)

#### Get TCGA data ####

#DLBCL 

qry.rna_DLBCL <- GDCquery(project = "NCICCR-DLBCL",
                          data.category= "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")
GDCdownload(qry.rna_DLBCL)
DLBCL <- GDCprepare(qry.rna_DLBCL, summarizedExperiment = TRUE) 

dir.create("Results")
saveRDS(as.data.frame(rowData(DLBCL)), "Results/row_data.RDS")

#### Match clinical information with COO classification ####
#The classification of samples based on the cell of origin framework was manually obtained from the supplementary material
#of original publication of the data ("Genetics and Pathogenesis of Diffuse Large B-Cell Lymphoma", DOI: 10.1056/NEJMoa1801445)

cases_by_subtype <- read.csv("Inputs/cases_by_subtype.csv")
cases_by_subtype <- cases_by_subtype[cases_by_subtype$dbGaP.subject.ID %in% DLBCL$submitter_id,]

CI <- as.data.frame(colData(DLBCL))
cases_by_subtype <- cases_by_subtype %>% dplyr::select(dbGaP.subject.ID, Gene.Expression.Subgroup)
CI <- CI %>% inner_join(cases_by_subtype, by = c("submitter_id" = "dbGaP.subject.ID"))

saveRDS(CI, "Results/clinical_info.RDS")

factors_Lymph <- data.frame(Group = CI$Gene.Expression.Subgroup, Sample = CI$sample_submitter_id)
Ready_factors_Lymph <- as.data.frame(factors_Lymph$Group)

#Get counts assay
rnas_Lymph <- assay(DLBCL)
rownames(rnas_Lymph) <- rowData(rnas_Lymph)$gene_id
rnas_Lymph <- rnas_Lymph[!duplicated(rownames(rnas_Lymph)),]

#### Filtering ####

rnas_Lymph <- filter_TCGA(rnas_Lymph)

##### Get annot files ####

Lymph <- get_annot(rnas_Lymph, annot)
Lymph_annot <- Lymph[[2]]
saveRDS(Lymph_annot, "Results/Lymph_annot.RDS")

##### DEG analysis using DESeq2 ####

lymph_raw <- Lymph[[1]]
dir.create("Results/DEG")
#GCB vs ABC
dds <- DESeqDataSetFromMatrix(countData = round(lymph_raw[,factors_Lymph$Group=="GCB" | 
                                                            factors_Lymph$Group=="ABC"]),
                              colData = factors_Lymph[factors_Lymph$Group=="GCB" | 
                                                        factors_Lymph$Group=="ABC",],
                              design = ~ Group)


dds <- DESeq(dds)
res <-  results(dds)
res
summary(res)
resLFC <- lfcShrink(dds, coef="Group_GCB_vs_ABC", type="apeglm")
write.table(resLFC, file = "Results/DEG/resLFC_Group_GCB_vs_ABC.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#Unclass vs GCB
dds <- DESeqDataSetFromMatrix(countData = round(lymph_raw[,factors_Lymph$Group=="GCB" | 
                                                            factors_Lymph$Group=="Unclass"]),
                              colData = factors_Lymph[factors_Lymph$Group=="GCB" | 
                                                        factors_Lymph$Group=="Unclass",],
                              design = ~ Group)


dds <- DESeq(dds)
res <-  results(dds)
res
summary(res)
resLFC <- lfcShrink(dds, coef="Group_Unclass_vs_GCB", type="apeglm")
write.table(resLFC, file = "Results/DEG/resLFC_Group_Unclass_vs_GCB.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#Unclass vs ABC
dds <- DESeqDataSetFromMatrix(countData = round(lymph_raw[,factors_Lymph$Group=="ABC" | 
                                                            factors_Lymph$Group=="Unclass"]),
                              colData = factors_Lymph[factors_Lymph$Group=="ABC" | 
                                                        factors_Lymph$Group=="Unclass",],
                              design = ~ Group)


dds <- DESeq(dds)
res <-  results(dds)
res
summary(res)
resLFC <- lfcShrink(dds, coef="Group_Unclass_vs_ABC", type="apeglm")
write.table(resLFC, file = "Results/DEG/resLFC_Group_Unclass_vs_ABC.tsv", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)

#### Get norm data ####

Lymph_norm <- norm(Lymph[[1]], Lymph[[2]], Ready_factors_Lymph)
saveRDS(Lymph_norm, "Results/Norm_data.RDS")

#### Get norm count matrices ####
dir.create("Results/NormalizedData")
ABC_DLBCL_Norm <- Lymph_norm[, factors_Lymph$Group=="ABC"]
ABC_DLBCL_Norm <- ABC_DLBCL_Norm %>% dplyr::mutate(Gene = rownames(ABC_DLBCL_Norm)) %>% dplyr::relocate(Gene)
write.table(ABC_DLBCL_Norm, file = "Results/NormalizedData/rnas_norm_ABC_DLBC.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

GCB_DLBCL_Norm <- Lymph_norm[, factors_Lymph$Group=="GCB"]
GCB_DLBCL_Norm <- GCB_DLBCL_Norm %>% dplyr::mutate(Gene = rownames(ABC_DLBCL_Norm)) %>% dplyr::relocate(Gene)
write.table(GCB_DLBCL_Norm, file = "Results/NormalizedData/rnas_norm_GCB_DLBCL.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)

Unclass_DLBCL_Norm <- Lymph_norm[, factors_Lymph$Group=="Unclass"]
Unclass_DLBCL_Norm <- Unclass_DLBCL_Norm %>% dplyr::mutate(Gene = rownames(ABC_DLBCL_Norm)) %>% dplyr::relocate(Gene)
write.table(Unclass_DLBCL_Norm, file = "Results/NormalizedData/rnas_norm_Unclass_DLBCL.tsv", row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)









