library(clusterProfiler)
library(GOSemSim)
library(DOSE)
library(dplyr)

communities <- read.csv("communities_ABC_4.csv", row.names = 1)
communities_t <- t(communities)
communities_t <- gsub("\\.[0-9]*", "", communities_t)


for(i in 1:ncol(communities_t)) {

  FuncEnrich <- enrichGO(gene = communities_t[,i], OrgDb = "org.Hs.eg.db", ont = "BP",
                         keyType = "ENSEMBL", pvalueCutoff=0.00001)

  #We applied the redundancy reduction function simplify to remove similar biological processes

  if(!is.null(FuncEnrich)) {

    if(!is.na(FuncEnrich[1]$ID)) {
      simpl_FuncEnrich <- simplify(FuncEnrich, cutoff=0.7, by="p.adjust", select_fun=min) %>% as.data.frame()


      network <- data.frame(Community = rep(colnames(communities_t)[i], each=nrow(simpl_FuncEnrich)),
                            Interaction = rep(1, each=nrow(simpl_FuncEnrich)),
                            GO_BP = simpl_FuncEnrich$Description,
                            ID = simpl_FuncEnrich$ID,
                            GeneRatio = simpl_FuncEnrich$GeneRatio,
                            BgRatio = simpl_FuncEnrich$BgRatio,
                            pvalue = simpl_FuncEnrich$pvalue,
                            p.adjust = simpl_FuncEnrich$p.adjust,
                            qvalue = simpl_FuncEnrich$qvalue,
                            geneID = simpl_FuncEnrich$geneID,
                            Count = simpl_FuncEnrich$Count,
                            row.names=NULL)

      if(!exists("x")){
        x <- network
      } else {
        x <- rbind(x, network)

      }

    }

  }}

write.table(x, file = "ABC_Communities_Enriched_4.sif", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
