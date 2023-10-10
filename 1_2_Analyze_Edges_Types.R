
library(dplyr)
library(ggplot2)

######## Load Networks

ABC <- read.delim("ABC_10M.sif")
GCB <- read.delim("GCB_10M.sif")
Unclass <- read.delim("Unclass_10M.sif")

######## Add Chromosomes

Add_Chromosomes_Types <- function(x,y) {
  
  x <- x %>% dplyr::inner_join(y, by = c("source"="Ensembl_ID_Version")) %>% 
    dplyr::select("source", "mi", "target", "Chr", "Type") %>% dplyr::rename("Chr_source"="Chr", "Type_source"="Type") %>%
    dplyr::inner_join(y, by = c("target"="Ensembl_ID_Version")) %>% 
    dplyr::select("source", "mi", "target", "Chr_source", "Type_source", "Chr", "Type") %>% dplyr::rename("Chr_target"="Chr", "Type_target"="Type") %>%
    arrange(desc(mi))
  
  return(x)
}

Lymph_annot <- readRDS("Results/Lymph_annot.RDS")

ABC <- Add_Chromosomes_Types(ABC, Lymph_annot)
dim(ABC)
GCB  <- Add_Chromosomes_Types(GCB, Lymph_annot)
dim(GCB)
Unclass   <- Add_Chromosomes_Types(Unclass, Lymph_annot)
dim(Unclass)

######## Count #########

cis_counting <- function(x) {
  
  h <- 1
  j <- 1
  y <- 1
  cis_count <- 0
  
  while(h*10^(j) <= 1e+07) {
    for(i in y:nrow(x[1:(h*10^(j)),])) {
      if(x$Chr_source[i] == x$Chr_target[i]) { 
        if(h*10^(j) == 10) {
          cis_count = 1
        } else {
          cis_count = 1 + cis_count
        }
      } else {
        cis_count = cis_count
      }
    }
    if(h*10^(j) == 10) {
      cis_matrix <- tibble(total_edges = h*10^(j), cis = (cis_count), cis_proportion = (cis_count/(h*10^(j))))
    } else {
      a <- tibble(total_edges = h*10^(j), cis = (cis_count), cis_proportion = (cis_count/(h*10^(j))))
      cis_matrix <- rbind(cis_matrix, a)
      
    }
    if(h != 9) {
      y = h*10^(j) + 1
      h = h + 1
    } else {
      y = h*10^(j) + 1
      h = 1
      j = j + 1
    }
  }
  return(cis_matrix)
}

ABC_count <- cis_counting(ABC)
GCB_count <- cis_counting(GCB)
Unclass_count  <- cis_counting(Unclass)

global_count <- data.frame(total_edges = ABC_count$total_edges, 
                           ABC_count = ABC_count$cis_proportion,
                           GCB_count = GCB_count$cis_proportion,
                           Unclass_count = Unclass_count$cis_proportion)
write.csv(global_count, file="Results/global_cis_count_lymphomas.csv")

##### Count Gene types linked #####

change_gene_types <- function(x) {
  print(table(as.factor(x$Type_source)))
  x$Type_source <- gsub(".*pseudogene", "pseudogene", x$Type_source)
  x$Type_source <- gsub("IG.*", "IG", x$Type_source)
  x$Type_source <- gsub("TR.*", "TR", x$Type_source)
  print(table(as.factor(x$Type_source)))
  
  print(table(as.factor(x$Type_target)))
  x$Type_target <- gsub(".*pseudogene", "pseudogene", x$Type_target)
  x$Type_target <- gsub("IG.*", "IG", x$Type_target)
  x$Type_target <- gsub("TR.*", "TR", x$Type_target)
  print(table(as.factor(x$Type_target)))
  
  return(x)
  
}

ABC <- change_gene_types(ABC)
dim(ABC)
GCB  <- change_gene_types(GCB)
dim(GCB)
Unclass  <- change_gene_types(Unclass)
dim(Unclass)

link_types <- function(x) {
  for(i in 1:nrow(x)) {
    x$link_type[i] <- paste(pmin(x$Type_source[i],x$Type_target[i]), 
                            pmax(x$Type_source[i],x$Type_target[i]),sep="-")
  }
  return(x)
}

ABC_10k <- link_types(ABC[1:10000,])
dim(ABC_10k)
GCB_10k  <- link_types(GCB[1:10000,])
dim(GCB_10k)
Unclass_10k   <- link_types(Unclass[1:10000,])
dim(Unclass_10k)

ABC_10k_EdgesTypes <- as.data.frame((table(as.factor(ABC_10k$link_type))))
ABC_10k_EdgesTypes <- ABC_10k_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

GCB_10k_EdgesTypes <- as.data.frame((table(as.factor(GCB_10k$link_type))))
GCB_10k_EdgesTypes <- GCB_10k_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

Unclass_10k_EdgesTypes <- as.data.frame((table(as.factor(Unclass_10k$link_type))))
Unclass_10k_EdgesTypes <- Unclass_10k_EdgesTypes %>% mutate(fraction = Freq/10000) %>% arrange(desc(fraction))

Global_Edges_Types <- rbind(Unclass_10k_EdgesTypes[1:5,],
                            GCB_10k_EdgesTypes[1:5,], ABC_10k_EdgesTypes[1:5,])

Global_Edges_Types[1:5,4] <- "ABC"
Global_Edges_Types[6:10,4] <- "GCB"
Global_Edges_Types[11:15,4] <- "Unclass"

Global_Edges_Types$V4 <- as.factor(Global_Edges_Types$V4)
Global_Edges_Types$V4 <- factor(Global_Edges_Types$V4, levels=c('ABC', 'GCB', 'Unclass'))

biotypes_counts <- ggplot(Global_Edges_Types, aes(fill=Var1, x=V4, y=fraction)) + geom_bar(position="stack", stat="identity") + ylim(0,1) +
  scale_fill_brewer(palette = "Paired") +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle("Fractions of Top 5 of interactions between RNA biotypes in the networks") +
  labs(y = "Fraction of initeractions by biotype", x = "Cancer and respective control") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.text = element_text(size=12))


ggsave(filename = "Results/biotypes_counts_DLBCL.tiff", plot = biotypes_counts, units="in", 
       width=10, height=10,  dpi =300)

