#' Basic gene enrichment analysis 
#' instruction = past a gene list of differential gene expression to https://maayanlab.cloud/Enrichr/
library(tidyverse)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)
library(ggpubr)

#' ==========================================================================================================
                                         # 3 hour post radiotherapy 
#============================================================================================================
hallmark_1d<- read.csv("hallmark.1d.csv", stringsAsFactors = TRUE)
hallmark_1d<- hallmark_1d[c(1:10),] 
hall<- ggplot(hallmark_1d, aes(reorder(HallMark.Pathway, log.p.adjust), log.p.adjust)) + 
  theme_classic() + geom_bar(stat = "identity", fill= "#FF0000") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (up_genes)", y="-log10(p.adjust)", title = "Hallmark pathway")  +
  scale_fill_gradient(low = "red", high = "blue") +
  coord_flip() 

hall
#' G0 Ontology=== Biological process at 3 hrs
GO_biological_process_1d<- read.csv("GO_Biological_1d.csv", stringsAsFactors = TRUE)
GO_biological_process_1d<- GO_biological_process_1d[c(1:10),] 

# filter top 20 
# biological_process<- GO_biological_process_3hrs[c(1:20),] filter row== good chunk
Bio<-ggplot(GO_biological_process_1d, aes(reorder(GO.Pathways, log_pvalue), log_pvalue)) + 
  theme_classic()+ 
  geom_bar(stat = "identity", fill= "#FF0000")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (up_genes)", y="-log10(p.adjust)", title = "GO Biological Process") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
Bio

#' Go Molecular process 
GO_molecular_process_1d<- read.csv("GO_Molecular_1d.csv", stringsAsFactors = TRUE)

# filter top 20 
Molecular_process<- GO_molecular_process_1d[c(1:20),] 
Molec<- ggplot(Molecular_process, aes(reorder(GO.Pathways, log_pvalue), log_pvalue)) + 
  theme_classic()+ geom_bar(stat = "identity", fill= "#16A085")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 8), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways", y="-log10(p.adjust)", title = "GO Molecular function") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
Molec

#' Downregulated genes 
library(tidyverse)
hallmark_1d_downgene<- read.csv("hallmark_1d_downgene.csv", stringsAsFactors = TRUE)
hallmark_1d_downgene<- hallmark_1d_downgene[c(1:9),]
hallmark_down<- ggplot(hallmark_1d_downgene, aes(reorder(HallMark.Pathway, log.p.adjust), log.p.adjust)) + 
  theme_classic() + geom_bar(stat = "identity", fill= "#0000FF") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (down_genes)", y="-log10(p.adjust)", title = "Hallmark pathway")  +
  scale_fill_gradient(low = "red", high = "blue") +
  coord_flip() 

#' new ggplot
ggplot(hallmark_1d_downgene, aes(reorder(HallMark.Pathway, log.p.adjust), log.p.adjust)) + 
  theme_classic() + geom_bar(stat = "identity", fill= "#0000FF") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (down_genes)", y="-log10(p.adjust)", title = "Hallmark pathawy")  +
  scale_fill_gradient(low = "red", high = "blue") +
  coord_flip() 
hallmark_down

#' GO biological process
GO_biological_process_1d.downgene<- read.csv("GO_biological_process_1d.downgene.csv", stringsAsFactors = TRUE)
GO_biology_1d<- GO_biological_process_1d.downgene[c(1:10),]
Go_down<- ggplot(GO_biology_1d, aes(reorder(GO.Pathways, log_pvalue), log_pvalue)) + 
  theme_classic()+ geom_bar(stat = "identity", fill= "#0000FF")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (down_genes)", y="-log10(p.adjust)", title = "GO Biological Process") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
Go_down

#'arrange plot 
deg<-     ggarrange(hall, hallmark_down, Bio, Go_down, ncol = 2, labels = c("A", "B", "C", "D"), 
          nrow = 2, font.label = list(size=16)) 
deg

deg_down<-     ggarrange(Bio, Go_down, ncol = 1, labels = c("A", "B"), 
                    nrow = 2, font.label = list(size=16)) 

deg_down
ggsave("deg.tiff", plot = deg, height=12, width = 12, units = "in", dpi = 300)

ggarrange(bxp, dp, dens, ncol = 2, nrow = 2,
          labels = c("A", "B", "C"),
          font.label = list(size = 16, color = "red"))


# ========================================================================================================================
                                     # one day post radiotherapy 
# ========================================================================================================================
hallmark_1d<- read.csv("hallmark.1d.csv", stringsAsFactors = TRUE)
hallmark_1d<- hallmark_1d[c(1:10),] 

#' ggplot 
hall<- ggplot(hallmark_1d, aes(reorder(HallMark.Pathway, log.p.adjust), log.p.adjust)) + 
  theme_classic() + geom_bar(stat = "identity", fill= "#FF0000") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (up_genes)", y="-log10(p.adjust)", title = "Hallmark pathway")  +
  scale_fill_gradient(low = "red", high = "blue") +
  coord_flip() 
hall

#' G0 Ontology=== Biological process at 3 hrs
GO_biological_process_1d<- read.csv("GO_Biological_1d.csv", stringsAsFactors = TRUE)
GO_biological_process_1d<- GO_biological_process_1d[c(1:10),] 

#' filter top 20 
# biological_process<- GO_biological_process_3hrs[c(1:20),] filter row== good chunk
Bio<- ggplot(GO_biological_process_1d, aes(reorder(GO.Pathways, log_pvalue), log_pvalue)) + 
  theme_classic()+ 
  geom_bar(stat = "identity", fill= "#FF0000")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (up_genes)", y="-log10(p.adjust)", title = "GO Biological Process") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
Bio

#' Go Molecular process 
GO_molecular_process_1d<- read.csv("GO_Molecular_1d.csv", stringsAsFactors = TRUE)
# filter top 20 
Molecular_process<- GO_molecular_process_1d[c(1:20),] 
Molec<- ggplot(Molecular_process, aes(reorder(GO.Pathways, log_pvalue), log_pvalue)) + 
  theme_classic()+ geom_bar(stat = "identity", fill= "#16A085")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 8), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways", y="-log10(p.adjust)", title = "GO Molecular function") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
Molec

#' Downregulated genes 
library(tidyverse)
hallmark_1d_downgene<- read.csv("hallmark_1d_downgene.csv", stringsAsFactors = TRUE)
hallmark_1d_downgene<- hallmark_1d_downgene[c(1:9),]
hallmark_down<- ggplot(hallmark_1d_downgene, aes(reorder(HallMark.Pathway, log.p.adjust), log.p.adjust)) + 
  theme_classic() + geom_bar(stat = "identity", fill= "#0000FF") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text = element_text(face = "bold", size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (down_genes)", y="-log10(p.adjust)", title = "Hallmark pathway")  +
  scale_fill_gradient(low = "red", high = "blue") +
  coord_flip() 
hallmark_down

#' GO biological process
GO_biological_process_1d.downgene<- read.csv("GO_biological_process_1d.downgene.csv", stringsAsFactors = TRUE)
GO_biology_1d<- GO_biological_process_1d.downgene[c(1:10),]
Go_down<- ggplot(GO_biology_1d, aes(reorder(GO.Pathways, log_pvalue), log_pvalue)) + 
  theme_classic()+ geom_bar(stat = "identity", fill= "#0000FF")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size = 10), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) + 
  labs(x="Pathways (down_genes)", y="-log10(p.adjust)", title = "GO Biological Process") + 
  coord_flip() + 
  scale_fill_gradient(low = "red", high = "blue")
Go_down

deg<-     ggarrange(hall, hallmark_down, Bio, Go_down, ncol = 2, labels = c("A", "B", "C", "D"), 
          nrow = 2, font.label = list(size=16)) 
deg
deg_down<-     ggarrange(Bio, Go_down, ncol = 1, labels = c("A", "B"), 
                    nrow = 2, font.label = list(size=16)) 
deg_down
ggsave("deg.tiff", plot = deg, height=12, width = 12, units = "in", dpi = 300)
