# Weighted approach, quantifying the metabolic activities
## Extracting Seurat data
library(Seurat)
library(tidyverse)

seurat <- readRDS(file = 'C:/Users/David/Desktop/Dissertation/Project_Work/Seurat/100720_L4_all_cells_Seurat.rds')

seurat.genes <- seurat@assays$SCT@data
seurat.genes <- as.matrix(seurat.genes) 


## Get a list of genes in a particular pathway
gene.list <- read.table(file = 'C:/Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/Gene_List_With_weight_Normalised.tsv', header = T) 
gene.list <- filter(gene.list, Pathway == 'AminoAcids_Metabolism') # Do this for all 14 pathways
gene.list <- gene.list %>% filter(in_object == T)
ID_Weight <- select(gene.list, 'Gene', 'Percent_Weight')


selected.matrix <- seurat.genes[ID_Weight$Gene,]
selected.df <- as.data.frame(selected.matrix)


weight.df <- selected.df[,1:100955]*ID_Weight$Percent_Weight

Weighted.Expression.AA <- colSums(weight.df[,1:100955], na.rm = F)
## Quality check of scRNA-seq data
DimPlot(subset(seurat, modality != 'n/a'), reduction = "umap", group.by = 'modality')

plot <- VlnPlot(subset(seurat, modality != 'n/a'), group.by = 'modality', feature = "nCount_RNA", pt.size = 0)
plot + scale_y_continuous(trans="log10")

count <- table(subset(seurat, modality != 'n/a')$modality)
barplot(count, main = 'Cell count in different modalities', xlab = 'Modality', ylim = c(0,30000))

VlnPlot(modality, 'pct_counts_Mito', group.by = 'modality', pt.size = 0)


## Loading expression into seurat meta data. This step is done for all 14 pathways and imported back into the meta data.
seurat@meta.data <- mutate(seurat@meta.data, AA_Metab_Expression = Weighted.Expression.AA, .after = birthtime)

# Violin plot of the expression
VlnPlot(subset(seurat, modality != 'n/a'), group.by = 'modality', features = 'Bio_Oxi_Expression', pt.size = 0)
VlnPlot(subset(seurat, modality != 'n/a'), group.by = 'modality', features = 'Carb_Expression', pt.size = 0)
VlnPlot(subset(seurat, modality != 'n/a'), group.by = 'modality', features = 'AA_Metab_Expression', pt.size = 0)

