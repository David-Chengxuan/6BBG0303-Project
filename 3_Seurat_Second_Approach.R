#Extracting gene
library(Seurat)
library(tidyverse)


seurat <- readRDS(file = 'C:/Users/David/Desktop/Dissertation/Project_Work/Seurat/100720_L4_all_cells_Seurat.rds')
seurat.genes <- seurat@assays$SCT@data
seurat.genes <- as.matrix(seurat.genes)

# Get a list of genes in a particular pathway
gene.list <- read.table(file = 'C:/Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/Gene_List_With_weight_Normalised.tsv', header = T) 
gene.list <- filter(gene.list, Pathway == 'Vitamins,Cofactors_Metabolism')
gene.list <- gene.list %>% filter(in_object == T)
Pathway.genes <- select(gene.list, 'Gene')

#Get a random set of genes that is not in the selected pathway
all.genes <- rownames(seurat.genes)
Control.genes <- all.genes[!all.genes %in% Pathway.genes$Gene]
Control.genes <- as.data.frame(Control.genes)
Control.genes <- sample_n(Control.genes, 128)


# Check for overlapping
Pathway.genes$Gene %in% Control.genes$Control.genes


Pathway.matrix <- seurat.genes[Pathway.genes$Gene,]
Pathway.df <- as.data.frame(Pathway.matrix)
Mean.Pathway.Exp <- colMeans(Pathway.df, na.rm = T)

Control.matrix <- seurat.genes[Control.genes$Control.genes,]
Control.df <- as.data.frame(Control.matrix)
Mean.Control.Exp <- colMeans(Control.df, na.rm = T)

Vitamins_Cofactors_Metabolism_Expression <- (2^Mean.Pathway.Exp)/(2^Mean.Control.Exp)

seurat@meta.data <- mutate(seurat@meta.data, Vitamins_Cofactors_Metabolism_Expression = Vitamins_Cofactors_Metabolism_Expression, .after = TCA_Cycle_Respiratory_Electron_transport_Expression)

VlnPlot(subset(seurat, modality != 'n/a'), features = 'AA_Expression', pt.size = 0, group.by = 'modality')

length(Bio.Oxi.Expression) == ncol(seurat.genes)  ## if this returns FALSE, it means that some cells are missing from Average1

saveRDS(seurat,file = 'C:/Users/David/Desktop/Dissertation/Project_Work/Seurat_Expression_inc.rds')
