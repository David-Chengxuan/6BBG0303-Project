# Downstream analysis
## Correlation between different pathways
library(tidyverse)

seurat <- readRDS(file = 'C:/Users/David/Desktop/Dissertation/Project_Work/Seurat_Expression_inc.rds')

correlation <- seurat@meta.data[,c("AA_Expression","Bio_Oxi_Expression","Carb_Expression","Energy_Metab_Expression","Inositol_Phosphate_Expression","Lipids_Metabolism_Expression","Mitochondrial_FeS_Biosynthesis_Expression","NO_Metab_Expression","Nucleotide_Metabolism_Expression","Reversible_CO2_Hydration_Expression","Porphyrins_Metabolism_Expression","Pyrophosphate_Hydrolysis_Expression","TCA_Cycle_Respiratory_Electron_transport_Expression","Vitamins_Cofactors_Metabolism_Expression")]

cor <- cor(correlation)

# Plotting correlation
library(ggcorrplot)
ggcorrplot(cor)

# Dendrogram
dd <- dist(scale(cor), method = 'euclidean')
hc <- hclust(dd, method = "ward.D2")
plot(hc)

# Dendrogram in ggcorrplot
ggcorrplot(cor, hc.order = T)

# Draw a dot plot for the selected cells based on cell type
DotPlot(subset(seurat, modality != 'n/a'), features = c("AA_Expression","Bio_Oxi_Expression","Carb_Expression","Energy_Metab_Expression","Inositol_Phosphate_Expression","Lipids_Metabolism_Expression","Mitochondrial_FeS_Biosynthesis_Expression","NO_Metab_Expression","Nucleotide_Metabolism_Expression","Reversible_CO2_Hydration_Expression","Porphyrins_Metabolism_Expression","Pyrophosphate_Hydrolysis_Expression","TCA_Cycle_Respiratory_Electron_transport_Expression","Vitamins_Cofactors_Metabolism_Expression"), group.by = 'modality', scale.by = 'size') + xlab('Pathway') + theme(axis.text.x = element_text(angle = 90))
DotPlot(subset(seurat, neurotransmitter != 'n/a'), features = c("AA_Expression","Bio_Oxi_Expression","Carb_Expression","Energy_Metab_Expression","Inositol_Phosphate_Expression","Lipids_Metabolism_Expression","Mitochondrial_FeS_Biosynthesis_Expression","NO_Metab_Expression","Nucleotide_Metabolism_Expression","Reversible_CO2_Hydration_Expression","Porphyrins_Metabolism_Expression","Pyrophosphate_Hydrolysis_Expression","TCA_Cycle_Respiratory_Electron_transport_Expression","Vitamins_Cofactors_Metabolism_Expression"), group.by = 'neurotransmitter', scale.by = 'size') + xlab('Pathway') + theme(axis.text.x = element_text(angle = 90))