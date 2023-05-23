# Compiling reactome

library(tidyverse)
library(vroom)


metalist <- list.files("C://Users/David/Desktop/Dissertation/Project_Work/Metabolism_Worm_Genesets", recursive = T, pattern = ".tsv$",full.names = T)
metadf <- vroom(metalist, comment = '#', id = 'file')

nrow(metadf)
ncol(metadf)
head(metadf)

simple.metadf <- select(metadf, 'file', Name = '...1')
simple.metadf$file <- str_remove(simple.metadf$file, "Metabolism_Worm_Genesets/")
simple.metadf$file <- str_remove(simple.metadf$file, ".tsv$")
final.metadf <- separate(simple.metadf, col = file, into = c("Pathway", "Group", "Subgroup"), sep = "/")

# Assigning Gene ID

df <- read.table(file = 'C://Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/GeneID.tsv', fill = T, header = T) 
df2 <- read.table(file = 'C://Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/Metabolic_Genes_Full.tsv', header = T) 

df <- select(df, Gene, Name = ID)
df3 <- left_join(df2, df, by= "Name")

filter(df3, Gene %in% NA) # fill some gene ID manually

df3[148,5] <- 'WBGene00015540'
df3[226,5] <- 'WBGene00009186'
df3[227,5] <- 'WBGene00046697'
df3[228,5] <- 'WBGene00000096'
df3[231,5] <- 'WBGene00000095'
df3[312,5] <- 'WBGene00045379'
df3[323,5] <- 'WBGene00018758'
df3[324,5] <- 'WBGene00018925'
df3[325,5] <- 'WBGene00009049'
df3[327,5] <- 'WBGene00020268'
df3[1222,5] <- 'WBGene00000281'
df3[1224,5] <- 'WBGene00000284'
df3[1226,5] <- 'WBGene00000281'
df3[1228,5] <- 'WBGene00000284'
df3[1230,5] <- 'WBGene00000281'
df3[1232,5] <- 'WBGene00000284'
df3[1234,5] <- 'WBGene00000281'
df3[1236,5] <- 'WBGene00000284'

# Extracting genes from SCT assay
library(Seurat)
seurat <- readRDS(file = "C://Users/David/Desktop/Dissertation/Project_Work/Seurat/100720_L4_all_cells_Seurat.rds")
Genes <- rownames(seurat@assays$SCT@counts)

# Testing if the genes are found in the seurat object
df3$Gene %in% Genes
df3 <- mutate(df3, in_object = Gene %in% Genes)

write.table(df3, file = 'C://Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/Annotated_Gene_List_with_Seurat_Existence.tsv', quote = F, row.names = F)

# Assigning gene weight
metalist <- list.files("C://Users/David/Desktop/Dissertation/Project_Work/Metabolism_Worm_Genesets", recursive = T, pattern = ".tsv$",full.names = T)
metadf <- vroom(metalist, comment = '#', id = 'file')

# Emerge the columns and convert it into a matrix, Run the first 7 lines, then convert to matrix
# Row name = Gene, Column = Larval Stages
# Find the weight of the genes in different pathways

calculation.df <- select(metadf, -c('...1','file'))
matrix <- data.matrix(calculation.df)
Mean.Weight <- rowMeans(matrix, na.rm = T)
Mean.Weight <- round(Mean.Weight, digits = 2)

full.list <- read.table(file = 'C://Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/Annotated_Gene_List_with_Seurat_Existence.tsv', header = T)
full.list <- mutate(full.list, Weight= Mean.Weight)

full.list <- full.list %>% group_by(Pathway) %>% mutate(Percent_Weight = Weight/sum(Weight)) 

write.table(full.list,file = 'C://Users/David/Desktop/Dissertation/Project_Work/Gene_List_with_ID/Gene_List_With_weight_Normalised.tsv')




