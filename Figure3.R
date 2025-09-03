#Figure 3C, Extended Data 5B

library(ggrepel)
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(clusterProfiler)






brain <- readRDS('CLEANED.merged.clustered.data.rds')
mural <- readRDS("myoVC.mural.clustered.data.rds")


meta.brain <- brain@meta.data
meta.mural <- mural@meta.data


common_cells <- intersect(rownames(meta.brain), rownames(meta.mural))


meta.brain$celltype <- as.character(meta.brain$celltype)


meta.brain[common_cells, "celltype"] <- as.character(meta.mural[common_cells, "celltype"])

meta.brain$celltype[meta.brain$celltype == "myoVC"] <- "Myofibroblast"


meta.brain$celltype <- factor(meta.brain$celltype)


brain@meta.data <- meta.brain

Idents(brain) <- brain$celltype

meta <- brain@meta.data


expr <- brain@assays$RNA$data


meta$ACTA2 <- expr["ACTA2", ]
meta$FN1 <- expr["FN1", ]

brain@meta.data = meta

library(dplyr)
library(tidyr)
library(readr)


meta <- meta %>%
  mutate(celltype = case_when(
    celltype %in% c("Pericyte_1", "Pericyte_2") ~ "Pericyte",
    celltype %in% c("SMC_1", "SMC_2") ~ "SMC",
    TRUE ~ celltype
  ))


pseudobulk_df <- meta %>%
  group_by(subject, celltype) %>%
  summarize(mean_FN1 = mean(FN1, na.rm = TRUE), .groups = "drop")


subjects_with_4types <- pseudobulk_df %>%
  count(subject) %>%
  filter(n >= 4) %>%
  pull(subject)

pseudobulk_df <- pseudobulk_df %>%
  filter(subject %in% subjects_with_4types)


fn1_matrix <- pseudobulk_df %>%
  pivot_wider(names_from = celltype, values_from = mean_FN1)


fn1_matrix_clean <- fn1_matrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))


write_csv(fn1_matrix_clean, "FN1_expression_subjects_with_4_celltypes.csv", na = "")







brain=readRDS('CLEANED.merged.clustered.data.rds')
meta=brain@meta.data

meta.yang <- subset(meta, subset = dataset == "yang")
meta.sun <- subset(meta, subset = dataset == "sun")

yang.subjects <- unique(meta.yang$subject)
sun.subjects <- unique(meta.sun$subject)



brain <- readRDS("myoVC.mural.clustered.data.rds")

brain$celltype[brain$celltype == "myoVC"] <- "Myofibroblast"


Idents(brain) <- brain$celltype


meta = brain@meta.data

meta.myofibs <- subset(meta, subset = celltype == "Myofibroblast")

meta.myofibs.yang <- subset(meta.myofibs, subset = dataset == "yang")
meta.myofibs.sun <- subset(meta.myofibs, subset = dataset == "sun")


myofibs.list.yang <- rownames(meta.myofibs.yang)
myofibs.list.sun <- rownames(meta.myofibs.sun)





brain=readRDS('Yang.clustered.all.data.rds')
meta = brain@meta.data


brain <- subset(brain, subset = subject %in% yang.subjects)
meta = brain@meta.data

brain$dataset <- "yang"


DimPlot(brain, label = T)
FeaturePlot(brain, features = "ABCA10")


Idents(brain) <- "seurat_clusters"


brain$celltype <- as.character(Idents(brain))  


brain$celltype[Idents(brain) == 8] <- "Fibroblast"


brain$celltype[rownames(brain@meta.data) %in% myofibs.list.yang] <- "Myofibroblast"


Idents(brain) <- brain$celltype

brain.yang <- subset(brain, subset = celltype == "Fibroblast" | celltype == "Myofibroblast")

DimPlot(brain.yang)






brain <- readRDS("ROSMAP.VascularCells.seurat.harmony.final.rds")
meta <- brain@meta.data


brain <- subset(brain, subset = subject %in% sun.subjects)
meta.sun <- brain@meta.data


meta.sun$dataset <- "sun"
colnames(meta.sun)[colnames(meta.sun) == "age_death"] <- "age"


colnames(meta.sun)[colnames(meta.sun) == "msex"] <- "sex"


meta.sun$sex <- ifelse(meta.sun$sex == "Male", "M",
                       ifelse(meta.sun$sex == "Female", "F", meta.sun$sex))


colnames(meta.sun)[colnames(meta.sun) == "apoe_genotype"] <- "apoe"


meta.sun$apoe <- ifelse(meta.sun$apoe == "33", "e3",
                        ifelse(meta.sun$apoe %in% c("34", "44"), "e4", meta.sun$apoe))

colnames(meta.sun)[colnames(meta.sun) == "ADdiag2types"] <- "AD.status"


brain@meta.data = meta.sun



brain$celltype[rownames(brain@meta.data) %in% myofibs.list.sun] <- "Myofibroblast"


Idents(brain) <- brain$celltype

brain.sun <- subset(brain, subset = celltype == "Fib" | celltype == "Myofibroblast")

DimPlot(brain.sun)



seurat_list <- list(brain.sun, brain.yang)



brain <- Reduce(merge, seurat_list)
meta=brain@meta.data


brain[["RNA"]] <- JoinLayers(brain[["RNA"]])
meta = brain@meta.data

brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain)

brain$celltype[brain$celltype == "Fib"] <- "Fibroblast"
brain$celltype <- factor(brain$celltype)  
Idents(brain) <- brain$celltype

meta = brain@meta.data



expr <- brain@assays$RNA$data


meta$ACTA2 <- expr["ACTA2", ]
meta$FN1 <- expr["FN1", ]

brain@meta.data = meta

library(dplyr)
library(tidyr)
library(readr)




pseudobulk_df <- meta %>%
  group_by(subject, celltype) %>%
  dplyr::summarize(mean_FN1 = mean(FN1, na.rm = TRUE), .groups = "drop")


subjects_with_4types <- pseudobulk_df %>%
  dplyr::count(subject) %>%
  dplyr::filter(n >= 2) %>%
  dplyr::pull(subject)

pseudobulk_df <- pseudobulk_df %>%
  filter(subject %in% subjects_with_4types)


fn1_matrix <- pseudobulk_df %>%
  pivot_wider(names_from = celltype, values_from = mean_FN1)


fn1_matrix_clean <- fn1_matrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "", .)))


write_csv(fn1_matrix_clean, "FN1_expression_myofibs_fibs.csv", na = "")









brain.yang <- subset(brain.yang, subset = celltype == "Fibroblast")

brain.sun <- subset(brain.sun, subset = celltype == "Fib")

seurat_list <- list(brain.sun, brain.yang)


brain <- Reduce(merge, seurat_list)
meta=brain@meta.data


brain[["RNA"]] <- JoinLayers(brain[["RNA"]])
meta = brain@meta.data



brain$celltype[brain$celltype == "Fib"] <- "Fibroblast"
brain$celltype <- factor(brain$celltype)  
Idents(brain) <- brain$celltype

meta = brain@meta.data


brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain)


Idents(brain) <- brain$apoe

result <- FindMarkers(object = brain, 
                      ident.1 = "e4", 
                      ident.2 = "e3",
                      test.use = "MAST",
                      latent.vars = c("age", "sex", "AD.status", "dataset")) 

result$gene <- rownames(result)

write.csv(result, file = "fibroblast.apoe.MAST.csv")


library(msigdbr)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(tibble)
library(readxl)
library(dplyr)
library(fgsea)

set.seed(42) 

cgp_gene_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
unique_gene_sets = split(x = cgp_gene_sets$gene_symbol, f = cgp_gene_sets$gs_name)



all_GSEA_results <- data.frame()





mast <- result


rownames(mast) <- mast$gene


transformed_value <- -log10(mast$p_val)
transformed_value <- sign(mast$avg_log2FC) * abs(transformed_value)


new_mast <- data.frame(transformed_value, row.names = rownames(mast))


named_vector <- deframe(new_mast)
named_vector <- setNames(named_vector, rownames(mast))
named_vector <- sort(named_vector, decreasing = TRUE)


GSEA_result <- fgseaMultilevel(unique_gene_sets, named_vector)


GSEA_result <- GSEA_result %>%
  mutate(across(where(is.list), ~ sapply(., toString)))


all_GSEA_results <- rbind(all_GSEA_results, GSEA_result)


write.csv(all_GSEA_results, file = "fibroblast.MAST.GSEA.csv")


#######################################################################################


#Extended Data 5C





brain=readRDS('CLEANED.merged.clustered.data.rds')
meta=brain@meta.data


brain <- subset(brain, subset = dataset == "haney")


DimPlot(brain)

meta = brain@meta.data

brain <- subset(brain, subset = celltype == "Pericyte" | celltype == "SMC")

meta = brain@meta.data


library(readxl)
haney.meta <- read_excel('haney.meta.xlsx')


library(dplyr)


meta <- meta %>%
  left_join(haney.meta, by = c("subject" = "CaseID"))

brain@meta.data =  meta





expr <- brain[["RNA"]]$data


meta <- brain@meta.data



fn1_expr <- expr["FN1", ]  


meta$FN1 <- fn1_expr




fn1_summary <- meta %>%
  group_by(subject) %>%
  dplyr::summarize(
    mean_FN1 = mean(FN1, test = TRUE),
    age.at.diagnosis = unique(neurological_dx_age)  
  )

library(ggplot2)


meta_subset <- meta[, c("subject", "age", "sex", "apoe")]  



cor.test(fn1_summary$mean_FN1, fn1_summary$age.at.diagnosis, method = "pearson")



p <- ggplot(fn1_summary, aes(x = mean_FN1, y = age.at.diagnosis)) +  
  geom_point(size = 3, shape = 21, fill = "black") +  
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", fill = "#0072B2", alpha = 0.2) +
  labs(
    x = "Mean mural cell FN1 Expression",
    y = "Age at AD Diagnosis",
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "plain"),     
    axis.text = element_text(face = "plain")       
  )


ggsave("fn1.diagnosis.age.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")


meta$global.pathology.score <- rowMeans(meta[, c("sum_lb_density", "PlaqueTotal", "TangleTotal", "infarct_cerebral_total_volume")], na.rm = TRUE)


library(dplyr)
library(ggplot2)


meta_clean <- meta %>%
  filter(!is.na(global.pathology.score), !is.na(FN1))




q75 <- quantile(meta_clean$global.pathology.score, probs = 0.75, na.rm = TRUE)


meta_clean$global.pathology.score_bin <- ifelse(meta_clean$global.pathology.score > q75, "Top 25%", "Lower 75%")



fn1_summary_binned <- meta_clean %>%
  dplyr::group_by(subject, global.pathology.score_bin, age, sex, apoe) %>%
  dplyr::summarize(mean_FN1 = mean(FN1), .groups = "drop")


write.csv(fn1_summary_binned, file = 'fn1.pathology.csv')








