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
library(MatchIt)

set.seed(123)

#Figure 1A, Extended Data 1A, B

brain=readRDS('Yang.clustered.all.data.rds')

brain <- subset(brain, subset = seurat_clusters == 2 | seurat_clusters == 5 | seurat_clusters == 3 | seurat_clusters == 11 | seurat_clusters == 13 | seurat_clusters == 15 | seurat_clusters == 0 | seurat_clusters == 9)

cluster_to_celltype <- c(
  "0" = "Endo",
  "9" = "Endo",
  "2" = "Mural",
  "5" = "Mural",
  "11" = "Astro",
  "3" = "Astro",
  "13" = "Astro",
  "15" = "Astro"
)

brain$celltype <- recode(brain$seurat_clusters, !!!cluster_to_celltype)

brain$celltype[is.na(brain$celltype)] <- "NA"

table(brain$celltype)

Idents(brain) <- brain$celltype

brain.yang <- brain
meta.yang <- brain.yang@meta.data


brain.yang@meta.data = meta.yang
meta.yang = brain.yang@meta.data

Idents(brain.yang) <- meta.yang$celltype

meta.yang$dataset <- "yang"
meta.yang$brain_region <- "Hippocampus"
meta.yang$batch <- "yang"
brain.yang@meta.data = meta.yang


brain <- brain.yang

cell_counts <- brain@meta.data %>%
  group_by(subject) %>%
  summarise(total_cells = n(), .groups = "drop")

yang_subjects <- brain@meta.data %>%
  distinct(subject, apoe, age, sex, AD.status) %>%
  left_join(cell_counts, by = "subject")

yang_subjects <- yang_subjects %>%
  mutate(
    apoe_binary = ifelse(apoe == "e3", 1, 0),  
    sex_binary = ifelse(sex == "M", 1, 0),     
    AD_status_binary = ifelse(AD.status == "AD", 1, 0)  
  )

yang_subjects <- yang_subjects %>%
  mutate(across(c(apoe_binary, sex_binary, AD_status_binary, total_cells, age), as.numeric))

table(yang_subjects$apoe_binary)



match_model <- matchit(apoe_binary ~ age + sex_binary + AD_status_binary + total_cells, 
                       data = yang_subjects, 
                       method = "nearest",  
                       ratio = 1)  

matched_data <- match.data(match_model)
matched_subjects <- as.character(matched_data$subject[matched_data$apoe_binary == 0])  

valid_e3_subjects <- as.character(yang_subjects$subject[yang_subjects$apoe == "e3"])

if (!"subject" %in% colnames(brain@meta.data)) {
  stop("Error: 'subject' column not found in brain@meta.data. Check column names.")
}

brain <- subset(brain, 
                cells = rownames(brain@meta.data[brain@meta.data$subject %in% c(valid_e3_subjects, matched_subjects), ]))



brain.yang <- brain






brain.sun <- readRDS("ROSMAP.VascularCells.seurat.harmony.final.rds")
meta.sun <- brain.sun@meta.data


brain.sun <- subset(brain.sun, cells = which(!is.na(brain.sun$subject)))
meta.sun = brain.sun@meta.data


brain.sun <- subset(brain.sun, subset = apoe_genotype == "33" | apoe_genotype == "34" | apoe_genotype == "44")
meta.sun = brain.sun@meta.data


DimPlot(brain.sun, label = T)

brain.sun <- subset(brain.sun, subset = celltype == "Endo" | celltype == "SMC" | celltype == "Per")
meta.sun = brain.sun@meta.data

Idents(brain.sun) <- meta.sun$celltype


meta.sun$dataset <- "sun"
colnames(meta.sun)[colnames(meta.sun) == "age_death"] <- "age"

colnames(meta.sun)[colnames(meta.sun) == "msex"] <- "sex"

meta.sun$sex <- ifelse(meta.sun$sex == "Male", "M",
                       ifelse(meta.sun$sex == "Female", "F", meta.sun$sex))

colnames(meta.sun)[colnames(meta.sun) == "apoe_genotype"] <- "apoe"

meta.sun$apoe <- ifelse(meta.sun$apoe == "33", "e3",
                        ifelse(meta.sun$apoe %in% c("34", "44"), "e4", meta.sun$apoe))

colnames(meta.sun)[colnames(meta.sun) == "ADdiag2types"] <- "AD.status"


brain.sun@meta.data = meta.sun

meta.sun = brain.sun@meta.data


brain <- brain.sun


cell_counts <- brain@meta.data %>%
  group_by(subject) %>%
  summarise(total_cells = n(), .groups = "drop")

sun_subjects <- brain@meta.data %>%
  distinct(subject, apoe, age, sex, AD.status) %>%
  left_join(cell_counts, by = "subject")

sun_subjects <- sun_subjects %>%
  mutate(
    apoe_binary = ifelse(apoe == "e4", 1, 0),  
    sex_binary = ifelse(sex == "M", 1, 0),     
    AD_status_binary = ifelse(AD.status == "AD", 1, 0)  
  )

sun_subjects <- sun_subjects %>%
  mutate(across(c(apoe_binary, sex_binary, AD_status_binary, total_cells, age), as.numeric))

table(sun_subjects$apoe_binary) 

set.seed(123)


match_model <- matchit(apoe_binary ~ age + sex_binary + AD_status_binary + total_cells, 
                       data = sun_subjects, 
                       method = "nearest",  
                       ratio = 1) 

matched_data <- match.data(match_model)
matched_subjects <- as.character(matched_data$subject[matched_data$apoe_binary == 1])  

valid_e3_subjects <- as.character(matched_data$subject[matched_data$apoe_binary == 0])

if (!"subject" %in% colnames(brain@meta.data)) {
  stop("Error: 'subject' column not found in brain@meta.data. Check column names.")
}

brain <- subset(brain, 
                cells = rownames(brain@meta.data[brain@meta.data$subject %in% c(valid_e3_subjects, matched_subjects), ]))


brain.sun <- brain


meta.sun <- brain.sun@meta.data







brain.haney <- readRDS("Haney.clustered.adonly.data.rds")
meta.haney <- brain.haney@meta.data
DimPlot(brain.haney, label = T)

brain.haney <- subset(brain.haney, subset = celltype == "Vascular" | celltype == "Astro")
meta.haney = brain.haney@meta.data

meta.haney$dataset <- "haney"

colnames(meta.haney)[colnames(meta.haney) == "ADstatus"] <- "AD.status"

colnames(meta.haney)[colnames(meta.haney) == "sample"] <- "subject"

meta.haney$brain_region <- "Prefrontal_cortex"
meta.haney$batch <- "haney"

brain.haney@meta.data = meta.haney

brain <- brain.haney

cell_counts <- brain@meta.data %>%
  group_by(subject) %>%
  summarise(total_cells = n(), .groups = "drop")


yang_subjects <- brain@meta.data %>%
  distinct(subject, apoe, age, sex, AD.status) %>%
  left_join(cell_counts, by = "subject")

yang_subjects <- yang_subjects %>%
  mutate(
    apoe_binary = ifelse(apoe == "e3", 1, 0),
    sex_binary = ifelse(sex == "M", 1, 0),     
    AD_status_binary = ifelse(AD.status == "AD", 1, 0)
  )

yang_subjects <- yang_subjects %>%
  mutate(across(c(apoe_binary, sex_binary, AD_status_binary, total_cells, age), as.numeric))


table(yang_subjects$apoe_binary)  

set.seed(123)

match_model <- matchit(apoe_binary ~ age + sex_binary + AD_status_binary + total_cells, 
                       data = yang_subjects, 
                       method = "nearest",  
                       ratio = 1)  

matched_data <- match.data(match_model)
matched_subjects <- as.character(matched_data$subject[matched_data$apoe_binary == 0])  


valid_e3_subjects <- as.character(yang_subjects$subject[yang_subjects$apoe == "e3"])

if (!"subject" %in% colnames(brain@meta.data)) {
  stop("Error: 'subject' column not found in brain@meta.data. Check column names.")
}

brain <- subset(brain, 
                cells = rownames(brain@meta.data[brain@meta.data$subject %in% c(valid_e3_subjects, matched_subjects), ]))




brain.haney <- brain
meta.haney <- brain.haney@meta.data




brain.yang_individuals <- brain.yang@meta.data %>%
  dplyr::group_by(subject) %>%
  dplyr::summarize(
    apoe = dplyr::first(apoe),
    age = mean(age),
    sex = dplyr::first(sex),
    AD.status = dplyr::first(AD.status),  
    avg_cell_contrib = n()  
  )

brain.haney_individuals <- brain.haney@meta.data %>%
  dplyr::group_by(subject) %>%
  dplyr::summarize(
    apoe = dplyr::first(apoe),
    age = mean(age),
    sex = dplyr::first(sex),
    AD.status = dplyr::first(AD.status),
    avg_cell_contrib = n()  
  )

brain.sun_individuals <- brain.sun@meta.data %>%
  dplyr::group_by(subject) %>%
  dplyr::summarize(
    apoe = dplyr::first(apoe),
    age = mean(age),
    sex = dplyr::first(sex),
    AD.status = dplyr::first(AD.status),
    avg_cell_contrib = n()  
  )



yang_results <- list()
yang_results$age <- t.test(age ~ apoe, data = brain.yang_individuals)
yang_results$sex <- fisher.test(table(brain.yang_individuals$sex, brain.yang_individuals$apoe))
yang_results$AD.status <- fisher.test(table(brain.yang_individuals$AD.status, brain.yang_individuals$apoe))
yang_results$avg_cell_contrib <- t.test(avg_cell_contrib ~ apoe, data = brain.yang_individuals)

haney_results <- list()
haney_results$age <- t.test(age ~ apoe, data = brain.haney_individuals)
haney_results$sex <- fisher.test(table(brain.haney_individuals$sex, brain.haney_individuals$apoe))
haney_results$avg_cell_contrib <- t.test(avg_cell_contrib ~ apoe, data = brain.haney_individuals)

sun_results <- list()
sun_results$age <- t.test(age ~ apoe, data = brain.sun_individuals)
sun_results$sex <- chisq.test(table(brain.sun_individuals$sex, brain.sun_individuals$apoe))
sun_results$AD.status <- chisq.test(table(brain.sun_individuals$AD.status, brain.sun_individuals$apoe))
sun_results$avg_cell_contrib <- t.test(avg_cell_contrib ~ apoe, data = brain.sun_individuals)



library(dplyr)
library(knitr)
library(kableExtra)
library(tibble)

format_age <- function(x) {
  sprintf("%.0f (%.0f–%.0f)", mean(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

format_n_pct <- function(vec, value = NULL) {
  if (!is.null(value)) vec <- vec == value
  n <- sum(vec, na.rm = TRUE)
  pct <- round(100 * n / sum(!is.na(vec)))
  sprintf("%d (%.0f%%)", n, pct)
}

format_p <- function(p) {
  if (p < 0.001) return("<0.001")
  else sprintf("%.3f", p)
}

make_summary_table <- function(df, results) {
  df3 <- df %>% filter(apoe == "e3")
  df4 <- df %>% filter(apoe == "e4")
  
  values <- list(
    Characteristic = c("Mean age (range)",
                       "Female, n (%)",
                       "AD diagnosis, n (%)",
                       "Average cell yield, mean (range)"),
    
    `APOE3/3` = c(
      format_age(df3$age),
      format_n_pct(df3$sex, "F"),
      format_n_pct(df3$AD.status, "AD"),
      format_age(df3$avg_cell_contrib)
    ),
    
    `APOE4 carriers` = c(
      format_age(df4$age),
      format_n_pct(df4$sex, "F"),
      format_n_pct(df4$AD.status, "AD"),
      format_age(df4$avg_cell_contrib)
    ),
    
    `p-value` = c(
      format_p(results$age$p.value),
      format_p(results$sex$p.value),
      if (!is.null(results$AD.status)) format_p(results$AD.status$p.value) else NA,
      format_p(results$avg_cell_contrib$p.value)
    )
  )
  
  as.data.frame(values, stringsAsFactors = FALSE)
}

add_dataset_label <- function(label, n_e3, n_e4) {
  tibble(
    Characteristic = paste0(label, " (n=", n_e3, " APOE3/3, n=", n_e4, " APOE4)"),
    `APOE3/3` = "",
    `APOE4 carriers` = "",
    `p-value` = ""
  )
}

yang_summary  <- make_summary_table(brain.yang_individuals, yang_results)
haney_summary <- make_summary_table(brain.haney_individuals, haney_results)
sun_summary   <- make_summary_table(brain.sun_individuals, sun_results)

colnames(yang_summary) <- colnames(haney_summary) <- colnames(sun_summary) <-
  c("Characteristic", "APOE3/3", "APOE4 carriers", "p-value")

yang_summary[]  <- lapply(yang_summary, as.character)
haney_summary[] <- lapply(haney_summary, as.character)
sun_summary[]   <- lapply(sun_summary, as.character)

combined_summary <- bind_rows(
  add_dataset_label("Yang Dataset",  nrow(filter(brain.yang_individuals,  apoe == "e3")), nrow(filter(brain.yang_individuals,  apoe == "e4"))),
  yang_summary,
  add_dataset_label("Haney Dataset", nrow(filter(brain.haney_individuals, apoe == "e3")), nrow(filter(brain.haney_individuals, apoe == "e4"))),
  haney_summary,
  add_dataset_label("Sun Dataset",   nrow(filter(brain.sun_individuals,   apoe == "e3")), nrow(filter(brain.sun_individuals,   apoe == "e4"))),
  sun_summary
)

yang_row  <- 1
haney_row <- yang_row + nrow(yang_summary) + 1
sun_row   <- haney_row + nrow(haney_summary) + 1


kable_out <- kable(combined_summary, align = "lccr", escape = FALSE,
                   caption = "Summary of APOE3/3 vs APOE4 carriers across datasets") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed")) %>%
  column_spec(1, bold = TRUE) %>%
  add_header_above(c(" " = 1, "APOE Genotype" = 2, " " = 1)) %>%
  row_spec(0, extra_css = "border-bottom: 2px solid black;") %>%  
  row_spec(yang_row,  bold = TRUE, background = "#f2f2f2", extra_css = "border-top: 2px solid black; border-bottom: 2px solid black;") %>%
  row_spec(haney_row, bold = TRUE, background = "#f2f2f2", extra_css = "border-top: 2px solid black; border-bottom: 2px solid black;") %>%
  row_spec(sun_row,   bold = TRUE, background = "#f2f2f2", extra_css = "border-top: 2px solid black; border-bottom: 2px solid black;")

kable_out


save_kable(kable_out, file = "demographics.html")
webshot2::webshot("demographics.html", "demographics.pdf", zoom = 2, vwidth = 1200, vheight = 800, selector = "table")


############################################################################################

#Figure 1A, B, Extended Data 1C, D


seurat_list <- list(brain.haney, brain.sun, brain.yang)


brain <- Reduce(merge, seurat_list)
meta=brain@meta.data


brain[["RNA"]] <- JoinLayers(brain[["RNA"]])
meta = brain@meta.data


id_data='merged.processed.filtered.data'
saveRDS(brain,file=paste(id_data,'.rds',sep=''))


brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain)

brain <- RunPCA(brain, features = VariableFeatures(object = brain))
DimPlot(brain, reduction = "pca")
ElbowPlot(brain,ndims = 50)
k=1:15


brain <- brain %>%
  RunHarmony(c("subject", "dataset"), plot_convergence = TRUE, theta = c(2, 4))


brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

meta = brain@meta.data

DimPlot(brain, label = T)


id_data='uncleaned.merged.clustered.data'
saveRDS(brain,file=paste(id_data,'.rds',sep=''))


brain = readRDS("uncleaned.merged.clustered.data.rds")
meta = brain@meta.data

FeaturePlot(brain, features = "MBP", label = T) #8, 7 has high expression of oligo markers
FeaturePlot(brain, features = "PDGFRB", label = T) #11 has high expression of peri + astro markers
FeaturePlot(brain, features = "PECAM1", label = T) #12 has high expression of astro + endo markers, #14 isn't highly expressing vascular markers


brain <- subset(brain, subset = seurat_clusters == 8 | seurat_clusters == 7 | seurat_clusters == 11 | seurat_clusters == 12 | seurat_clusters == 14, invert = T)


meta = brain@meta.data

DimPlot(brain, label = T)



brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain)

brain <- RunPCA(brain, features = VariableFeatures(object = brain))
DimPlot(brain, reduction = "pca")
ElbowPlot(brain,ndims = 50)
k=1:15



brain <- brain %>% 
  RunHarmony(c("subject", "dataset"), plot_convergence = TRUE, theta = c(2, 4))


brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()

meta = brain@meta.data

DimPlot(brain, label = T)

brain <- subset(brain, subset = seurat_clusters == 7, invert = T) #Lower pericyte marker expression


meta = brain@meta.data

DimPlot(brain, label = T)


brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
brain <- ScaleData(brain)

brain <- RunPCA(brain, features = VariableFeatures(object = brain))
DimPlot(brain, reduction = "pca")
ElbowPlot(brain,ndims = 50)
k=1:15



brain <- brain %>% 
  RunHarmony(c("subject", "dataset"), plot_convergence = TRUE, theta = c(2, 4))


brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()

meta = brain@meta.data

DimPlot(brain, label = T)


p <- DimPlot(brain, split.by = "dataset")
ggsave("Integration.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



FeaturePlot(brain, features = "PECAM1", label = T) #2, 3, 4 = endos
FeaturePlot(brain, features = "PDGFRB", label = T) #1 = pericytes
FeaturePlot(brain, features = "MYH11", label = T) #5 = SMCs
FeaturePlot(brain, features = "AQP4", label = T) #0= astros

p <- FeaturePlot(brain, features = c("PECAM1", "PDGFRB", "MYH11", "AQP4"), label = F) #2, 3, 4 = endos
ggsave("Clusterannotations.svg", plot = p, width = 6, height = 4.5, dpi = 300, device = "svg")



DotPlot(brain, features = c("PECAM1", "PDGFRB", "MYH11", "AQP4")) + RotatedAxis()


brain$celltype <- "Unknown" 

brain$celltype[brain$seurat_clusters %in% c("2", "3", "4")] <- "Endothelial"
brain$celltype[brain$seurat_clusters %in% c("1")] <- "Pericyte"
brain$celltype[brain$seurat_clusters %in% c("5")] <- "SMC"
brain$celltype[brain$seurat_clusters %in% c("0")] <- "Astrocyte"

brain$celltype <- factor(brain$celltype)

table(brain$celltype)

DimPlot(brain, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()

Idents(brain) <- brain$celltype


results <- FindAllMarkers(brain, min.pct = 0.25, only.pos = T)
write.csv(results, file = "major.cluster.markers.csv")


colors <- c("#006884", "#058799", "#E74C2E", "#EE9E80")

p <- DimPlot(brain, cols = colors)


ggsave("UMAP.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



id_data='CLEANED.merged.clustered.data'
saveRDS(brain,file=paste(id_data,'.rds',sep=''))

####################################################################################

#Figure 1C


rstudioapi::restartSession()

library(ggrepel)
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(destiny)
library(slingshot); library(SingleCellExperiment)
library(RColorBrewer); library(scales)
library(viridis); library(UpSetR)
library(pheatmap); library(msigdbr)
library(fgsea); library(knitr)
library(ggplot2); library(gridExtra)
library(tradeSeq); library(Seurat)
library(plot3D); library(fields); library(rgl)
library(DAseq); library(SeuratExtend)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(ggrepel)
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)

brain=readRDS('CLEANED.merged.clustered.data.rds')
meta=brain@meta.data

brain_sce <- as.SingleCellExperiment(brain)
brain_milo <- Milo(brain_sce)

brain_milo <- buildGraph(brain_milo, k = 50, d = 30, reduced.dim = "HARMONY")
brain_milo <- makeNhoods(brain_milo, prop = 0.1, k = 50, d = 30, refined = TRUE, reduced_dims = "HARMONY")
plotNhoodSizeHist(brain_milo)
brain_milo <- countCells(brain_milo, meta.data = data.frame(colData(brain_milo)), sample="subject")



brain_milo <- calcNhoodDistance(brain_milo, d=30, reduced.dim = "HARMONY")

id_data='CLEANED.merged.clustered.miloR'
saveRDS(brain_milo,file=paste(id_data,'.rds',sep=''))




brain_milo=readRDS('CLEANED.merged.clustered.miloR.rds')





brain_design <- data.frame(colData(brain_milo))[,c("subject", "age", "sex", "AD.status", "apoe", "dataset")]

brain_design <- distinct(brain_design)
rownames(brain_design) <- brain_design$subject


da_results <- testNhoods(brain_milo, reduced.dim = "HARMONY", design = ~ age + sex + AD.status + dataset + apoe, design.df = brain_design)




ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) 

da_results <- annotateNhoods(brain_milo, da_results, coldata_col = "celltype")

plotDAbeeswarm(da_results, group.by = "celltype")

sig <- subset(da_results, subset = SpatialFDR < 0.05)

write.csv(sig, file = "MiloR.sig.neighborhoods.csv")




library(ggplot2)
library(ggbeeswarm)
library(colorspace)

darken <- function(color) adjustcolor(color, alpha.f = 0.8)

lim.val <- max(abs(da_results$logFC))
eps <- lim.val / 25

significant_results <- da_results[da_results$SpatialFDR <= 0.05, ]

mean_logFC <- aggregate(logFC ~ celltype, data = significant_results, mean)


wilcox_results <- aggregate(logFC ~ celltype, data = significant_results, function(x) wilcox.test(x, mu = 0)$p.value)
wilcox_results$adj_p <- p.adjust(wilcox_results$logFC, method = "BH")

wilcox_results$significance <- cut(wilcox_results$adj_p, 
                                   breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                   labels = c("****", "***", "**", "*", "n.s."), 
                                   right = FALSE)  

facet_labels <- setNames(
  paste0(wilcox_results$celltype, "\n(", wilcox_results$significance, ")"), 
  wilcox_results$celltype
)

da_bee <- ggplot(mapping = aes(x = celltype, y = logFC)) +
  geom_quasirandom(data = significant_results, aes(colour = logFC), size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size = 1) +
  geom_hline(data = mean_logFC, aes(yintercept = logFC), linetype = "solid", colour = "black", size = 2) +
  coord_flip() +
  facet_grid(celltype ~ ., scales = "free_y", space = "free", labeller = as_labeller(facet_labels)) +  
  scale_colour_gradient2(low = darken('blue'), mid = 'grey80', high = darken('#ff0000'),
                         midpoint = 0) +
  scale_y_continuous(limits = c(-lim.val - eps, lim.val + eps)) +
  theme_bw() +
  labs(y = "Log Fold Change") +
  theme(strip.text.y = element_text(angle = 0, colour = 'black', size = 20),  
        strip.background = element_rect(colour = 'white', fill = 'white'),
        axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  guides(colour = FALSE)



da_bee


p <- da_bee

ggsave("DAbee.svg", plot = p, width = 5, height = 4, dpi = 300, device = "svg")



brain_milo <- buildNhoodGraph(brain_milo)

test <- brain_milo


library(igraph)
library(patchwork)

significance_threshold <- 0.05

significant_nhoods <- which(da_results$SpatialFDR < significance_threshold)

da_results_filtered <- da_results[significant_nhoods, ]

graph <- brain_milo@nhoodGraph[[1]]  
filtered_graph <- induced_subgraph(graph, vids = significant_nhoods)
test@nhoodGraph[[1]] <- filtered_graph

umap_pl <- plotReducedDim(test, dimred = "UMAP", colour_by="celltype", text_by = "celltype", text_size = 3) +
  guides(fill="none")

nh_graph_pl <- plotNhoodGraphDA(test, da_results_filtered, layout="UMAP", alpha=0.3) +
  scale_size_continuous(range = c(1, 5))


p <- nh_graph_pl

ggsave("Neighborhoods.graph.svg", plot = p, width = 4, height = 4, dpi = 300, device = "svg")



combined_plot <- umap_pl + nh_graph_pl + plot_layout(guides = "collect")


ggsave("combined_umap_nh_graph.pdf", plot = combined_plot, width = 12, height = 10)

#############################################################################################


#Figure 1F, Extended Data 1E


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




brain=readRDS('CLEANED.merged.clustered.data.rds')
meta=brain@meta.data

brain <- subset(brain, subset = celltype == "SMC" | celltype == "Pericyte")



k=1:15
brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()

DimPlot(brain, label = T)

meta=brain@meta.data

meta$celltype <- NA 

meta$celltype[meta$seurat_clusters %in% c(4, 0)] <- "Pericyte_1"
meta$celltype[meta$seurat_clusters %in% c(1)] <- "Pericyte_2"
meta$celltype[meta$seurat_clusters %in% c(2)] <- "SMC_1"
meta$celltype[meta$seurat_clusters %in% c(3)] <- "SMC_2"

brain@meta.data = meta

Idents(brain) <- brain$celltype



results <- FindAllMarkers(brain, min.pct = 0.25, only.pos = T)
write.csv(results, file = "muralcell.cluster.markers.csv")



DimPlot(brain, label = T)

colors <- c("#006884", "#058799", "#E74C2E", "#EE9E80")

p <- DimPlot(brain, reduction = "umap", group.by = "celltype", cols = colors)

ggsave("MuralUMAP.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")


id_data='Mural.merged.clustered.data.rds'
saveRDS(brain,file=paste(id_data,'.rds',sep=''))



########################################################################################


#Figure 1G, H, Extended Data 1F, G


library(dplyr)
library(ggplot2)
library(grid)
library(readxl)
library(tidyr)
library(msigdbr)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(tibble)
library(fgsea)
library(limma)
library(edgeR)
library(matrixStats)
library(tidyverse)
library(Seurat)


brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data


brain$celltype.apoe <- paste(brain$celltype, brain$apoe, sep = "_")
Idents(brain) <- "celltype.apoe"
print(unique(Idents(brain)))


cell_types <- c("Pericyte_1", "Pericyte_2", "SMC_1", "SMC_2")



for (cell_type in cell_types) {
  
  ident_1 <- paste0(cell_type, "_e4")
  ident_2 <- paste0(cell_type, "_e3")
  
  
  result <- FindMarkers(object = brain, 
                        ident.1 = ident_1, 
                        ident.2 = ident_2,
                        test.use = "MAST",
                        latent.vars = c("age", "sex", "AD.status", "dataset"))
  
  
  filename <- paste0("MAST.apoe.", cell_type, ".age.sex.ADstatus.dataset.csv")
  
  
  write.csv(result, file = filename)
}




cell_types <- c("Pericyte_1", "Pericyte_2", "SMC_1", "SMC_2")


results_list <- list()


for (cell_type in cell_types) {
  
  filename <- paste0("MAST.apoe.", cell_type, ".age.sex.ADstatus.dataset.csv")
  
  
  df <- read.csv(filename, row.names = 1)
  
  
  df$cluster <- cell_type
  
  
  results_list[[cell_type]] <- df
}


merged_results <- do.call(rbind, results_list)


merged_results <- tibble::rownames_to_column(merged_results, var = "gene")


head(merged_results)

write.csv(merged_results, file = "merged.muralcell.MAST.results.csv")




process_data <- function(file) {
  read_excel(file, .name_repair = "unique") %>%        
    rename_with(~ "gene", 1) %>%                       
    mutate(gene = as.character(gene)) %>%              
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)   
}



Pericyte_1 <- process_data("MAST.apoe.Pericyte_1.age.sex.ADstatus.dataset.xlsx")
Pericyte_2 <- process_data("MAST.apoe.Pericyte_2.age.sex.ADstatus.dataset.xlsx")
SMC_1 <- process_data("MAST.apoe.SMC_1.age.sex.ADstatus.dataset.xlsx")
SMC_2 <- process_data("MAST.apoe.SMC_2.age.sex.ADstatus.dataset.xlsx")





count_genes <- function(df, cell_type) {
  df %>%
    mutate(direction = ifelse(avg_log2FC > 0, "Upregulated", "Downregulated")) %>%
    group_by(direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(Cell_Type = cell_type)
}


Pericyte_1_count <- count_genes(Pericyte_1, "Pericyte_1")
Pericyte_2_count <- count_genes(Pericyte_2, "Pericyte_2")
SMC_1_count <- count_genes(SMC_1, "SMC_1")
SMC_2_count <- count_genes(SMC_2, "SMC_2")


gene_counts <- bind_rows(Pericyte_1_count, Pericyte_2_count, SMC_1_count, SMC_2_count)


gene_counts_wide <- gene_counts %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = list(Upregulated = 0, Downregulated = 0)) %>%
  mutate(Downregulated = -Downregulated)  


p <- ggplot(gene_counts_wide, aes(x = reorder(Cell_Type, Upregulated))) +
  geom_bar(aes(y = Upregulated, fill = "Upregulated"), stat = "identity") +
  geom_bar(aes(y = Downregulated, fill = "Downregulated"), stat = "identity") +
  geom_text(aes(y = Upregulated + 10, label = abs(Upregulated)), color = "black", hjust = 0, size = 6, fontface = "bold") +  
  geom_text(aes(y = Downregulated - 10, label = abs(Downregulated)), color = "black", hjust = 1, size = 6, fontface = "bold") +  
  coord_flip() +  
  scale_fill_manual(values = c("Upregulated" = "#B22222", "Downregulated" = "#4682B4")) +
  labs(title = "APOE4 vs. APOE3 DEGs", x = NULL, y = NULL, fill = "Regulation") +  
  theme_minimal() +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.line = element_blank(),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),  
    axis.ticks.y = element_line(color = "black", size = 1.5),  
    axis.ticks.length = unit(0.4, "cm"),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  
  )


ggsave("DEGnumbers.svg", plot = p + 
         theme(
           axis.text.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 20)),  
           axis.ticks.length = unit(0.12, "cm"),  
           axis.ticks.y = element_line(size = 0.6),  
           plot.title = element_text(size = 11, face = "bold"),  
           legend.text = element_text(size = 7),  
           legend.title = element_text(size = 8),  
           legend.position = "bottom",  
           plot.margin = margin(10, 30, 10, 40)  
         ), 
       width = 4, height = 3, dpi = 300, device = "svg")





set.seed(42) 

cgp_gene_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
unique_gene_sets = split(x = cgp_gene_sets$gene_symbol, f = cgp_gene_sets$gs_name)



all_GSEA_results <- data.frame()



cell_types <- c("Pericyte_1", "Pericyte_2", "SMC_1", "SMC_2")



for (cell_type in cell_types) {
  
  filename <- paste0("MAST.apoe.", cell_type, ".age.sex.ADstatus.dataset.xlsx")
  
  
  mast <- read_excel(filename)
  
  
  rownames(mast) <- mast$...1
  
  
  transformed_value <- -log10(mast$p_val)
  transformed_value <- sign(mast$avg_log2FC) * abs(transformed_value)
  
  
  new_mast <- data.frame(transformed_value, row.names = rownames(mast))
  
  
  named_vector <- deframe(new_mast)
  named_vector <- setNames(named_vector, rownames(mast))
  named_vector <- sort(named_vector, decreasing = TRUE)
  
  
  GSEA_result <- fgseaMultilevel(unique_gene_sets, named_vector)
  
  
  GSEA_result <- GSEA_result %>%
    mutate(across(where(is.list), ~ sapply(., toString)))
  
  
  GSEA_result$cell_type <- cell_type
  
  
  all_GSEA_results <- rbind(all_GSEA_results, GSEA_result)
}

write.csv(all_GSEA_results, file = "GSEA_combined_results.csv", row.names = FALSE)



all_GSEA_results <- read.csv("GSEA_combined_results.csv")







filtered_data <- all_GSEA_results %>%
  filter(padj < 0.05)


filtered_data <- filtered_data %>%
  mutate(direction = ifelse(NES > 0, "Upregulated", "Downregulated"))


pathway_counts <- filtered_data %>%
  group_by(cell_type, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(Downregulated = -Downregulated)  




pathway_counts <- pathway_counts %>%
  mutate(Total = abs(Upregulated) + abs(Downregulated)) %>%
  arrange(desc(Total)) %>%  
  mutate(cell_type = factor(cell_type, levels = rev(unique(cell_type))))  

p <- ggplot(pathway_counts, aes(x = cell_type)) +
  geom_bar(aes(y = Upregulated, fill = "Upregulated"), stat = "identity") +
  geom_bar(aes(y = Downregulated, fill = "Downregulated"), stat = "identity") +
  geom_text(aes(y = Upregulated + 1, label = abs(Upregulated)), color = "black", hjust = 0, size = 6, fontface = "bold") +  
  geom_text(aes(y = Downregulated - 1, label = abs(Downregulated)), color = "black", hjust = 1, size = 6, fontface = "bold") +  
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#B22222", "Downregulated" = "#4682B4")) +
  labs(title = "Dysregulated Pathways per Cell Type", x = NULL, y = "Number of Dysregulated Pathways", fill = "Regulation") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.ticks.y = element_line(color = "black", size = 1),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)  
  )


ggsave("Dysregulated_Pathways.svg", plot = p + 
         theme(
           axis.text.y = element_text(size = 8, face = "bold", color = "black", margin = margin(r = 8)),  
           axis.ticks.length = unit(0.08, "cm"),  
           axis.ticks.y = element_line(size = 0.6), 
           plot.title = element_text(size = 10, face = "bold"), 
           legend.text = element_text(size = 6), 
           legend.title = element_text(size = 7),  
           legend.position = "bottom",  
           plot.margin = margin(5, 15, 5, 15) 
         ), 
       width = 4, height = 3, dpi = 300, device = "svg")


GSEA <- all_GSEA_results



significant_smc2 <- GSEA %>%
  filter(cell_type == "SMC_2", padj < 0.05) %>%
  mutate(
    pathway = gsub("^REACTOME_", "", pathway),  
    pathway = gsub("_", " ", pathway)           
  )

top_pathways <- significant_smc2 %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10) %>%   
  bind_rows(
    significant_smc2 %>%
      arrange(NES) %>%
      slice_head(n = 10)   
  )

top_pathways$pathway_factor <- factor(top_pathways$pathway, levels = top_pathways$pathway[order(top_pathways$NES)])

first_downregulated_y <- top_pathways %>%
  filter(NES < 0) %>%
  arrange(desc(NES)) %>%
  head(1) %>%
  pull(pathway)

first_downregulated_y_pos <- which(levels(top_pathways$pathway_factor) == first_downregulated_y)

p <- ggplot(top_pathways, aes(x = NES, y = pathway_factor, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("#4682B4", "#B22222")) +
  scale_x_continuous(limits = c(-2.5, 2.5)) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "solid") +
  geom_hline(yintercept = first_downregulated_y_pos + 0.5,  
             linetype = "dashed", color = "black", size = 0.4) +
  labs(
    x = "NES",
    y = NULL,
    title = "Top 10 Upregulated and Downregulated Pathways in SMC_2"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    plot.title = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "none"
  )

print(p)



selected_pathways_clean <- c(
  "REACTOME_MUSCLE_CONTRACTION",
  "REACTOME_RHO_GTPASE_EFFECTORS",
  "REACTOME_MAP2K_AND_MAPK_ACTIVATION",
  "REACTOME_ELASTIC_FIBRE_FORMATION",
  "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS",
  "REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND",
  "REACTOME_SCAVENGING_BY_CLASS_F_RECEPTORS"
)

pretty_labels <- c(
  "REACTOME_MUSCLE_CONTRACTION" = "Muscle contraction",
  "REACTOME_RHO_GTPASE_EFFECTORS" = "Rho GTPase signaling",
  "REACTOME_MAP2K_AND_MAPK_ACTIVATION" = "MAPK activation",
  "REACTOME_ELASTIC_FIBRE_FORMATION" = "Elastic fiber formation",
  "REACTOME_ECM_PROTEOGLYCANS" = "ECM proteoglycans",
  "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE" = "HSF1-mediated heat shock",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS" = "ETC and ATP synthesis",
  "REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND" = "HSP90 chaperone cycle",
  "REACTOME_SCAVENGING_BY_CLASS_F_RECEPTORS" = "Class F receptor scavenging"
)

top_pathways <- GSEA %>%
  filter(cell_type == "SMC_2", padj < 0.05) %>%
  filter(pathway %in% selected_pathways_clean) %>%
  mutate(pretty_pathway = pretty_labels[pathway]) %>%
  mutate(pretty_pathway = factor(pretty_pathway, levels = pretty_pathway[order(NES)]))

first_downregulated <- top_pathways %>%
  filter(NES < 0) %>%
  arrange(desc(NES)) %>%
  pull(pretty_pathway) %>%
  .[1]

first_downregulated_y_pos <- which(levels(top_pathways$pretty_pathway) == first_downregulated)

p <- ggplot(top_pathways, aes(x = NES, y = pretty_pathway, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("#4682B4", "#B22222")) +
  scale_x_continuous(limits = c(min(top_pathways$NES) - 0.5, max(top_pathways$NES) + 0.5)) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "solid") +
  geom_hline(yintercept = first_downregulated_y_pos + 0.5,
             linetype = "dashed", color = "black", size = 0.4) +
  labs(
    x = "NES",
    y = NULL,
    title = "Representative Enriched Pathways in SMC_2"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    plot.title = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "none"
  )

ggsave("Top_curated_Pathways_SMC2.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")




brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data



meta <- brain@meta.data
meta$cell_id <- rownames(meta)

cell_counts <- meta %>%
  group_by(subject, celltype) %>%
  tally(name = "n_cells") %>%
  filter(n_cells >= 5)

meta_filtered <- meta %>%
  inner_join(cell_counts, by = c("subject", "celltype"))

brain_filtered <- subset(brain, cells = meta_filtered$cell_id)

brain.pseudo <- AggregateExpression(brain_filtered, assays = "RNA", return.seurat = T, group.by = c("apoe", "subject", "celltype"))

tail(Cells(brain.pseudo))

library(dplyr)

meta_full <- brain_filtered@meta.data

subject_meta <- meta_full %>%
  dplyr::select(subject, age, sex, AD.status, dataset) %>%
  distinct()

pseudo_meta <- brain.pseudo@meta.data
pseudo_meta$sample <- rownames(pseudo_meta)


pseudo_meta <- left_join(pseudo_meta, subject_meta, by = "subject")

brain.pseudo@meta.data <- pseudo_meta

meta <- brain.pseudo@meta.data

meta_smc2 <- meta %>% filter(celltype == "SMC-2")

samples_to_keep <- meta_smc2$orig.ident


expr_mat <- brain.pseudo@assays$RNA$counts  

samples_to_keep <- intersect(samples_to_keep, colnames(expr_mat))

expr_mat<- expr_mat[, samples_to_keep, drop = FALSE]



dim(expr_mat)
head(colnames(expr_mat))

meta <- brain.pseudo@meta.data
meta$sample <- meta$orig.ident
head(meta)


samples_to_keep <- intersect(colnames(expr_mat), meta$sample)

meta_filtered <- meta %>%
  filter(sample %in% samples_to_keep)

meta_filtered <- meta_filtered[match(samples_to_keep, meta_filtered$sample), ]

stopifnot(all(meta_filtered$sample == colnames(expr_mat)))



deg_results <- list()
for (ct in unique(meta_filtered$celltype)) {
  message("Processing: ", ct)
  
  meta_ct <- meta_filtered %>% filter(celltype == ct)
  expr_ct <- expr_mat[, meta_ct$sample, drop = FALSE]
  
  meta_ct <- meta_ct[match(colnames(expr_ct), meta_ct$sample), ]
  stopifnot(ncol(expr_ct) == nrow(meta_ct))
  
  
  d <- DGEList(expr_ct)
  d <- calcNormFactors(d, method = 'upperquartile', p = 0.9)
  
  cpm_mat <- edgeR::cpm(d)
  k <- min(table(meta_ct$apoe))   
  keep <- rowSums(cpm_mat >= 1) >= 1   
  d <- d[keep,, keep.lib.sizes = FALSE]
  
  
  
  meta_ct$apoe <- factor(meta_ct$apoe, levels = c("e3", "e4"))
  meta_ct$sex <- factor(meta_ct$sex)
  meta_ct$dataset <- factor(meta_ct$dataset)
  meta_ct$AD.status <- factor(meta_ct$AD.status)
  design <- model.matrix(~ apoe + age + sex + dataset + AD.status, data = meta_ct)
  rownames(design) <- colnames(expr_ct)
  
  fit <-  voomLmFit(d, design = design, plot = T)
  fit <- eBayes(fit, trend = T)
  
  deg_table <- topTable(fit, coef = "apoee4", number = Inf, sort.by = "P")
  deg_table$gene <- rownames(deg_table)
  deg_results[[ct]] <- deg_table
}

results <- deg_results$`SMC-2`
results$gene <- rownames(results)

write.csv(test, file = "pseudobulk.DGE.SMC2.csv")



set.seed(10) 

cgp_gene_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
unique_gene_sets = split(x = cgp_gene_sets$gene_symbol, f = cgp_gene_sets$gs_name)


all_GSEA_results <- data.frame()



mast <- results


transformed_value <- -log10(mast$P.Value)
transformed_value <- sign(mast$logFC) * abs(transformed_value)

new_mast <- data.frame(transformed_value, row.names = rownames(mast))

named_vector <- deframe(new_mast)
named_vector <- setNames(named_vector, rownames(mast))
named_vector <- sort(named_vector, decreasing = TRUE)

GSEA_result <- fgseaMultilevel(unique_gene_sets, named_vector)

GSEA_result <- GSEA_result %>%
  mutate(across(where(is.list), ~ sapply(., toString)))


all_GSEA_results <- rbind(all_GSEA_results, GSEA_result)


write.csv(all_GSEA_results, file = "SMC2.pseudobulk.GSEA.csv")


all_GSEA_results <- read.csv('SMC2.pseudobulk.GSEA.csv')


filtered_data <- all_GSEA_results %>%
  filter(padj < 0.05)




GSEA <- all_GSEA_results



selected_pathways_clean <- c(
  "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_SMOOTH_MUSCLE_CONTRACTION",
  "REACTOME_SYNDECAN_INTERACTIONS",
  "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX"
)

pretty_labels <- c(
  "REACTOME_ECM_PROTEOGLYCANS" = "ECM proteoglycans",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION" = "ECM organization",
  "REACTOME_SMOOTH_MUSCLE_CONTRACTION" = "Smooth muscle contraction",
  "REACTOME_SYNDECAN_INTERACTIONS" = "Syndecan interactions",
  "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX" = "TGF-β signaling"
)

top_pathways <- GSEA %>%
  filter(padj < 0.05) %>%
  filter(pathway %in% selected_pathways_clean) %>%
  mutate(pretty_pathway = pretty_labels[pathway]) %>%
  mutate(pretty_pathway = factor(pretty_pathway, levels = pretty_pathway[order(NES)]))

p <- ggplot(top_pathways, aes(x = NES, y = pretty_pathway, fill = NES > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("#B22222", "#B22222")) +
  labs(
    x = "NES",
    y = NULL,
    title = "Upregulated APOE4 SMC 2 pathways via limma pseudobulk analysis"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    plot.title = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "none"
  )

print(p)


ggsave("E4SMC2pseudobulkGSEA.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



#######################################################################################

#Figure 1I, Extended Data 1H-J


brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data


brain <- AddModuleScore(brain, 
                        features = list(c("ACTA2", "CDH11", "ITGA1", "TPM2", "TAGLN", "PALLD", "SORBS1", "ITGA2", "ITGA8", "LMOD1"
                        )),
                        name = "contraction.genes",
                        ctrl = 20)
meta = brain@meta.data



brain <- AddModuleScore(brain,
                        features = list(c("FN1", "COL4A1", "COL6A2", "COL1A1", "FBN1", "COL3A1", "COL1A2", "DCN", "COL14A1", "COL5A2")),
                        name = "ecm.genes",
                        ctrl = 20)

meta = brain@meta.data




brain$coexpression.genes1 <- pmin(brain$contraction.genes1, brain$ecm.genes1)






p <- FeaturePlot(
  brain, 
  features = "JointMinScore", 
  reduction = "umap",
  min.cutoff = 'q90',
  max.cutoff = 1,
  order = T
)

ggsave("merge_scale.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



p <- FeaturePlot(
  brain, 
  features = "JointMinScore", 
  reduction = "umap", 
  split.by = "apoe",
  min.cutoff = 'q90',
  max.cutoff = 1,
  order = T
)

ggsave("merge.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



p <- FeaturePlot(
  brain, 
  features = "contraction.genes1", 
  reduction = "umap", 
  min.cutoff = 'q90',
  max.cutoff = 1.75,
  order = T
)

ggsave("contractionscale.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")


p <- FeaturePlot(
  brain, 
  features = "contraction.genes1", 
  reduction = "umap", 
  split.by = "apoe",
  min.cutoff = 'q90',
  max.cutoff = 1.75,
  order = T
)

ggsave("contraction.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")


p <- FeaturePlot(
  brain, 
  features = "ecm.genes1", 
  reduction = "umap", 
  min.cutoff = 'q90',
  max.cutoff = 1,
  order = T
)

ggsave("ecmscale.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")


p <- FeaturePlot(
  brain, 
  features = "ecm.genes1", 
  reduction = "umap", 
  split.by = "apoe",
  min.cutoff = 'q90',
  max.cutoff = 1,
  order = T
)

ggsave("ecm.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



meta = brain@meta.data

brain <- subset(brain, subset = celltype == "SMC_2")

meta = brain@meta.data

library(dplyr)
library(ggplot2)

min_cells <- 5

cell_counts <- brain@meta.data %>%
  dplyr::count(subject) %>%
  dplyr::rename_with(~ "n_cells", .cols = "n")


keep_subjects <- cell_counts %>%
  filter(n_cells >= min_cells) %>%
  pull(subject)

filtered_meta <- brain@meta.data %>%
  filter(subject %in% keep_subjects)

pseudobulk_df <- filtered_meta %>%
  group_by(subject, apoe, age, sex, dataset, AD.status) %>%
  summarise(mean_contraction_score = mean(contraction.genes1, na.rm = TRUE), .groups = "drop")



min_score <- min(pseudobulk_df$mean_contraction_score, na.rm = TRUE)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_contraction_score_shifted =mean_contraction_score - min_score)


mean_e3 <- pseudobulk_df %>%
  filter(apoe == "e3") %>%
  summarise(mean_val = mean(mean_contraction_score_shifted, na.rm = TRUE)) %>%
  pull(mean_val)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_contraction_score_normalized = mean_contraction_score_shifted / mean_e3)


mean_e4 <- pseudobulk_df %>%
  filter(apoe == "e4") %>%
  summarise(mean_val = mean(mean_contraction_score_normalized, na.rm = TRUE)) %>%
  pull(mean_val)

print(mean_e4)

write.csv(pseudobulk_df, file = "pseudobulk.contraction.csv", row.names = FALSE)







min_cells <- 5

cell_counts <- brain@meta.data %>%
  dplyr::count(subject) %>%
  dplyr::rename_with(~ "n_cells", .cols = "n")

keep_subjects <- cell_counts %>%
  filter(n_cells >= min_cells) %>%
  pull(subject)

filtered_meta <- brain@meta.data %>%
  filter(subject %in% keep_subjects)

pseudobulk_df <- filtered_meta %>%
  group_by(subject, apoe, age, sex, dataset, AD.status) %>%
  summarise(mean_ecm_score = mean(ecm.genes1, na.rm = TRUE), .groups = "drop")


min_score <- min(pseudobulk_df$mean_ecm_score, na.rm = TRUE)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_ecm_score_shifted = mean_ecm_score - min_score)


mean_e3 <- pseudobulk_df %>%
  filter(apoe == "e3") %>%
  summarise(mean_val = mean(mean_ecm_score_shifted, na.rm = TRUE)) %>%
  pull(mean_val)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_ecm_score_normalized = mean_ecm_score_shifted / mean_e3)


mean_e4 <- pseudobulk_df %>%
  filter(apoe == "e4") %>%
  summarise(mean_val = mean(mean_ecm_score_normalized, na.rm = TRUE)) %>%
  pull(mean_val)

print(mean_e4)

write.csv(pseudobulk_df, file = "pseudobulk.ecm.csv", row.names = FALSE)






min_cells <- 5

cell_counts <- brain@meta.data %>%
  dplyr::count(subject) %>%
  dplyr::rename_with(~ "n_cells", .cols = "n")


keep_subjects <- cell_counts %>%
  filter(n_cells >= min_cells) %>%
  pull(subject)

filtered_meta <- brain@meta.data %>%
  filter(subject %in% keep_subjects)

pseudobulk_df <- filtered_meta %>%
  group_by(subject, apoe, age, sex, AD.status, dataset) %>%
  summarise(mean_coexpression_score = mean(coexpression.genes1, na.rm = TRUE), .groups = "drop")



min_score <- min(pseudobulk_df$mean_coexpression_score, na.rm = TRUE)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_coexpression_score_shifted = mean_coexpression_score - min_score)


mean_e3 <- pseudobulk_df %>%
  filter(apoe == "e3") %>%
  summarise(mean_val = mean(mean_coexpression_score_shifted, na.rm = TRUE)) %>%
  pull(mean_val)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_coexpression_score_normalized = mean_coexpression_score_shifted / mean_e3)


mean_e4 <- pseudobulk_df %>%
  filter(apoe == "e4") %>%
  summarise(mean_val = mean(mean_coexpression_score_normalized, na.rm = TRUE)) %>%
  pull(mean_val)

print(mean_e4)


write.csv(pseudobulk_df, file = "pseudobulk.coexpression.csv", row.names = FALSE)







results <- read.csv('pseudobulk.DGE.SMC2.csv')

up.sig <- subset(results, subset = logFC > 0.25 & P.Value < 0.05)

gene.list.up <- up.sig$X


library(enrichR)


dbs <- c("CellMarker_Augmented_2021")


enrich_results <- enrichr(gene.list.up, dbs)[[1]]

parts <- do.call(rbind, strsplit(enrich_results$Overlap, "/"))
enrich_results$k <- as.numeric(parts[,1])
enrich_results$M <- as.numeric(parts[,2])

res_f <- subset(enrich_results, M >= 10 & k >= 3) 


CellMarker <- res_f


write.csv(CellMarker, file = "SMC2.pseudobulk.enrichr.csv")

library(ggplot2)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(stringr)

plot_df <- CellMarker %>%
  mutate(
    CellType = str_extract(Term, "^[^:]+"),
    CellClass = case_when(
      str_detect(CellType, "Myofibroblast") ~ "Myofibroblast",
      str_detect(CellType, "Smooth Muscle|Vascular Smooth Muscle") ~ "SMC",
      str_detect(CellType, "Pericyte") ~ "Pericyte",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(CellType) %>%
  slice_max(order_by = Combined.Score, n = 1) %>%
  ungroup() %>%
  arrange(desc(Combined.Score)) %>%
  slice_head(n = 5) %>%
  mutate(
    Combined_Score_Scaled = Combined.Score / 1000,
    Term_clean = str_wrap(CellType, width = 30),
    Term_clean = factor(Term_clean, levels = rev(Term_clean[order(Combined.Score, decreasing = TRUE)]))
  )

colors <- c("Myofibroblast" = "#B22222", "SMC" = "gray60", "Pericyte" = "gray60", "Other" = "gray60")

p <- ggplot(plot_df, aes(x = Combined_Score_Scaled, y = Term_clean, fill = CellClass)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Combined Score (×10³), Adjusted P < 0.0001",
    y = NULL,
    title = "Top Cell Type Annotations of APOE4 SMC2 Cluster"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  )

print(p)

ggsave("E4SMC2_cellmarker_enrich.svg", plot = p, width = 3, height = 3, dpi = 300, device = "svg")








brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data

brain <- subset(brain, subset = celltype == "SMC_2")



meta <- brain@meta.data


cell_counts <- meta %>%
  group_by(subject) %>%
  tally(name = "n_cells") %>%
  filter(n_cells >= 5)


brain_filtered <- subset(brain, subset = subject %in% cell_counts$subject)



brain.pseudo <- AggregateExpression(brain_filtered, assays = "RNA", return.seurat = T, group.by = c("apoe", "subject"))

tail(Cells(brain.pseudo))


original_meta <- brain@meta.data


pseudo_samples <- colnames(brain.pseudo)


pseudo_meta <- data.frame(sample = pseudo_samples) %>%
  tidyr::separate(sample, into = c("apoe", "subject"), sep = "_") %>%
  mutate(subject = as.character(subject))  


per_subject_meta <- original_meta %>%
  dplyr::select(subject, age, sex, dataset, AD.status) %>%
  distinct()  


pseudo_meta <- pseudo_meta %>%
  left_join(per_subject_meta, by = "subject") %>%
  column_to_rownames("subject")  


brain.pseudo <- AddMetaData(brain.pseudo, metadata = pseudo_meta)


library(Seurat)
library(GSVA)
library(tidyverse)


expr_mat <- brain.pseudo@assays$RNA$data  


dim(expr_mat)
head(colnames(expr_mat))  


library(readxl)


markers <- read_excel('MyofibDEGs.xlsx')


clusters_to_include <- c("c11")


all_top_expressed_genes <- character()


for (cluster_id in clusters_to_include) {
  
  cluster_markers <- subset(markers, cluster == cluster_id)
  
  
  cluster_genes <- cluster_markers$gene
  expressed_genes <- intersect(cluster_genes, rownames(expr_mat))
  
  
  mean_expr <- rowMeans(expr_mat[expressed_genes, , drop = FALSE])
  
  
  top_genes <- names(sort(mean_expr, decreasing = TRUE))[1:min(10, length(mean_expr))]
  
  
  all_top_expressed_genes <- c(all_top_expressed_genes, top_genes)
}


union_top_expressed_genes <- unique(all_top_expressed_genes)


gene_set_list <- list(myofibroblast_signature = union_top_expressed_genes)



library(GSVA)


test <- as.matrix(expr_mat)

param <- gsvaParam(test, gene_set_list)

gsva_scores <- gsva(param, verbose = T)



gsva_df <- as.data.frame(t(gsva_scores))  
gsva_df$subject <- rownames(gsva_df)


meta_df <- brain.pseudo@meta.data




meta_df$subject <- rownames(meta_df)



if (!"subject" %in% colnames(meta_df)) {
  meta_df <- meta_df %>%
    tibble::rownames_to_column("subject")
}

meta_df <- meta_df %>%
  dplyr::select(subject, apoe, age, sex, AD.status, dataset)



plot_df <- left_join(gsva_df, meta_df, by = "subject")



min_val <- min(plot_df$myofibroblast_signature, na.rm = TRUE)
plot_df$gsva_shifted <- plot_df$myofibroblast_signature - min_val


e3_mean <- mean(plot_df$gsva_shifted[plot_df$apoe == "e3"], na.rm = TRUE)
plot_df$gsva_norm <- plot_df$gsva_shifted / e3_mean



write.csv(plot_df, file = "pseudobulk.gsva.myofibroblast.csv")


######################################################################################

#Figure 1J, Extended Data 2A-D




library(DAseq)


brain = readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data


labels_res <- unique(meta[meta$apoe == "e3", "subject"])


labels_nonres <- unique(meta[meta$apoe == "e4", "subject"])






labels_res <- as.character(labels_res)
labels_nonres <- as.character(labels_nonres)

plot.embedding <- as.data.frame(Embeddings(brain, reduction = "umap"))

X_pca <- as.matrix(Embeddings(brain, reduction = "harmony"))

cell.labels <- as.character(meta$subject)  



da_cells <- getDAcells(
  X = X_pca,
  cell.labels = cell.labels,
  labels.1 = labels_res,
  labels.2 = labels_nonres,
  k.vector = seq(50, 500, 50),
  plot.embedding = plot.embedding
)


da_cells$pred.plot


da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(-0.7,0.7),
  plot.embedding = plot.embedding
)

da_cells$da.cells.plot







meta$cell_name <- rownames(meta)


upregulated_cells <- rownames(meta)[da_cells$da.up]
downregulated_cells <- rownames(meta)[da_cells$da.down]


meta$regulation_status <- "Neutral"


meta$regulation_status[meta$cell_name %in% upregulated_cells] <- "Upregulated"


meta$regulation_status[meta$cell_name %in% downregulated_cells] <- "Downregulated"


table(meta$regulation_status)

brain@meta.data = meta

library(ggplot2)


umap_coords <- as.data.frame(Embeddings(brain, reduction = "umap"))
umap_coords$regulation_status <- brain$regulation_status


colors <- c("Neutral" = "lightgray", "Upregulated" = "red", "Downregulated" = "blue")


p <- ggplot() +
  
  geom_point(data = umap_coords[umap_coords$regulation_status == "Neutral", ], 
             aes(x = umap_1, y = umap_2, color = regulation_status), 
             size = 0.5, alpha = 0.5) +  
  
  geom_point(data = umap_coords[umap_coords$regulation_status == "Upregulated", ], 
             aes(x = umap_1, y = umap_2, color = regulation_status), 
             size = 0.5) +  
  
  geom_point(data = umap_coords[umap_coords$regulation_status == "Downregulated", ], 
             aes(x = umap_1, y = umap_2, color = regulation_status), 
             size = 0.5) +  
  ggtitle("APOE4 DA Regions") +  
  scale_color_manual(values = colors, name = "Regulation Status") +  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    axis.line = element_line(color = "black"),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  
  )


ggsave("regulation_status_plot.svg", plot = p, width = 4, height = 3, dpi = 300)




counts <- brain[["RNA"]]$counts
m <- brain@meta.data

save(list = setdiff(ls(), c("brain", "meta", "k")), file = "mural.DAseq.RData")



#RESTART R HERE
unloadNamespace("DAseq")
library(Seurat, lib.loc="C:/seurat_v4_lib")
library(DAseq)

load("mural.DAseq.RData")




da_regions <- getDAregion(
  X = X_pca,
  da.cells = da_cells,
  cell.labels = cell.labels,
  labels.1 = labels_res,
  labels.2 = labels_nonres,
  resolution = 0.05,
  plot.embedding = plot.embedding
)


da_regions$da.region.plot




brain.old <- CreateSeuratObject(counts = counts, meta.data = m)

brain.old <- addDAslot(brain.old, da_regions, da.slot = "da", set.ident = F)


saveRDS(brain.old, file = "OLDSEURAT.mural.DAseq.rds")



#RESTART R HERE
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


brain.old <- readRDS("OLDSEURAT.mural.DAseq.rds")
meta.old = brain.old@meta.data

brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data

meta$da = meta.old$da


meta$celltype[meta.old$da == 1] <- "myoVC"

brain@meta.data = meta

Idents(brain) <- brain$celltype

DimPlot(brain, label = T)

saveRDS(brain, file = "myoVC.mural.clustered.data.rds")

results <- FindMarkers(brain, ident.1 = "myoVC", ident.2 = "SMC_2", test.use = "MAST")

results$gene <- rownames(results)

write.csv(results, file = "myoVCSMC2.mast.csv")


library(readxl)

results <- read_excel("myoVCSMC2.mast.xlsx")

up.sig <- subset(results, subset = avg_log2FC > 0.25 & p_val_adj < 0.05)

down.sig <- subset(results, subset = avg_log2FC < -0.25 & p_val_adj < 0.05)


gene.list.up <- up.sig$gene


library(enrichR)





dbs <- c("CellMarker_Augmented_2021")


enrich_results <- enrichr(gene.list.up, dbs)[[1]]

parts <- do.call(rbind, strsplit(enrich_results$Overlap, "/"))
enrich_results$k <- as.numeric(parts[,1])
enrich_results$M <- as.numeric(parts[,2])

res_f <- subset(enrich_results, M >= 10 & k >= 3) 


CellMarker <- res_f


write.csv(CellMarker, file = "DAseqE4region.enrichr.csv")




library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(readxl)


filtered_DEGs_UP <- subset(results, subset = avg_log2FC > 0.25 & p_val_adj < 0.05)
gene_list_up <- filtered_DEGs_UP$...1

cP_GO_UP <- enrichGO(
  gene = gene_list_up,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  readable = F,
  pool = T
)


cP_GO_UP <- clusterProfiler::simplify(cP_GO_UP)


filtered_DEGs_DOWN <- subset(results, subset = avg_log2FC < -0.25 & p_val_adj < 0.05)
gene_list_down <- filtered_DEGs_DOWN$...1

cP_GO_DOWN <- enrichGO(
  gene = gene_list_down,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  readable = F,
  pool = T
)


cP_GO_DOWN <- clusterProfiler::simplify(cP_GO_DOWN)


write.csv(results.all, file = "DAseqE4region.CP.csv")



library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(grid)


myoVC_df <- as.data.frame(cP_GO_UP) %>%
  arrange(p.adjust) %>%
  slice_head(n = 3)
myoVC_df$Cluster <- "MyoVC"


SMC2_df <- as.data.frame(cP_GO_DOWN) %>%
  arrange(p.adjust) %>%
  slice_head(n = 3)
SMC2_df$Cluster <- "SMC 2"


myoVC_df$signed_log10_p <- -log10(myoVC_df$p.adjust)  
SMC2_df$signed_log10_p <- log10(SMC2_df$p.adjust)     


combined_df <- bind_rows(myoVC_df, SMC2_df)


combined_df <- combined_df %>%
  arrange(signed_log10_p) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


p <- ggplot(combined_df, aes(x = signed_log10_p, y = Description, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("MyoVC" = "#B22222", "SMC 2" = "#4682B4")) +
  scale_x_continuous(limits = c(-10, 30), expand = c(0, 0)) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  labs(
    x = "Signed -Log10 Adjusted P-value",
    y = NULL,
    title = "Top 3 Enriched Pathways in MyoVC and SMC 2"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    plot.title = element_text(size = 10, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom"
  )


ggsave("MyoVC_vs_SMC2_Pathways.svg", plot = p, width = 8, height = 8, dpi = 300, device = "svg")




library(ggplot2)
library(dplyr)
library(grid)


myoVC_genes <- c("ACTA2", "FN1", "COL3A1", "COL8A1", "TGFB1I1")
smc2_genes  <- c("FLT1", "CLDN5", "TIMP3", "ADAMTS9", "EPAS1")
gene_list <- c(myoVC_genes, smc2_genes)


plot_df <- results %>%
  filter(gene %in% gene_list, p_val_adj < 0.05) %>%
  mutate(
    Direction = case_when(
      gene %in% myoVC_genes ~ "MyoVC",
      gene %in% smc2_genes ~ "SMC 2"
    )
  )




plot_df <- plot_df %>%
  arrange(avg_log2FC) %>%
  mutate(gene = factor(gene, levels = gene))  


gene_levels <- levels(plot_df$gene)
split_index <- sum(plot_df$avg_log2FC < 0)
y_sep <- split_index + 0.5  





p <- ggplot(plot_df, aes(x = avg_log2FC, y = gene, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("MyoVC" = "#B22222", "SMC 2" = "#4682B4")) +
  geom_vline(xintercept = 0, color = "black", size = 0.6) +
  annotate("segment", x = 0.1, xend = 0.8, y = 10.8, yend = 10.8,
           arrow = arrow(length = unit(0.18, "inches")), color = "#B22222", size = 0.8) +
  annotate("text", x = 0.45, y = 11.2, label = "MyoVC Enriched", size = 4.2, fontface = "bold", color = "#B22222") +
  annotate("segment", x = -0.1, xend = -0.8, y = 0.2, yend = 0.2,
           arrow = arrow(length = unit(0.18, "inches")), color = "#4682B4", size = 0.8) +
  annotate("text", x = -0.45, y = -0.3, label = "SMC 2 Enriched", size = 4.2, fontface = "bold", color = "#4682B4") +
  labs(
    x = "log2(FC) (FDR < 0.05)",
    y = NULL,
    subtitle = "Expressed in >10% of cells"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.ticks = element_line(color = "black", size = 0.8),  
    axis.ticks.length = unit(0.15, "cm"),
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),  
    axis.text.x = element_text(face = "bold", size = 9, color = "black"),
    axis.title.x = element_text(face = "bold", size = 10, color = "black"),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    legend.position = "none"
  ) +
  geom_hline(yintercept = y_sep, linetype = "dashed", color = "black", size = 1)



ggsave("MyoVC_SMC2_log2FC_final.svg", plot = p, width = 3, height = 3, dpi = 300, device = "svg")






library(ggplot2)
library(dplyr)
library(stringr)


plot_df <- CellMarker %>%
  mutate(
    CellType = str_extract(Term, "^[^:]+"),
    CellClass = case_when(
      str_detect(CellType, "Myofibroblast") ~ "Myofibroblast",
      str_detect(CellType, "Smooth Muscle|Vascular Smooth Muscle") ~ "SMC",
      str_detect(CellType, "Pericyte") ~ "Pericyte",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(CellType) %>%
  slice_max(order_by = Combined.Score, n = 1) %>%
  ungroup() %>%
  group_by(CellClass) %>%
  slice_max(order_by = Combined.Score, n = 1) %>%
  ungroup() %>%
  mutate(
    Combined_Score_Scaled = Combined.Score / 1000,
    Term_clean = str_wrap(CellType, width = 30),
    Term_clean = factor(Term_clean, levels = rev(Term_clean[order(Combined.Score, decreasing = TRUE)]))
  )


colors <- c("Myofibroblast" = "#B22222", "SMC" = "gray60", "Pericyte" = "gray60", "Other" = "gray60")


p <- ggplot(plot_df, aes(x = Combined_Score_Scaled, y = Term_clean, fill = CellClass)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(
    x = "Combined Score (×10³)",
    y = NULL,
    title = "Top Cell Types Matching MyoVC Signature (One per Class)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm")
  )

ggsave("MyoVC_cellmarker_enrich.svg", plot = p, width = 3, height = 3, dpi = 300, device = "svg")






library(ggrepel)
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)
library(destiny)
library(slingshot); library(SingleCellExperiment)
library(RColorBrewer); library(scales)
library(viridis); library(UpSetR)
library(pheatmap); library(msigdbr)
library(fgsea); library(knitr)
library(ggplot2); library(gridExtra)
library(tradeSeq); library(Seurat)
library(plot3D); library(fields); library(rgl)
library(DAseq); library(SeuratExtend)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(ggrepel)
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)


set.seed(1234)

brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data




brain_sce <- as.SingleCellExperiment(brain)
brain_milo <- Milo(brain_sce)


brain_milo <- buildGraph(brain_milo, k = 30, d = 30, reduced.dim = "HARMONY")
brain_milo <- makeNhoods(brain_milo, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "HARMONY")
plotNhoodSizeHist(brain_milo)
brain_milo <- countCells(brain_milo, meta.data = data.frame(colData(brain_milo)), sample="subject")




brain_milo <- calcNhoodDistance(brain_milo, d=30, reduced.dim = "HARMONY")

saveRDS(brain_milo, file = "mural.miloR.rds")



brain_design <- data.frame(colData(brain_milo))[,c("subject", "age", "sex", "AD.status", "apoe", "dataset")]

brain_design <- distinct(brain_design)
rownames(brain_design) <- brain_design$subject



da_results <- testNhoods(brain_milo, reduced.dim = "HARMONY", design = ~ age + sex + AD.status + dataset + apoe, design.df = brain_design)



library(Matrix)
nc <- nhoodCounts(brain_milo)  


sum_row0 <- sum(Matrix::rowSums(nc) == 0)  
sum_col0 <- sum(Matrix::colSums(nc) == 0)  
which_col0 <- which(Matrix::colSums(nc) == 0)
colnames(nc)[which_col0]


missing_in_design <- setdiff(colnames(nc), rownames(brain_design))
missing_in_counts <- setdiff(rownames(brain_design), colnames(nc))


keep_samples <- colnames(nc)[Matrix::colSums(nc) > 0]
nc <- nc[, keep_samples, drop = FALSE]
brain_design <- droplevels(brain_design[keep_samples, , drop = FALSE])


keep_nhoods <- Matrix::rowSums(nc) > 0
nc <- nc[keep_nhoods, , drop = FALSE]




da_results <- testNhoods(
  brain_milo,
  reduced.dim = "HARMONY",
  design = ~ age + sex + AD.status + dataset + apoe,
  design.df = brain_design
)







ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) 

da_results <- annotateNhoods(brain_milo, da_results, coldata_col = "celltype")


plotDAbeeswarm(da_results, group.by = "celltype")

sig <- subset(da_results, subset = SpatialFDR < 0.05)


brain_milo <- buildNhoodGraph(brain_milo)

test <- brain_milo


library(igraph)
library(patchwork)


significance_threshold <- 0.05


significant_nhoods <- which(da_results$SpatialFDR < significance_threshold)


da_results_filtered <- da_results[significant_nhoods, ]


graph <- brain_milo@nhoodGraph[[1]]  
filtered_graph <- induced_subgraph(graph, vids = significant_nhoods)
test@nhoodGraph[[1]] <- filtered_graph


umap_pl <- plotReducedDim(test, dimred = "UMAP", colour_by="celltype", text_by = "celltype", text_size = 3) +
  guides(fill="none")


nh_graph_pl <- plotNhoodGraphDA(test, da_results_filtered, layout="UMAP", alpha=0.3) +
  scale_size_continuous(range = c(1, 5))


p <- nh_graph_pl

ggsave("Mural.neighborhoods.graph.svg", plot = p, width = 4, height = 4, dpi = 300, device = "svg")












########################################################################################

#Figure 1K, Extended Data 2E


brain <- readRDS("OLDSEURAT.mural.DAseq.rds")


df <- brain@meta.data

Idents(brain) <- df$celltype

df$cluster <- Idents(brain)


cluster_of_interest <- c("SMC_2")
df_cluster <- df[df$cluster == cluster_of_interest, ]

df_cluster <- subset(df_cluster, subset = AD.status == "nonAD")


df_cluster$high_both <- df_cluster$da == 1





library(dplyr)

subject_level <- df_cluster %>%
  group_by(subject, apoe) %>%
  summarise(
    total_cells = n(),
    high_both_cells = sum(high_both),
    prop_high_both = high_both_cells / total_cells,
    age = dplyr::first(age),
    sex = dplyr::first(sex),
    ad = dplyr::first(AD.status),
    dataset = dplyr::first(dataset),
    region = dplyr::first(brain_region),
    .groups = "drop"
  )





subject_level$apoe <- factor(subject_level$apoe, levels = c("e3", "e4"))
subject_level$sex <- as.factor(subject_level$sex)
subject_level$dataset <- as.factor(subject_level$dataset)




library(splines)
model_random <- glm(
  cbind(high_both_cells, total_cells - high_both_cells) ~ apoe + ns(age, df = 3) + sex + dataset,
  family = binomial(),
  data = subject_level
)



summary(model_random)



model_covariate <- model_random


library(car)

vif(model_covariate)


summary(model_covariate)$deviance
summary(model_covariate)$null.deviance

1 - (model_covariate$deviance / model_covariate$null.deviance)

library(DHARMa)
p <- simulateResiduals(model_covariate, plot = TRUE)



library(broom.mixed)
library(broom)
library(ggplot2)
library(dplyr)
library(forcats)
library(grid)  

coef_df_random <- tidy(model_random, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)", term != "sd__(Intercept)") %>%
  mutate(
    term = dplyr::recode(term,
                         "apoee4" = "APOE4 allele",
                         "datasetyang" = "Dataset: Yang (vs. Sun)",
                         "sexM" = "Sex: Male (vs. Female)",
                         "ns(age, df = 3)1" = "Age (nonlinear)",
                         "ns(age, df = 3)2" = NA_character_,  
                         "ns(age, df = 3)3" = NA_character_
    )
  ) %>%
  filter(!is.na(term)) %>%  
  mutate(
    term = fct_reorder(term, estimate),
    significant = ifelse(conf.low > 1 | conf.high < 1, "Yes", "No")
  )


p_random <- ggplot(coef_df_random, aes(x = term, y = estimate, color = significant)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15, linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_y_log10() +
  annotation_logticks(
    sides = "b",
    short = unit(1.5, "mm"),
    mid = unit(2.5, "mm"),
    long = unit(3.5, "mm")
  ) +
  scale_color_manual(values = c("Yes" = "#D55E00", "No" = "gray40")) +
  labs(
    x = NULL,
    y = "Odds Ratio (log scale)",
    title = "Predictors of high contraction and ECM gene co-expression in SMC_2 cells (Random Effects Model)",
    color = "Significant (95% CI)"
  ) +
  coord_flip() +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.y = element_text(hjust = 0, size = 12),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 20, 10, 10),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )


print(p_random)


ggsave("forestplot_with_interaction.svg", plot = p_random, width = 6, height = 4.5, dpi = 300, device = "svg")



subject_level$pred <- predict(model_covariate, type = "response")

write.csv(subject_level, file = "myopred.csv")


#########################################################################################

#Extended Data 2F




df <- read_excel('myoVCSMC2.mast.xlsx')



filtered_df <- df[df$p_val_adj < 0.05, ]
filtered_df <- filtered_df[filtered_df$avg_log2FC > 0.25, ]

rownames(filtered_df) <- filtered_df$...1


ordered_df <- filtered_df[order(filtered_df$p_val_adj), ]

rownames(ordered_df) <- ordered_df$...1



gene_names <- rownames(ordered_df)


head(gene_names)


library(babelgene)


results <- orthologs(gene_names, species = "mouse", human = TRUE)

head(results)

ms.genes <- results$ensembl





ms <- read.table("4mo_vsd_wooutliers5.txt", header = TRUE)
rownames(ms) <- ms$GeneID
meta.ms <- read.table("hapoe34run_designfile_4mo_wPenLitter.txt", header = TRUE)


meta.ms$Sample_Name <- paste0("X", meta.ms$Sample_Name)


ms <- ms[, colnames(ms) %in% meta.ms$Sample_Name]


meta.ms <- meta.ms[match(colnames(ms), meta.ms$Sample_Name), ]


meta.ms$new_name <- ave(meta.ms$TERM_3, meta.ms$TERM_3, FUN = function(x) {
  paste0(x, "_", seq_along(x))
})


colnames(ms) <- meta.ms$new_name




gene.set <- list(ms.genes)


library(GSVA)

test <- as.matrix(ms)

param <- gsvaParam(test, gene.set)

gsva.results <- gsva(param, verbose = T)





gsva_t <- t(gsva.results)



apoe_dose <- as.numeric(sub("_.*", "", rownames(gsva_t)))  


apoe4_dose <- ifelse(apoe_dose == 33, 0, ifelse(apoe_dose == 34, 1, 2))


gsva_scores <- gsva_t[, 1]



cor.test(apoe4_dose, gsva_scores, method = "pearson")



library(ggplot2)

df <- data.frame(
  APOE4_dose = apoe4_dose,
  GSVA_score = gsva_scores
)

ggplot(df, aes(x = APOE4_dose, y = GSVA_score)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(title = "Correlation of GSVA score with APOE4 dosage",
       x = "APOE4 alleles (0, 1, 2)",
       y = "GSVA score")






df <- data.frame(
  GSVA = gsva_scores,
  APOE4_dose = apoe4_dose,
  Sed = meta.ms$TERM_4,
  Sex = meta.ms$TERM_1  
)


model <- lm(GSVA ~ APOE4_dose + Sex + Sed, data = df)


summary(model)




library(ggplot2)  






gsva_matrix <- as.data.frame(gsva.results)
metadata_df <- meta.ms



library(ggplot2)
library(ggridges)


gsva_df <- data.frame(
  new_name = colnames(gsva_matrix),
  GSVA_score  = as.numeric(gsva_matrix[1, ])  
)

plot_df <- merge(gsva_df, metadata_df, by = "new_name")


plot_df$TERM_3 <- factor(plot_df$TERM_3, levels = c("33", "34", "44"))

plot_df$APOE_genotype <- plot_df$TERM_3   
plot_df$Sex <- plot_df$TERM_1             





ggplot(plot_df, aes(x = GSVA_score, y = APOE_genotype)) +
  
  geom_density_ridges(
    fill = "gray70",        
    alpha = 0.5,            
    color = "black",        
    scale = 0.8             
  ) +
  
  geom_boxplot(
    width = 0.15,           
    color = "black",        
    fill = NA,              
    size = 0.8,             
    outlier.shape = NA      
  ) +
  
  geom_jitter(
    aes(color = Sex),       
    size = 2,               
    height = 0.1, width = 0 
  ) +
  
  scale_color_manual(
    name = "Sex", 
    values = c("F" = "#E64B35", "M" = "#4DBBD5") 
  ) +
  
  scale_y_discrete(limits = c("33", "34", "44")) +
  
  labs(x = "GSVA score", y = "APOE Genotype") +
  
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),        
    axis.text  = element_text(size = 12),        
    legend.position = "right",                   
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )






library(ggplot2)
library(ggridges)


plot_df$APOE_genotype <- factor(plot_df$APOE_genotype, levels = c("33", "34", "44"))


label_map <- c("33" = "APOE3/3", "34" = "APOE3/4", "44" = "APOE4/4")


pval <- 0.038


p <- ggplot(plot_df, aes(x = GSVA_score, y = APOE_genotype)) +
  
  
  geom_density_ridges(
    fill = "gray70",
    scale = 0.6,
    rel_min_height = 0.01,
    color = "black",
    alpha = 0.7,
    show.legend = FALSE
  ) +
  
  
  geom_jitter(
    aes(fill = Sex),
    shape = 21,
    size = 2.5,
    stroke = 0.4,
    color = "black",
    width = 0,
    height = 0.1
  ) +
  
  
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    size = 0.7,
    outlier.shape = NA
  ) +
  
  
  annotate("text", x = -0.6, y = 1, label = "APOE3/3", fontface = "bold", size = 4, hjust = 0) +
  annotate("text", x = -0.6, y = 2, label = "APOE3/4", fontface = "bold", size = 4, hjust = 0) +
  annotate("text", x = -0.6, y = 3, label = "APOE4/4", fontface = "bold", size = 4, hjust = 0) +
  
  
  labs(
    title = bquote("Correlation of APOE4 genotype dosage with myoVC gene activity in APOE-TR mice," ~~ "pcc = 0.25," ~~ italic("P") ~ "=" ~ .(format(pval, digits = 2))),
    x = "Gene activity (GSVA score)",
    y = "Density"
  ) +
  
  scale_fill_manual(
    values = c("F" = "#E64B35", "M" = "#4DBBD5"),
    labels = c("Female", "Male"),
    name = "Sex"
  ) +
  
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )



ggsave("MyoVCmousescore.svg", plot = p, width = 6.5, height = 3, dpi = 300, device = "svg")




