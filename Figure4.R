#Figure 4A, B, Extended Data 5F-I



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
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(readxl)




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

lr_network <- readRDS("lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")





receiver = "Myofibroblast"
expressed_genes_receiver <- get_expressed_genes(receiver, brain, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("Astrocyte", "Endothelial", "Pericyte_1", "Pericyte_2", "SMC_1", "SMC_2", "Myofibroblast")

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, brain, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)


results <- read_excel("myoVCSMC2.mast.xlsx")


DE_table_receiver <- results

DE_table_receiver$X <- results$...1


geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(X)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities


p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity


best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)


vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank())) 



active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()


active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()

order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")



ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both")  

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))








ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank()) +
  scale_fill_gradientn(colors = c("white", "darkorange"),
                       breaks = c(0.02, 0.024, 0.028))  

p_ligand_aupr


ggsave("p_ligand_aupr.svg", plot = p_ligand_aupr, width = 1.5, height = 4, dpi = 300, device = "svg")



active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target


ggsave("p_ligand_target.svg", plot = p_ligand_target, width = 6, height = 3, dpi = 300, device = "svg")



ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

ggsave("p_ligand_receptor.svg", plot = p_ligand_receptor, width = 4, height = 3, dpi = 300, device = "svg")



p_dotplot <- DotPlot(subset(brain, celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

ggsave("p_dotplot.svg", plot = p_dotplot, width = 5, height = 3, dpi = 300, device = "svg")



celltype_order <- levels(Idents(brain)) 



DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = brain,
  condition_colname = "apoe",
  condition_oi = "e4",
  condition_reference = "e3",
  celltype_col = "celltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands,
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ]) 

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc


ggsave("p_lfc.svg", plot = p_lfc, width = 2, height = 4, dpi = 300, device = "svg")



figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

ggsave("combined_plot.pdf", plot = combined_plot, width = 14, height = 10)


brain$celltype <- as.character(brain$celltype)



ligand_type_indication_df <- assign_ligands_to_celltype(brain,
                                                        best_upstream_ligands,
                                                        celltype_col = "celltype") 



active_ligand_target_links_df$target_type <- "Myofibroblast" 


circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df,
                                           cutoff = 0.3) 




ligand_colors <- c("General" = "#377EB8", "Astrocyte" = "#4DAF4A", "Endothelial" = "#984EA3",
                   "Myofibroblast" = "#FF7F00", "Pericyte_1" = "#FFFF33", "Pericyte_2" = "#F781BF",
                   "SMC_1"= "#E41A1C", "SMC_2" = "black") 

target_colors <- c("Myofibroblast" = "#999999") 


vis_circos_obj <- prepare_circos_visualization(circos_links,
                                               ligand_colors = ligand_colors,
                                               target_colors = target_colors,
                                               celltype_order = NULL) 


make_circos_plot(vis_circos_obj, transparency = FALSE,  args.circos.text = list(cex = 0.5)) 


make_circos_plot(vis_circos_obj, transparency = TRUE,  args.circos.text = list(cex = 0.5)) 







par(bg = "transparent")


celltype_order <- unique(circos_links$ligand_type) %>% sort() %>% .[. != "General"] %>% c(., "General")


circos_legend <- ComplexHeatmap::Legend(
  labels = celltype_order,
  background = ligand_colors[celltype_order],
  type = "point",
  grid_height = unit(3, "mm"),
  grid_width = unit(3, "mm"),
  labels_gp = grid::gpar(fontsize = 8)
)

circos_legend_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(circos_legend))

make_circos_plot(vis_circos_obj, transparency = TRUE, args.circos.text = list(cex = 0.5))
p_circos_no_legend <- recordPlot()

cowplot::plot_grid(p_circos_no_legend, circos_legend_grob, rel_widths = c(1, 0.1))


library(cowplot)

library(cowplot)


p_circos_combined <- plot_grid(p_circos_no_legend, circos_legend_grob, rel_widths = c(1, 0.1))


p_with_margins <- ggdraw() +
  draw_plot(p_circos_no_legend)


ggsave("ligand_target_circos.pdf",
       p_with_margins,
       width = 4, height = 3, device = cairo_pdf)


p_legends <- ggdraw() +
  draw_plot(p_circos_combined)


ggsave("circos_legend.pdf",
       p_legends,
       width = 4, height = 3, device = cairo_pdf)



#####################################################################################

#Figure 4C

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


library(SeuratExtend)



brain.da = readRDS("myoVC.mural.clustered.data.rds")
meta = brain.da@meta.data

meta$celltype[meta$celltype == "myoVC"] <- "Myofibroblast"

Idents(brain.da) <- meta$celltype



sce <- as.SingleCellExperiment(brain.da)


sce <- slingshot(sce, clusterLabels = 'celltype', reducedDim = 'UMAP', start.clus = "Pericyte_2")



lin1 <- getLineages(sce, "celltype", start.clus = "Pericyte_2", reducedDim = "UMAP")



dim <- brain.da@reductions$umap@cell.embeddings


crv1 <- getCurves(lin1)



curve <- SlingshotDataSet(crv1)


pseudotime <- slingPseudotime(curve)
cellWeights <- slingCurveWeights(curve)



library(mgcv)
library(slingshot)

curve <- SlingshotDataSet(crv1)




brain <- readRDS("Mural.merged.clustered.data.rds.rds")
meta = brain@meta.data


library(GSEABase)

gene_sets <- getGmt("HALLMARK_TGF_BETA_SIGNALING.v2024.1.Hs.gmt")

gene.ids <- gene_sets[["HALLMARK_TGF_BETA_SIGNALING"]]@geneIds

gene.ids <- list(gene.ids)


brain <- AddModuleScore(brain, 
                        features = gene.ids,
                        name = "tgfb.genes",
                        ctrl = 20)

meta = brain@meta.data




plot_score_gg <- function(score_colname, title, ylab = score_colname) {
  score_df <- brain@meta.data[, score_colname, drop = FALSE]
  score_df$cell <- rownames(score_df)
  
  shared_cells <- intersect(colnames(sce), score_df$cell)
  score_df <- score_df[match(colnames(sce), score_df$cell), ]
  stopifnot(all(score_df$cell == colnames(sce)))
  
  colData(sce)$score <- score_df[[score_colname]]
  
  pseudotime <- slingPseudotime(curve)
  weights <- slingCurveWeights(curve)
  
  df1 <- data.frame(
    pseudotime = pseudotime[, 1],
    score = colData(sce)$score,
    weight = weights[, 1]
  )
  
  df2 <- data.frame(
    pseudotime = pseudotime[, 2],
    score = colData(sce)$score,
    weight = weights[, 2]
  )
  
  df1 <- df1[!is.na(df1$pseudotime), ]
  df2 <- df2[!is.na(df2$pseudotime), ]
  
  gam1 <- gam(score ~ s(pseudotime, k = 7), weights = weight, data = df1)
  gam2 <- gam(score ~ s(pseudotime, k = 7), weights = weight, data = df2)
  
  df1$pred <- predict(gam1, newdata = df1)
  df2$pred <- predict(gam2, newdata = df2)
  
  df1$lineage <- "Pericyte-to-Myofibroblast"
  df2$lineage <- "Pericyte-to-Pericyte"
  
  df_all <- rbind(df1, df2)
  
  p <- ggplot(df_all, aes(x = pseudotime, y = pred, color = lineage)) +
    geom_line(size = 1.5) +
    labs(title = title, x = "Pseudotime", y = ylab, color = "Lineage") +
    scale_color_manual(values = c(
      "Pericyte-to-Myofibroblast" = "#d62728",
      "Pericyte-to-Pericyte" = "#1f77b4"
    )) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black")
    )
  
  return(list(plot = p, df1 = df1, df2 = df2))
}

library(zoo)
library(mgcv)

compare_auc_permutation <- function(df1, df2, label = "", n_perm = 1000, plot_hist = TRUE) {
  
  df1 <- df1[!is.na(df1$pseudotime), ]
  df2 <- df2[!is.na(df2$pseudotime), ]
  
  compute_auc <- function(df) {
    df <- df[order(df$pseudotime), ]
    pt_diff <- diff(df$pseudotime)
    pred_avg <- rollmean(df$pred, 2, align = "left")
    sum(pt_diff * pred_avg)
  }
  
  
  df1$lineage <- "myofibroblast"
  df2$lineage <- "pericyte"
  df_all <- rbind(df1, df2)
  
  
  gam1 <- gam(score ~ s(pseudotime, k = 7), weights = weight, data = df1)
  gam2 <- gam(score ~ s(pseudotime, k = 7), weights = weight, data = df2)
  
  df1$pred <- predict(gam1, newdata = df1)
  df2$pred <- predict(gam2, newdata = df2)
  
  auc1 <- compute_auc(df1)
  auc2 <- compute_auc(df2)
  delta_obs <- auc1 - auc2
  
  
  set.seed(42)
  delta_null <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    perm_labels <- sample(df_all$lineage)
    df_all$perm_lineage <- perm_labels
    
    b1 <- df_all[df_all$perm_lineage == "myofibroblast", ]
    b2 <- df_all[df_all$perm_lineage == "pericyte", ]
    
    gam_b1 <- gam(score ~ s(pseudotime, k = 7), weights = weight, data = b1)
    gam_b2 <- gam(score ~ s(pseudotime, k = 7), weights = weight, data = b2)
    
    b1$pred <- predict(gam_b1, newdata = b1)
    b2$pred <- predict(gam_b2, newdata = b2)
    
    delta_null[i] <- compute_auc(b1) - compute_auc(b2)
  }
  
  
  p_val <- mean(abs(delta_null) >= abs(delta_obs))
  ci <- quantile(delta_null, c(0.025, 0.975))
  
  cat("=== Permutation AUC Comparison for:", label, "===\n")
  cat("Observed Δ-AUC:", round(delta_obs, 3), "\n")
  cat("95% null CI:", round(ci[1], 3), "to", round(ci[2], 3), "\n")
  cat("Permutation-based p-value:", round(p_val, 4), "\n\n")
  
  if (plot_hist) {
    hist(delta_null, breaks = 50, main = paste("Permutation Null Δ-AUC:", label),
         xlab = "Δ-AUC (shuffled lineages)", col = "gray", border = "white")
    abline(v = delta_obs, col = "red", lwd = 2)
  }
  
  return(list(
    delta_obs = delta_obs,
    auc_myoVC = auc1,
    auc_pericyte = auc2,
    delta_null = delta_null,
    null_ci = ci,
    p_value = p_val
  ))
}


out_tgfb <- plot_score_gg("tgfb.genes1", "TGF-β Signature Across Pseudotime")
print(out_tgfb$plot)

perm_tgfb <- compare_auc_permutation(out_tgfb$df1, out_tgfb$df2, label = "TGF-β Signature")





########################################################################################


#Figure 4D




brain <- readRDS("myoVC.mural.clustered.data.rds")

brain$celltype[brain$celltype == "myoVC"] <- "Myofibroblast"


Idents(brain) <- brain$celltype

brain <- subset(brain, subset = celltype == "Pericyte_1" | celltype == "Pericyte_2")



expr <- brain[["RNA"]]$data



top10_tgfb_genes <- c(
  "TGFBR1",
  "TGFBR2",
  "TGFB1",
  "TGFB2")


brain <- AddModuleScore(brain,
                        features = list(top10_tgfb_genes),
                        name = "tgfb.genes",
                        ctrl = 20)  


meta = brain@meta.data

library(dplyr)
library(ggplot2)


pseudobulk_df <- meta %>%
  group_by(subject, apoe) %>%
  summarise(mean_tgfb_score = mean(tgfb.genes1, na.rm = TRUE), .groups = "drop")


min_score <- min(pseudobulk_df$mean_tgfb_score, na.rm = TRUE)


pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_tgfb_score_shifted = mean_tgfb_score - min_score)


median_e3 <- pseudobulk_df %>%
  filter(apoe == "e3") %>%
  summarise(med = median(mean_tgfb_score_shifted, na.rm = TRUE)) %>%
  pull(med)

pseudobulk_df <- pseudobulk_df %>%
  mutate(mean_tgfb_score_normalized = mean_tgfb_score_shifted / median_e3)



median_e4 <- pseudobulk_df %>%
  filter(apoe == "e4") %>%
  summarise(med = median(mean_tgfb_score_normalized, na.rm = TRUE)) %>%
  pull(med)

print(median_e4)

write.csv(pseudobulk_df, file = "pseudobulk.tgfbscore.csv")


##########################################################################

#Figure 4F


library(limma)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)




fpkm_data <-read.csv("FPKM.csv")


cleaned_data <- fpkm_data[!grepl("^(\\d{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec))$", fpkm_data$Genes), ]



fpkm_data <- cleaned_data

write.table(
  fpkm_data,
  "bulk.cibersort.txt",   
  sep = "\t",           
  row.names = FALSE,    
  col.names = FALSE,     
  quote = FALSE        
)


fpkm_data <- as.data.frame(fpkm_data)


rownames(fpkm_data) <- fpkm_data$Genes


fpkm_data <- fpkm_data[, -which(colnames(fpkm_data) == "Genes")]



metadata <- data.frame(
  sample = colnames(fpkm_data),
  condition = c("e3", "e3", "e3", "e4", "e4", "e4")
)



log_fpkm <- log2(fpkm_data + 0.1)



design <- model.matrix(~ condition, data = metadata)


rownames(design) <- colnames(log_fpkm)


weights <- arrayWeights(log_fpkm, design)


fit <- lmFit(log_fpkm, design, weights = weights)


fit <- eBayes(fit, trend = TRUE, robust = TRUE)


results <- topTable(fit, coef = "conditione4", number = Inf)  

results$X <- rownames(results)





library(clusterProfiler)
library(ggplot2)
library(enrichplot)


filtered_DEGs_UP <- subset(results, subset = adj.P.Val < 0.05 & logFC > 0.25)
gene_list_up <- filtered_DEGs_UP$X

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


filtered_DEGs_DOWN <- subset(results, subset = adj.P.Val < 0.05 & logFC < -0.25)
gene_list_down <- filtered_DEGs_DOWN$X

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



library(ggplot2)
library(dplyr)
library(grid)


tgf_up_df <- as.data.frame(cP_GO_UP) %>%
  filter(grepl("transforming", Description, ignore.case = TRUE)) %>%
  arrange(p.adjust)


tgf_up_df$Direction <- "Upregulated"
tgf_up_df$signed_log10_p <- -log10(tgf_up_df$p.adjust)


tgf_up_df$Description <- factor(tgf_up_df$Description, levels = rev(tgf_up_df$Description))

p_tgf <- ggplot(tgf_up_df, aes(x = signed_log10_p, y = Description, fill = Direction == "Upregulated")) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#D73027")) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  labs(
    title = "Top Upregulated GO Terms Involving TGF-β",
    x = expression(Signed~ -log[10]~(Adjusted~italic(p))),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 9, face = "bold", color = "black", lineheight = 1.1),
    axis.text.x = element_text(size = 9, face = "bold", color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks.x = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    plot.title = element_text(size = 11, face = "bold", color = "black", hjust = 0.5),
    legend.position = "none"
  )


ggsave("TGFb_GO_Upregulated_Clean.svg", plot = p_tgf, width = 6, height = 3.5, dpi = 600, device = "svg")
