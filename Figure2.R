#Figure 2C-E

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

colors <- c("#006884", "#058799", "#EE9E80", "#E74C2E", "red")


cell_order <- c(
  WhichCells(brain.da, idents = "Pericyte_2"),
  WhichCells(brain.da, idents = "Pericyte_1"),
  WhichCells(brain.da, idents = "SMC_2"),
  WhichCells(brain.da, idents = "SMC_1"),
  WhichCells(brain.da, idents = "Myofibroblast")  
)


p <- DimPlot(
  brain.da,
  cells = cell_order,
  cols = colors,
  shuffle = FALSE  
)




curve1 <- SlingshotDataSet(crv1)@curves$Lineage1$s
curve2 <- SlingshotDataSet(crv1)@curves$Lineage2$s


curve1_df <- data.frame(x = curve1[,1], y = curve1[,2])
curve2_df <- data.frame(x = curve2[,1], y = curve2[,2])




p <- p + 
  geom_path(data = curve1_df, aes(x = x, y = y), color = "black", size = 2) +
  geom_path(data = curve2_df, aes(x = x, y = y), color = "black", size = 2)


ggsave("LineageDimPlot.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")



sling <- slingPseudotime(sce)
sling <- as.data.frame(sling)
brain.da$Lineage1 <- sling$Lineage1
brain.da$Lineage2 <- sling$Lineage2


p <- FeaturePlot(brain.da, features = "Lineage1")






p <- p + 
  geom_path(data = curve1_df, aes(x = x, y = y), color = "black", size = 2)


ggsave("Lineage1.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")




p <- FeaturePlot(brain.da, features = "Lineage2")




p <- p + 
  geom_path(data = curve2_df, aes(x = x, y = y), color = "black", size = 2)

ggsave("Lineage2.svg", plot = p, width = 4, height = 3, dpi = 300, device = "svg")








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


brain$myofibroblast.genes1 <- pmin(brain$contraction.genes1, brain$ecm.genes1)



brain <- AddModuleScore(brain, 
                        features = list(c("PDGFRB", "ANPEP", "CSPG4", "KCNJ8", "SLC20A2", "SLC6A1")),
                        name = "pericyte.genes",
                        ctrl = 20)

meta = brain@meta.data







curve <- SlingshotDataSet(crv1)


pseudotime <- slingPseudotime(curve)
cellWeights <- slingCurveWeights(curve)



library(mgcv)
library(slingshot)

curve <- SlingshotDataSet(crv1)





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





out_joint <- plot_score_gg("myofibroblast.genes1", "Myofibroblast Signature Across Pseudotime")
print(out_joint$plot)

perm_joint <- compare_auc_permutation(out_joint$df1, out_joint$df2, label = "Myofibroblast Signature")



out_peri <- plot_score_gg("pericyte.genes1", "Pericyte Signature Across Pseudotime")
print(out_peri$plot)

perm_peri <- compare_auc_permutation(out_peri$df1, out_peri$df2, label = "Pericyte Signature")




library(condiments)

df <- NULL
df <- as.data.frame(df)

library(ggplot2)
library(cowplot)  


df <- slingCurveWeights(sce, as.probs = TRUE)[, 1]
df <- as.data.frame(df)
df$weight_1 <- df$df
df$conditions <- meta$apoe


custom_colors <- c("e3" = "blue", "e4" = "red")


p1 <- ggplot(df, aes(x = weight_1, fill = conditions)) +
  geom_density(alpha = 0.5, color = "black") +  
  scale_fill_manual(
    values = custom_colors, 
    aesthetics = "fill",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(x = "Curve weight for the first lineage", y = "Density") +
  theme_classic() +  
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.text = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 14, face = "bold"),  
    legend.title = element_blank(),  
    legend.position = "top"  
  )


df$weight_2 <- slingCurveWeights(sce, as.probs = TRUE)[, 2]

p2 <- ggplot(df, aes(x = weight_2, fill = conditions)) +
  geom_density(alpha = 0.5, color = "black") +  
  scale_fill_manual(
    values = custom_colors, 
    aesthetics = "fill",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(x = "Curve weight for the second lineage", y = "Density") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  )


print(p1)
print(p2)

ggsave("Lineage1.weights.svg", plot = p1, width = 4, height = 3, dpi = 300, device = "svg")
ggsave("Lineage2.weights.svg", plot = p2, width = 4, height = 3, dpi = 300, device = "svg")


set.seed(12)
dif_res <- differentiationTest(sce, conditions = df$condition, global = FALSE, pairwise = TRUE)

knitr::kable(dif_res)

##########################################################################################

#Figure 2G, Extended Data 4A-C





library(Seurat)
library(ggplot2)
library(dplyr)


brain <- readRDS("Mural.merged.clustered.data.rds.rds")
brain <- subset(brain, subset = AD.status == 'nonAD')
meta = brain@meta.data



expr <- FetchData(brain, vars = c("CSPG4", "ACTA2"))
meta <- brain@meta.data
meta$high_both <- expr$CSPG4 > 1.5 & expr$ACTA2 > 1.5


subject_level <- meta %>%
  group_by(subject) %>%
  summarise(
    n_transitioning = sum(high_both),
    total_cells = n(),
    prop_transitioning = n_transitioning / total_cells,
    apoe = dplyr::first(apoe),
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
  cbind(n_transitioning, total_cells - n_transitioning) ~ apoe + ns(age, df = 3) + sex + dataset,
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
    title = "Predictors of transitioning mural cells (CSPG4+/ACTA2+) (Random Effects Model)",
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


ggsave("transitioning_forestplot_with_interaction.svg", plot = p_random, width = 6, height = 4.5, dpi = 300, device = "svg")




library(dplyr)
library(ggplot2)




subject_level$pred <- predict(model_covariate, type = "response")


write.csv(subject_level, file = "transitioningpred.csv")






library(Seurat)
library(ggplot2)
library(dplyr)


meta <- brain@meta.data
meta$CSPG4 <- FetchData(brain, vars = "CSPG4")[,1]
meta$ACTA2 <- FetchData(brain, vars = "ACTA2")[,1]


meta_nonAD <- meta %>% filter(AD.status == "nonAD")


meta_nonAD$apoe <- factor(meta_nonAD$apoe, levels = c("e3", "e4"), labels = c("APOE3", "APOE4"))



library(ggplot2)
library(dplyr)


meta_plot <- meta_nonAD %>%
  filter((CSPG4 > 0 | ACTA2 > 0) & is.finite(CSPG4) & is.finite(ACTA2))

library(ggplot2)

p <- ggplot(meta_plot, aes(x = CSPG4, y = ACTA2, color = apoe)) +
  stat_density_2d(aes(fill = apoe), geom = "polygon", alpha = 0.2, contour = TRUE, bins = 20) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  scale_color_manual(
    values = c("APOE3" = "#46B5A9", "APOE4" = "#D45176")
  ) +
  scale_fill_manual(
    values = c("APOE3" = "#46B5A9", "APOE4" = "#D45176")
  ) +
  labs(
    title = "CSPG4 (NG2) vs ACTA2 (SMA) in Non-AD Individuals",
    x = "CSPG4 expression (log-normalized)",
    y = "ACTA2 expression (log-normalized)",
    color = "APOE Genotype",
    fill = "APOE Genotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),          
    axis.ticks = element_line(),           
    axis.line = element_line(),            
    text = element_text(color = "black"),  
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "plain"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )


ggsave("transitioning.scatterplot.densities.svg", plot = p, width = 6, height = 4, dpi = 300, device = "svg")


p <- ggplot(meta_plot, aes(x = CSPG4, y = ACTA2, color = apoe)) +
  geom_point(size = 2.5, alpha = 0.3) +  
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.5, linetype = "dashed") +
  scale_color_manual(
    values = c("APOE3" = "#46B5A9", "APOE4" = "#D45176")
  ) +
  labs(
    title = "CSPG4 (NG2) vs ACTA2 (SMA) in Non-AD Individuals",
    x = "CSPG4 expression (log-normalized)",
    y = "ACTA2 expression (log-normalized)",
    color = "APOE Genotype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    axis.line = element_line(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "plain"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )

ggsave("transitioning.scatterplot.points.svg", plot = p, width = 6, height = 4, dpi = 300, device = "svg")

#############################################################################################

#Extended Data 4E-I



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

write.csv(results, file = 'results.csv')



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


cP_GO_UP <- simplify(cP_GO_UP)


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


cP_GO_DOWN <- simplify(cP_GO_DOWN)


up_df <- as.data.frame(cP_GO_UP) %>%
  arrange(p.adjust) %>%
  dplyr::slice_head(n = 5)
up_df$Direction <- "Upregulated"


down_df <- as.data.frame(cP_GO_DOWN) %>%
  arrange(p.adjust) %>%
  dplyr:: slice_head(n = 5)
down_df$Direction <- "Downregulated"



up_df$signed_log10_p <- -log10(up_df$p.adjust)  
down_df$signed_log10_p <- log10(down_df$p.adjust)  


combined_df <- rbind(up_df, down_df)


combined_df <- combined_df[order(combined_df$signed_log10_p, decreasing = TRUE), ]
combined_df$Description <- factor(combined_df$Description, levels = combined_df$Description)



combined_df$Description <- factor(
  combined_df$Description,
  levels = rev(combined_df$Description)
)

library(ggplot2)
library(grid)


p <- ggplot(combined_df, aes(x = signed_log10_p, y = Description, fill = Direction == "Upregulated")) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("FALSE" = "#4682B4", "TRUE" = "#B22222")) +
  scale_x_continuous(limits = c(-14, 35), expand = c(0, 0)) +
  geom_vline(xintercept = 0, color = "black", size = 1, linetype = "solid") +
  geom_hline(yintercept = 5.5, linetype = "dashed", color = "black", size = 0.4) +  
  labs(
    x = "Signed -Log10 Adjusted P-value",
    y = NULL,
    title = "Top 10 Upregulated and Downregulated Pathways"
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
    legend.position = "none"
  )


ggsave("Pathways.svg", plot = p, width = 6, height = 4.5, dpi = 300, device = "svg")





df <- results


library(ggplot2)
library(ggrepel)


df$significance <- ifelse(abs(df$logFC) > 0.25 & df$adj.P.Val < 0.05,
                          ifelse(df$logFC > 0, "Up", "Down"),
                          "Not Sig")


genes_of_interest <- c("FN1", "COL1A1", "COL4A1", "COL1A2", "COL3A1", "COL4A2", "ACTA2", "CDH11", "S100A6", "ITGA1", "ITGA11")  


df$label <- ifelse(df$X %in% genes_of_interest, df$X, NA)


volcano_plot <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  
  geom_point(alpha = 0.8, size = 1.5) +
  
  geom_point(data = df[!is.na(df$label), ], aes(x = logFC, y = -log10(adj.P.Val)), 
             color = "black", size = 6, shape = 21, fill = "yellow") +
  
  geom_label_repel(aes(label = label),
                   size = 5, 
                   fontface = "bold",
                   color = "black", 
                   fill = "white", 
                   box.padding = 0.6, 
                   point.padding = 0.3, 
                   segment.color = "black", 
                   segment.size = 1.2, 
                   max.overlaps = Inf, 
                   na.rm = TRUE) +
  
  scale_color_manual(values = c("Up" = "#B22222", "Down" = "#4682B4", "Not Sig" = "gray80")) +
  
  labs(
    x = "logFC",
    y = "-log10(Adjusted P-Value)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),  
    axis.ticks = element_line(color = "black", linewidth = 0.4), 
    axis.text = element_text(color = "black", size = 12),  
    axis.title = element_text(size = 14)  
  ) 


print(volcano_plot)


ggsave(filename = "final_volcano_plot_thinner_lines.svg", 
       plot = volcano_plot, 
       device = "svg", 
       width = 6, 
       height = 4.5)








