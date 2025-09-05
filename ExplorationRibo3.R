# Cargar librerías necesarias
library(tidyverse)
library(ggplot2)
library(readr)
library(ggrepel)

# Seleccionar archivos CSV de forma interactiva
cat("Selecciona el archivo de ribo-seq\n")
ribo_file <- file.choose()

cat("Selecciona el archivo de RNA-seq\n")
rna_file <- file.choose()

cat("Selecciona el archivo de anotación genética\n")
meta_file <- file.choose()

# Cargar los datasets
ribo <- read_csv("riboseq.csv")
rna <- read_csv("RNAseq_NpvsSomata.csv")
meta <- read_csv("SFARI.csv")

# Asegurar que los nombres de los genes estén en minúsculas
ribo <- ribo %>% mutate(gene = tolower(gene))
rna  <- rna  %>% mutate(gene = tolower(gene))
meta <- meta %>%
  rename(gene = `gene-symbol`) %>%
  mutate(gene = tolower(gene))
meta<- meta %>%
  filter(!is.na(gene) & !is.na(`gene-score`))

# Filtrar por significancia estadística (padj < 0.05)
ribo_sig <- ribo %>% filter(ribo_padj < 0.05)
rna_sig  <- rna  %>% filter(rna_padj  < 0.05)

# Unir con datos de anotación
ribo_sig <- ribo_sig %>% 
  inner_join(meta %>% select(gene, gene_score = `gene-score`, genetic_category = `genetic-category`, syndromic), 
             by = "gene")

rna_sig <- rna_sig %>% 
  inner_join(meta %>% select(gene, gene_score = `gene-score`, genetic_category = `genetic-category`, syndromic), 
             by = "gene")

# Identificar los 5 genes extremos por gene_score (ribo)
ribo_labels <- ribo_sig %>%
  group_by(gene_score) %>%
  slice_max(ribo_log2FC, n = 5, with_ties = FALSE) %>%
  bind_rows(
    ribo_sig %>%
      group_by(gene_score) %>%
      slice_min(ribo_log2FC, n = 5, with_ties = FALSE)
  )

# Gráfico ribo-seq con etiquetas
ggplot(ribo_sig, aes(x = gene_score, y = ribo_log2FC)) +
  geom_point(alpha = 0.7, aes(color=syndromic)) +
  geom_text_repel(data = ribo_labels, aes(label = gene), size = 3, max.overlaps = 100) +
  facet_wrap(~ localization) +
  labs(title = "Ribo-seq: log2FC vs gene-score",
       x = "Gene Score",
       y = "Ribo Log2FC") +
  theme_minimal()

# Gráfico RNA-seq: log2FC vs gene_score, facetado por localization
ggplot(rna_sig, aes(x = gene_score, y = rna_log2FC)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  labs(title = "RNA-seq: log2FC vs gene-score",
       x = "Gene Score",
       y = "RNA Log2FC") +
  theme_minimal()

# Genes en común significativos en ambos datasets
common_genes <- inner_join(ribo_sig, rna_sig, by = "gene", suffix = c("_ribo", "_rna"))


# Crear columna binaria de syndromic
common_genes <- common_genes %>%
  mutate(syndromic_binary = if_else(tolower(syndromic_ribo) == "yes", 1, 0))


# Etiquetas para ribo: top 5 max y min ribo_log2FC por syndromic_ribo
labels_ribo <- common_genes %>%
  group_by(syndromic_ribo) %>%
  slice_max(ribo_log2FC, n = 5, with_ties = FALSE) %>%
  bind_rows(
    common_genes %>%
      group_by(syndromic_ribo) %>%
      slice_min(ribo_log2FC, n = 5, with_ties = FALSE)
  ) %>%
  mutate(type = "ribo")

# Etiquetas para rna: top 5 max y min rna_log2FC por syndromic_ribo
labels_rna <- common_genes %>%
  group_by(syndromic_ribo) %>%
  slice_max(rna_log2FC, n = 5, with_ties = FALSE) %>%
  bind_rows(
    common_genes %>%
      group_by(syndromic_ribo) %>%
      slice_min(rna_log2FC, n = 5, with_ties = FALSE)
  ) %>%
  mutate(type = "rna")

# Unir etiquetas
labels_combined <- bind_rows(labels_ribo, labels_rna)

ggplot(common_genes) +
  geom_point(aes(x = gene_score_ribo, y = ribo_log2FC), alpha = 0.7, size = 2, shape = 16, color = "firebrick") +
  geom_point(aes(x = gene_score_rna, y = rna_log2FC), alpha = 0.7, size = 2, shape = 17, color = "steelblue") +
  geom_text_repel(data = labels_combined %>% filter(type == "ribo"),
                  aes(x = gene_score_ribo, y = ribo_log2FC, label = gene),
                  color = "firebrick", size = 3, max.overlaps = 100) +
  facet_wrap(~ syndromic_ribo) +
  labs(title = "Comparación RNA-seq y Ribo-seq\nFacetas por clasificación 'syndromic'",
       x = "Gene Score",
       y = "Log2FC",
       caption = "Símbolos: ● Ribo-seq (rojo) | ▲ RNA-seq (azul)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))


# Crear columnas que indican si un gen es top 5 max o min para ribo y rna, por syndromic
common_genes_summary <- common_genes %>%
  group_by(syndromic_ribo) %>%
  mutate(
    top5_max_ribo = ribo_log2FC >= quantile(ribo_log2FC, 0.95),
    top5_min_ribo = ribo_log2FC <= quantile(ribo_log2FC, 0.05),
    top5_max_rna  = rna_log2FC  >= quantile(rna_log2FC, 0.95),
    top5_min_rna  = rna_log2FC  <= quantile(rna_log2FC, 0.05)
  ) %>%
  ungroup() %>%
  filter(top5_max_ribo | top5_min_ribo | top5_max_rna | top5_min_rna) %>%
  select(gene, syndromic = syndromic_ribo, gene_score_ribo, ribo_log2FC, rna_log2FC,
         top5_max_ribo, top5_min_ribo, top5_max_rna, top5_min_rna) %>%
  arrange(syndromic, gene_score_ribo, desc(ribo_log2FC))

# Mostrar tabla resumen
print(common_genes_summary)


# Unir con meta para obtener gene_score y syndromic para ribo (si no lo tienes ya)
ribo_all <- ribo %>%
  inner_join(meta %>% select(gene, gene_score = `gene-score`), by = "gene")

# Identificar top 5 genes con mayor gene_score por localization
top_genes <- ribo_all %>%
  group_by(localization) %>%
  slice_max(gene_score, n = 5, with_ties = FALSE)

# Seleccionar top 5 genes por cada combinación localization x gene_score, basados en ribo_log2FC
top_genes <- ribo_all %>%
  group_by(localization, gene_score) %>%
  slice_max(ribo_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

# Crear variable de color según significancia (ribo_padj <= 0.05 es significativo)
ribo_all <- ribo_all %>%
  mutate(significance = if_else(ribo_padj <= 0.05, "Significant", "Not Significant"))

# Plot with jitter to separate points, color by significance, and gene name labels for top genes
ggplot(ribo_all, aes(x = factor(gene_score), y = ribo_log2FC, color = significance)) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.3) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 2, color = "black", max.overlaps = 70) +
  facet_wrap(~ localization) +
  scale_color_manual(values = c("Significant" = "darkolivegreen4", "Not Significant" = "red")) +
  labs(title = "Ribo-seq: log2FC vs gene_score (unfiltered p-value)",
       subtitle = "Red points: adjusted p-value > 0.05; Labels: top 5 ribo_log2FC per localization and gene_score",
       x = "Gene Score",
       y = "Ribo Log2 Fold Change",
       color = "Significance") +
  theme_minimal()

#### ESTE SI QUE SI #####

plots <- list()

for (loc in unique(ribo_all$localization)) {
  ribo_subset <- ribo_all %>% filter(localization == loc)
  top_subset <- top_genes %>% filter(localization == loc)
  
  n <- nrow(top_subset)
  if (n > 0) {
    if (loc == "SMT-translation-up") {
      # Nudges generales
      top_subset <- top_subset %>%
        mutate(
          nudge_x = 0.15,
          nudge_y = 0.4,
          segment.angle = 45,
          segment.curvature = 0.2
        )
      
      # Ajuste específico para hnrnph2 y ank3 para separarlos más
      for (gene_name in c("hnrnph2", "ank3")) {
        idx <- which(tolower(top_subset$gene) == tolower(gene_name))
        if (length(idx) == 1) {
          if (gene_name == "hnrnph2") {
            top_subset$nudge_x[idx] <- 0.05
            top_subset$nudge_y[idx] <- 0.3
            top_subset$segment.angle[idx] <- 60
            top_subset$segment.curvature[idx] <- 0.3
          } else if (gene_name == "ank3") {
            top_subset$nudge_x[idx] <- 0.02   # Muy cerca horizontalmente
            top_subset$nudge_y[idx] <- 0.1    # Muy cerca verticalmente
            top_subset$segment.angle[idx] <- 75
            top_subset$segment.curvature[idx] <- 0.1
          }
        }
      }
    } else {
      # Nudges estándar para otras localizations
      top_subset <- top_subset %>%
        mutate(
          nudge_y = rep(c(0.5, -0.5, 1, -1, 0.7, -0.7), length.out = n),
          nudge_x = rep(c(0.3, -0.3, 0.5, -0.5, 0.2, -0.2), length.out = n),
          segment.angle = 90,
          segment.curvature = 0
        )
    }
  }
  
  if (loc == "other") {
    p <- ggplot(ribo_subset, aes(x = factor(gene_score), y = ribo_log2FC, color = ribo_padj)) +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_text_repel(
        data = top_subset,
        aes(label = gene),
        size = 3,
        color = "black",
        box.padding = 0.8,
        point.padding = 0.5,
        nudge_x = top_subset$nudge_x,
        nudge_y = top_subset$nudge_y,
        segment.size = 0.4,
        segment.color = "gray40",
        segment.angle = top_subset$segment.angle,
        segment.curvature = top_subset$segment.curvature,
        min.segment.length = 0.05,
        max.overlaps = Inf,
        direction = "both"
      ) +
      scale_color_viridis_c(option = "plasma", direction = -1, name = "Adjusted p-value") +
      labs(
        title = paste("Ribo-seq in", loc),
        subtitle = "Top 5 ribo_log2FC per gene score; color = adjusted p-value",
        x = "Gene Score",
        y = "Ribo Log2 Fold Change"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
      )
  } else {
    p <- ggplot(ribo_subset, aes(x = factor(gene_score), y = ribo_log2FC, color = significance)) +
      geom_point(alpha = 0.6, size = 1.5) +
      geom_text_repel(
        data = top_subset,
        aes(label = gene),
        size = 3,
        color = "black",
        box.padding = 0.8,
        point.padding = 0.5,
        nudge_x = top_subset$nudge_x,
        nudge_y = top_subset$nudge_y,
        segment.size = 0.4,
        segment.color = "gray40",
        segment.angle = top_subset$segment.angle,
        segment.curvature = top_subset$segment.curvature,
        min.segment.length = 0.05,
        max.overlaps = Inf,
        direction = "both"
      ) +
      scale_color_manual(values = c("Significant" = "black", "Not Significant" = "red")) +
      labs(
        title = paste("Ribo-seq in", loc),
        subtitle = "Top 5 ribo_log2FC per gene score (red = adj. p-value > 0.05)",
        x = "Gene Score",
        y = "Ribo Log2 Fold Change",
        color = "Significance"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
      )
  }
  
  plots[[loc]] <- p
}

# Mostrar gráficos
for (plot in plots) {
  print(plot)
}


#### otras visuals ####


# Suponiendo que ya tienes:
# ribo_all (riboseq dataframe con columna 'gene')
# rna (rnaseq dataframe con columna 'gene')

# Asegúrate que 'gene' esté en minúscula en ambos
ribo_all$gene <- tolower(ribo_all$gene)
rna$gene <- tolower(rna$gene)

# Obtener top 10 genes por |ribo_log2FC| en cada gene_score
top10_by_gene_score <- ribo_all %>%
  group_by(gene_score) %>%
  slice_max(order_by = abs(ribo_log2FC), n = 10, with_ties = FALSE) %>%
  ungroup()

# Combinar con RNA-seq (sin gene_info)
top10_combined <- top10_by_gene_score %>%
  left_join(rna, by = "gene", suffix = c("_ribo", "_rna")) %>%
  mutate(
    ribo_signif = ifelse(ribo_padj < 0.05, "Significant", "Not Significant"),
    rna_signif = ifelse(rna_padj < 0.05, "Significant", "Not Significant"),
    gene_score = factor(gene_score)
  ) %>%
  select(gene, gene_score, localization, ribo_log2FC, rna_log2FC, ribo_signif, rna_signif)

# Pivotar para gráfico combinado
plot_data <- top10_combined %>%
  pivot_longer(
    cols = c(ribo_log2FC, rna_log2FC),
    names_to = "assay",
    values_to = "log2FC"
  ) %>%
  mutate(
    signif = ifelse(assay == "ribo_log2FC", ribo_signif, rna_signif),
    assay = recode(assay, ribo_log2FC = "Ribo-seq", rna_log2FC = "RNA-seq")
  )

# Gráfico propuesto
ggplot(plot_data, aes(x = gene_score, y = log2FC, color = signif)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.3) +
  facet_grid(assay ~ localization) +
  scale_color_manual(values = c("Significant" = "black", "Not Significant" = "red")) +
  labs(
    title = "Top 10 Ribo-seq Genes by Gene Score Compared with RNA-seq",
    x = "Gene Score",
    y = "Log2 Fold Change",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



### volcano ####

ribo_all <- ribo_all %>%
  left_join(meta %>% select(gene, syndromic), by = "gene")

# Etiquetamos genes significativos y si son sindrómicos
ribo_volcano <- ribo_all %>%
  mutate(
    syndromic = factor(syndromic, labels = c("Non-Syndromic", "Syndromic")),
    signif = case_when(
      ribo_padj < 0.05 & abs(ribo_log2FC) > 1 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Selección opcional de genes para etiquetar
genes_to_label <- ribo_volcano %>%
  filter(signif == "Significant", syndromic == "Syndromic") %>%
  arrange(ribo_padj) %>%
  slice_head(n = 15)  # top 15 síndromicos más significativos

# Volcano plot
ggplot(ribo_volcano, aes(x = ribo_log2FC, y = -log10(ribo_padj))) +
  geom_point(aes(color = syndromic, shape = signif), alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_text_repel(data = genes_to_label, aes(label = gene), size = 3, box.padding = 0.3) +
  scale_color_manual(values = c("Non-Syndromic" = "steelblue", "Syndromic" = "firebrick")) +
  labs(
    title = "Volcano Plot of Ribo-seq Genes",
    x = "Log2 Fold Change (Ribo-seq)",
    y = "-log10(adj p-value)",
    color = "Gene Category",
    shape = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


library(ggplot2)
library(dplyr)
library(ggrepel)

# Etiquetamos genes significativos y si son sindrómicos
ribo_volcano <- ribo_all %>%
  mutate(
    syndromic = factor(syndromic, labels = c("Non-Syndromic", "Syndromic")),
    signif = case_when(
      ribo_padj < 0.05 & abs(ribo_log2FC) > 1 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Selección de genes síndromicos significativos para etiquetar (opcional)
genes_to_label <- ribo_volcano %>%
  filter(signif == "Significant", syndromic == "Syndromic") %>%
  arrange(ribo_padj) %>%
  slice_head(n = 15)

# Volcano plot con facet por syndromic y localization
ggplot(ribo_volcano, aes(x = ribo_log2FC, y = -log10(ribo_padj))) +
  geom_point(aes(color = syndromic, shape = signif), alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c("Non-Syndromic" = "steelblue", "Syndromic" = "firebrick")) +
  labs(
    title = "Volcano Plot of Ribo-seq Genes",
    x = "Log2 Fold Change (Ribo-seq)",
    y = "-log10(adj p-value)",
    color = "Gene Category",
    shape = "Significance"
  ) +
  facet_grid(syndromic ~ localization) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

library(ggplot2)
library(dplyr)
library(ggrepel)

ribo_plot <- ribo_all %>%
  mutate(
    syndromic = factor(syndromic, labels = c("Non-Syndromic", "Syndromic")),
    signif = case_when(
      ribo_padj < 0.05 & abs(ribo_log2FC) > 1 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    signif = factor(signif, levels = c("Significant", "Not Significant"))
  )

# Etiquetamos los 20 genes más relevantes entre los significativos
genes_to_label <- ribo_plot %>%
  filter(signif == "Significant") %>%
  arrange(desc(gene_score)) %>%
  slice_head(n = 20)

# Gráfico sin tamaño variable
ggplot(ribo_plot, aes(x = gene_score, y = ribo_log2FC)) +
  geom_point(aes(color = syndromic), alpha = 0.75, size = 2.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 30
  ) +
  scale_color_manual(values = c("Non-Syndromic" = "steelblue", "Syndromic" = "firebrick")) +
  facet_grid(rows = vars(signif), cols = vars(localization), scales = "free") +
  labs(
    title = "Gene Score vs Ribo-seq Fold Change by Significance and Localization",
    subtitle = "Color = Syndromic",
    x = "Gene Score",
    y = "Log2 Fold Change (Ribo-seq)",
    color = "Syndromic Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 11)
  )
#### highlighting por syndromic o no ####

# Filtrar genes syndromic más destacados
genes_to_label_syndromic <- ribo_plot %>%
  filter(signif == "Significant", syndromic == "Syndromic") %>%
  arrange(desc(gene_score)) %>%
  slice_head(n = 20)

ggplot(ribo_plot, aes(x = gene_score, y = ribo_log2FC)) +
  geom_point(aes(color = syndromic), alpha = 0.75, size = 2.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_text_repel(
    data = genes_to_label_syndromic,
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 30
  ) +
  scale_color_manual(values = c("Non-Syndromic" = "steelblue", "Syndromic" = "firebrick")) +
  facet_grid(rows = vars(signif), cols = vars(localization), scales = "free") +
  labs(
    title = "Gene Score vs Ribo-seq FC (Syndromic genes highlighted)",
    x = "Gene Score",
    y = "Log2 Fold Change (Ribo-seq)",
    color = "Syndromic Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")


# Filtrar genes non-syndromic más destacados
genes_to_label_non_syndromic <- ribo_plot %>%
  filter(signif == "Significant", syndromic == "Non-Syndromic") %>%
  arrange(desc(gene_score)) %>%
  slice_head(n = 20)

ggplot(ribo_plot, aes(x = gene_score, y = ribo_log2FC)) +
  geom_point(aes(color = syndromic), alpha = 0.75, size = 2.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_text_repel(
    data = genes_to_label_non_syndromic,
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 30
  ) +
  scale_color_manual(values = c("Non-Syndromic" = "steelblue", "Syndromic" = "firebrick")) +
  facet_grid(rows = vars(signif), cols = vars(localization), scales = "free") +
  labs(
    title = "Gene Score vs Ribo-seq FC (Non-Syndromic genes highlighted)",
    x = "Gene Score",
    y = "Log2 Fold Change (Ribo-seq)",
    color = "Syndromic Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")


###correlation####
# Combinar todos los genes ribo con RNA-seq (sin gene_info)
combined_all <- ribo_all %>%
  left_join(rna, by = "gene", suffix = c("_ribo", "_rna")) %>%
  mutate(
    ribo_signif = ifelse(ribo_padj < 0.05, "Significant", "Not Significant"),
    rna_signif = ifelse(rna_padj < 0.05, "Significant", "Not Significant"),
    gene_score = factor(gene_score)
  )

# Correlation plot con todos los genes
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)  # para stat_cor()

p_by_score_all <- ggplot(combined_all, aes(x = ribo_log2FC, y = rna_log2FC)) +
  geom_point(aes(color = gene_score), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  facet_wrap(~ gene_score) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Correlation between Ribo-seq and RNA-seq by Gene Score (All Genes)",
    x = "Ribo-seq Log2 Fold Change",
    y = "RNA-seq Log2 Fold Change"
  ) +
  theme_minimal(base_size = 14)

print(p_by_score_all)

####combined significant pVal <0.05 rnaseq and riboseq####
#correlagrams
library(ggpubr)  # para stat_cor()

p_by_score <- ggplot(top10_combined, aes(x = ribo_log2FC, y = rna_log2FC)) +
  geom_point(aes(color = gene_score), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.3) +
  facet_wrap(~ gene_score) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Correlation between Ribo-seq and RNA-seq by Gene Score",
    x = "Ribo-seq Log2 Fold Change",
    y = "RNA-seq Log2 Fold Change"
  ) +
  theme_minimal(base_size = 14)

print(p_by_score)

library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)  # para stat_cor()

p_by_score <- ggplot(top10_combined, aes(x = ribo_log2FC, y = rna_log2FC)) +
  geom_point(aes(color = gene_score, shape = localization), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "gray70", size = 0.7, alpha = 0.5) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.3) +
  facet_wrap(~ gene_score) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Correlation between Ribo-seq and RNA-seq by Gene Score",
    x = "Ribo-seq Log2 Fold Change",
    y = "RNA-seq Log2 Fold Change",
    shape = "Localization"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom"
  )

print(p_by_score)



# Filtrar genes significativos en ambos datasets
combined_sig_both <- ribo_all %>%
  left_join(rna, by = "gene", suffix = c("_ribo", "_rna")) %>%
  mutate(
    ribo_signif = ifelse(ribo_padj < 0.05, "Significant", "Not Significant"),
    rna_signif  = ifelse(rna_padj  < 0.05, "Significant", "Not Significant"),
    gene_score = factor(gene_score)
  ) %>%
  filter(ribo_padj < 0.05, rna_padj < 0.05)

# Correlation plot con solo genes significativos en ambos experimentos
library(ggplot2)
library(ggrepel)
library(ggpubr)

p_corr_sig <- ggplot(combined_sig_both, aes(x = ribo_log2FC, y = rna_log2FC)) +
  geom_point(aes(color = gene_score), size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  facet_wrap(~ gene_score) +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(
    title = "Correlation of Significant Genes (RNA-seq & Ribo-seq)",
    x = "Ribo-seq Log2 Fold Change",
    y = "RNA-seq Log2 Fold Change"
  ) +
  theme_minimal(base_size = 14)

print(p_corr_sig)




# Unir los datos RNA-seq y Ribo-seq, y filtrar por significancia
combined_sig_both <- ribo_all %>%
  left_join(rna, by = "gene", suffix = c("_ribo", "_rna")) %>%
  inner_join(meta %>% select(gene, syndromic), by = "gene") %>%
  filter(ribo_padj < 0.05, rna_padj < 0.05) %>%
  mutate(
    gene_score = factor(gene_score),
    syndromic = factor(syndromic, levels = c(1, 0), labels = c("Syndromic", "Non-Syndromic"))
  )

# Función para graficar por grupo syndromic
plot_correlation_by_syndromic <- function(data, label_title) {
  ggplot(data, aes(x = ribo_log2FC, y = rna_log2FC)) +
    geom_point(aes(color = gene_score, shape = localization), size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    stat_cor(
      aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
      method = "pearson",
      label.x.npc = "left",
      label.y.npc = "top"
    ) +
    facet_wrap(~ gene_score) +
    geom_text_repel(aes(label = gene), size = 3, max.overlaps = Inf, box.padding = 0.3) +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    labs(
      title = paste("Correlation of Significant Genes -", label_title),
      x = "Ribo-seq Log2 Fold Change",
      y = "RNA-seq Log2 Fold Change"
    ) +
    theme_minimal(base_size = 14)
}

# Dividir los datos
syndromic_data <- combined_sig_both %>% filter(syndromic == "Syndromic")
nonsyndromic_data <- combined_sig_both %>% filter(syndromic == "Non-Syndromic")

# Generar gráficos
p_syndromic <- plot_correlation_by_syndromic(syndromic_data, "Syndromic Genes")
p_nonsyndromic <- plot_correlation_by_syndromic(nonsyndromic_data, "Non-Syndromic Genes")

# Mostrar los gráficos
p_syndromic
p_nonsyndromic



ggplot(combined_all, aes(x = factor(gene_score), y = ribo_log2FC, color = localization)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_point(data = subset(combined_all, gene %in% c("rph3a", "caskin1")),
             aes(x = factor(gene_score), y = ribo_log2FC),
             color = "red", size = 5) +
  ggrepel::geom_text_repel(data = subset(combined_all, gene %in% c("rph3a", "caskin1")),
                           aes(label = gene),
                           color = "black", size = 4) +
  theme_classic() +
  labs(title = "Ribo-seq log2FC by Gene Score and Localization",
       x = "Gene Score",
       y = "Ribo-seq log2FC") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


library(ggplot2)
library(ggrepel)
library(dplyr)

# Filtrar solo npl-translation-up
plot_data <- combined_all %>%
  filter(localization == "npl-translation-up")

highlight_genes <- c("rph3a", "caskin1")

ggplot(plot_data, aes(x = factor(gene_score), y = ribo_log2FC)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  geom_point(data = subset(plot_data, gene %in% highlight_genes),
             aes(x = factor(gene_score), y = ribo_log2FC),
             color = "red", size = 4) +
  geom_text_repel(data = subset(plot_data, gene %in% highlight_genes),
                  aes(label = gene),
                  color = "red", size = 4) +
  theme_classic() +
  labs(title = "Ribo-seq log2FC by Gene Score (NPL Translation Up)",
       x = "SFARI Gene Score",
       y = "Ribo-seq log2FC") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

