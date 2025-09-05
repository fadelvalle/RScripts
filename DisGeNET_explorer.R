# Cargar librerías necesarias
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggprism)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)

# === Selección del archivo ===
cat("Selecciona el archivo de Ribo-seq (CSV con columnas: gene, ribo_log2FC, ribo_padj):\n")
ribo_file <- file.choose()

# === Cargar y procesar Ribo-seq ===
ribo <- read_csv(ribo_file)

# Asegurar nombres de genes en minúscula
ribo <- ribo %>% mutate(gene = toupper(gene))

# Filtrar genes significativos (padj < 0.05)
ribo_sig <- ribo %>% filter(ribo_padj < 0.05)

# Separar en positivos y negativos según log2FC
ribo_pos <- ribo_sig %>% filter(ribo_log2FC > 0)
ribo_neg <- ribo_sig %>% filter(ribo_log2FC < 0)

# Obtener IDs de Entrez para cada grupo
entrez_pos <- bitr(ribo_pos$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_neg <- bitr(ribo_neg$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Unir Entrez IDs al dataframe original
ribo_pos <- ribo_pos %>% inner_join(entrez_pos, by = c("gene" = "SYMBOL"))
ribo_neg <- ribo_neg %>% inner_join(entrez_neg, by = c("gene" = "SYMBOL"))

# === Enriquecimiento con DisGeNET ===
cat("Ejecutando enrichDGN para genes con log2FC positivo...\n")
edo_pos <- enrichDGN(ribo_pos$ENTREZID)

cat("Ejecutando enrichDGN para genes con log2FC negativo...\n")
edo_neg <- enrichDGN(ribo_neg$ENTREZID)

# === Visualizaciones ===

# Dotplots
dotplot(edo_pos, showCategory = 20) + ggtitle("DisGeNET enrichment - Positive log2FC")
dotplot(edo_neg, showCategory = 20) + ggtitle("DisGeNET enrichment - Negative log2FC")

# Barplots
barplot(edo_pos, showCategory = 20) + ggtitle("Barplot - Positive log2FC")
barplot(edo_neg, showCategory = 20) + ggtitle("Barplot - Negative log2FC")

# Network plots
edox_pos <- setReadable(edo_pos, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
edox_neg <- setReadable(edo_neg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

listFC_pos <- ribo_pos$ribo_log2FC
names(listFC_pos) <- ribo_pos$gene

listFC_neg <- ribo_neg$ribo_log2FC
names(listFC_neg) <- ribo_neg$gene

cnetplot(edox_pos, foldChange = listFC_pos) + ggtitle("Cnetplot - Positive")
cnetplot(edox_neg, foldChange = listFC_neg) + ggtitle("Cnetplot - Negative")

# Heatmaps & emap
edo2_pos <- pairwise_termsim(edox_pos)
edo2_neg <- pairwise_termsim(edox_neg)

emapplot(edo2_pos, layout = "kk") + ggtitle("emapplot - Positive")
emapplot(edo2_neg, layout = "kk") + ggtitle("emapplot - Negative")

heatplot(edo2_pos, foldChange = listFC_pos)
heatplot(edo2_neg, foldChange = listFC_neg)

# === Exportar resultados ===
write_csv(as.data.frame(edo_pos), "disgenet_enrichment_positive.csv")
write_csv(as.data.frame(edo_neg), "disgenet_enrichment_negative.csv")

cat("Listo. Resultados guardados como:\n")
cat("- disgenet_enrichment_positive.csv\n")
cat("- disgenet_enrichment_negative.csv\n")

# === Coincidencia con genes SFARI ===
cat("Selecciona el archivo de anotación genética SFARI (con columna 'gene-symbol'):\n")
sfari_file <- file.choose()
sfari <- read_csv(sfari_file)

# Renombrar y asegurar que los nombres de los genes estén en mayúsculas
sfari <- sfari %>%
  rename(gene = `gene-symbol`) %>%
  mutate(gene = toupper(gene))

# Seleccionar solo columnas necesarias
sfari_clean <- sfari %>%
  select(gene, `gene-score`, syndromic, `genetic-category`) %>%
  rename(
    gene_score = `gene-score`,
    category = `genetic-category`
  )

# Filtrar coincidencias entre genes significativos en Ribo-seq y SFARI
ribo_sfari_overlap <- ribo_sig %>%
  filter(gene %in% sfari_clean$gene) %>%
  mutate(gene = toupper(gene)) %>%  # Asegurar formato consistente antes del join
  left_join(sfari_clean, by = "gene")

cat("Número de genes Ribo-seq significativos también en SFARI:", nrow(ribo_sfari_overlap), "\n")

# === Visualización: log2FC vs gene-score ===
if(nrow(ribo_sfari_overlap) > 0){
  ggplot(ribo_sfari_overlap, aes(x = gene_score, y = ribo_log2FC, color = category)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = gene), max.overlaps = 15, size = 3.5) +
    labs(
      title = "Genes Ribo-seq anotados en SFARI",
      x = "Gene Score (SFARI)",
      y = "log2FC (Ribo-seq)"
    ) +
    theme_prism(base_size = 12) +
    scale_color_brewer(palette = "Set1")
}


ggplot(ribo_sfari_overlap, aes(x = category, y = ribo_log2FC, fill = category)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(
    title = "Distribución de log2FC por categoría genética (SFARI)",
    x = "Categoría genética (SFARI)",
    y = "log2FC (Ribo-seq)"
  ) +
  theme_prism(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1")


# === Enriquecimiento DisGeNET para genes Ribo-seq anotados en SFARI ===
cat("Ejecutando enrichDGN para genes Ribo-seq anotados en SFARI...\n")

# Obtener IDs Entrez para genes anotados en SFARI
sfari_entrez <- bitr(ribo_sfari_overlap$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Unir con ribo_sfari_overlap para obtener log2FC
ribo_sfari_entrez <- ribo_sfari_overlap %>%
  inner_join(sfari_entrez, by = c("gene" = "SYMBOL"))

# Enriquecimiento
edo_sfari <- enrichDGN(ribo_sfari_entrez$ENTREZID)

# Visualización: Dotplot
dotplot(edo_sfari, showCategory = 20) + ggtitle("DisGeNET enrichment - SFARI genes")

# Visualización: Barplot
barplot(edo_sfari, showCategory = 20) + ggtitle("Barplot - SFARI genes")

# Visualización: Cnetplot
edox_sfari <- setReadable(edo_sfari, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
listFC_sfari <- ribo_sfari_entrez$ribo_log2FC
names(listFC_sfari) <- ribo_sfari_entrez$gene

cnetplot(edox_sfari, foldChange = listFC_sfari) + ggtitle("Cnetplot - SFARI genes")

# Heatmap y emapplot
edo2_sfari <- pairwise_termsim(edox_sfari)

emapplot(edo2_sfari, layout = "kk") + ggtitle("emapplot - SFARI genes")
heatplot(edo2_sfari, foldChange = listFC_sfari)

# Exportar resultados
write_csv(as.data.frame(edo_sfari), "disgenet_enrichment_sfari.csv")
cat("Resultados de enriquecimiento para genes SFARI guardados como: disgenet_enrichment_sfari.csv\n")

# === Enriquecimiento DisGeNET por log2FC en genes SFARI ===
cat("Procesando enriquecimiento DisGeNET para genes SFARI con log2FC positivo y negativo...\n")

# Filtrar por log2FC positivo y negativo
ribo_sfari_pos <- ribo_sfari_overlap %>% filter(ribo_log2FC > 0)
ribo_sfari_neg <- ribo_sfari_overlap %>% filter(ribo_log2FC < 0)

# Convertir a ENTREZID
entrez_sfari_pos <- bitr(ribo_sfari_pos$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_sfari_neg <- bitr(ribo_sfari_neg$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Unir para tener log2FC y anotaciones
ribo_sfari_pos <- ribo_sfari_pos %>% inner_join(entrez_sfari_pos, by = c("gene" = "SYMBOL"))
ribo_sfari_neg <- ribo_sfari_neg %>% inner_join(entrez_sfari_neg, by = c("gene" = "SYMBOL"))

# Enriquecimiento
edo_sfari_pos <- enrichDGN(ribo_sfari_pos$ENTREZID)
edo_sfari_neg <- enrichDGN(ribo_sfari_neg$ENTREZID)

# === Visualización ===

# Dotplots
dotplot(edo_sfari_pos, showCategory = 20) + ggtitle("DisGeNET - SFARI Positive log2FC")
dotplot(edo_sfari_neg, showCategory = 20) + ggtitle("DisGeNET - SFARI Negative log2FC")


# Barplots
barplot(edo_sfari_pos, showCategory = 20) + ggtitle("Barplot - SFARI Positive log2FC")
barplot(edo_sfari_neg, showCategory = 20) + ggtitle("Barplot - SFARI Negative log2FC")

# Cnetplots
edox_sfari_pos <- setReadable(edo_sfari_pos, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
edox_sfari_neg <- setReadable(edo_sfari_neg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

listFC_sfari_pos <- ribo_sfari_pos$ribo_log2FC
names(listFC_sfari_pos) <- ribo_sfari_pos$gene

listFC_sfari_neg <- ribo_sfari_neg$ribo_log2FC
names(listFC_sfari_neg) <- ribo_sfari_neg$gene

cnetplot(edox_sfari_pos, foldChange = listFC_sfari_pos) + ggtitle("Cnetplot - SFARI Positive")
cnetplot(edox_sfari_neg, foldChange = listFC_sfari_neg) + ggtitle("Cnetplot - SFARI Negative")

# Emapplot & Heatplot
edo2_sfari_pos <- pairwise_termsim(edox_sfari_pos)
edo2_sfari_neg <- pairwise_termsim(edox_sfari_neg)

emapplot(edo2_sfari_pos, layout = "kk") + ggtitle("emapplot - SFARI Positive")
emapplot(edo2_sfari_neg, layout = "kk") + ggtitle("emapplot - SFARI Negative")

heatplot(edo2_sfari_pos, foldChange = listFC_sfari_pos) + ggtitle("Heatplot - SFARI Positive")
heatplot(edo2_sfari_neg, foldChange = listFC_sfari_neg) + ggtitle("Heatplot - SFARI Negative")

# Exportar resultados
write_csv(as.data.frame(edo_sfari_pos), "disgenet_enrichment_sfari_positive.csv")
write_csv(as.data.frame(edo_sfari_neg), "disgenet_enrichment_sfari_negative.csv")

cat("Resultados guardados como:\n")
cat("- disgenet_enrichment_sfari_positive.csv\n")
cat("- disgenet_enrichment_sfari_negative.csv\n")



# === Integración con ClinGen haploinsuficiencia y triplosensibilidad desde Clingen1.csv ===
cat("Cargando archivo Clingen1.csv...\n")
clingen <- read_csv("Clingen1.csv", show_col_types = FALSE)
# Estandarizar nombres de columnas y contenido
clingen_clean <- clingen %>%
  rename(
    gene = `GENE/REGION`,
    hi_score = `HAPLOINSUFFICIENCY`,
    ts_score = `TRIPLOSENSITIVITY`,
    online_report = `ONLINE REPORT`
  ) %>%
  mutate(gene = toupper(gene)) %>%
  select(gene, hi_score, ts_score, online_report)
# Asegurar genes de ribo_sig estén en mayúscula
ribo_sig <- ribo_sig %>% mutate(gene = toupper(gene))

# Cruce con genes significativos del Ribo-seq
ribo_clingen_overlap <- ribo_sig %>%
  filter(gene %in% clingen_clean$gene) %>%
  left_join(clingen_clean, by = "gene")

cat("Genes significativos de Ribo-seq con anotaciones ClinGen:", nrow(ribo_clingen_overlap), "\n")

# === Visualización ===
if(nrow(ribo_clingen_overlap) > 0){
  ggplot(ribo_clingen_overlap, aes(x = hi_score, y = ribo_log2FC)) +
    geom_jitter(width = 0.1, size = 2.5, alpha = 0.7) +
    geom_text_repel(aes(label = gene), max.overlaps = 15, size = 3.5) +
    labs(
      title = "log2FC de genes Ribo-seq con anotaciones de haploinsuficiencia (ClinGen)",
      x = "HI Score (ClinGen)",
      y = "log2FC (Ribo-seq)"
    ) +
    theme_prism(base_size = 12)
}

# Exportar resultados
write_csv(ribo_clingen_overlap, "ribo_clingen_overlap.csv")
cat("Resultados de genes anotados en ClinGen guardados como: ribo_clingen_overlap.csv\n")


library(ggplot2)
library(dplyr)

# Aseguramos que la variable localization esté presente y sin NA
ribo_hi_plot <- ribo_clingen_overlap %>%
  filter(!is.na(hi_score) & hi_score != "" & !is.na(localization)) %>%
  mutate(
    gene = factor(gene),
    log2fc_cat = cut(ribo_log2FC, breaks = c(-Inf, -2, -1.3, 0, 1.3, 2, Inf),
                     labels = c("< -2", "-2 to -1.3", "-1.3 to 0", "0 to 1.3", "1.3 to 2", "> 2"))
  )

hi_categories <- unique(ribo_hi_plot$hi_score)

for (cat in hi_categories) {
  sub_df <- ribo_hi_plot %>%
    filter(hi_score == cat) %>%
    arrange(desc(ribo_log2FC)) %>%
    group_by(localization) %>%
    mutate(
      index = row_number(),
      col = (index - 1) %% 3 + 1,
      row = (index - 1) %/% 3 + 1
    ) %>%
    ungroup()
  
  p <- ggplot(sub_df, aes(x = col, y = -row, fill = ribo_log2FC)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(min(ribo_hi_plot$ribo_log2FC, na.rm = TRUE),
                 max(ribo_hi_plot$ribo_log2FC, na.rm = TRUE)),
      name = "log2FC"
    ) +
    geom_text(aes(label = gene), size = 3) +
    facet_wrap(~localization) +
    labs(
      title = paste("HI Score:", cat),
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right",
      strip.text = element_text(size = 12)
    )
  
  print(p)
}


# === Integración combinada SFARI + ClinGen ===

# Unir ribo_sfari_overlap con clingen_clean por 'gene'
ribo_sfari_clingen <- ribo_sfari_overlap %>%
  left_join(clingen_clean, by = "gene")

cat("Número de genes en el cruce combinado SFARI + ClinGen:", nrow(ribo_sfari_clingen), "\n")

# Exportar resultados combinados
write_csv(ribo_sfari_clingen, "ribo_sfari_clingen_combined.csv")
cat("Archivo combinado SFARI + ClinGen guardado como: ribo_sfari_clingen_combined.csv\n")

# Visualización simple: log2FC vs gene_score (SFARI) coloreado por haploinsuficiencia (ClinGen)
if(nrow(ribo_sfari_clingen) > 0){
  ggplot(ribo_sfari_clingen, aes(x = gene_score, y = ribo_log2FC, color = hi_score)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = gene), max.overlaps = 30, size = 3.5) +
    labs(
      title = "Genes Ribo-seq con anotaciones combinadas SFARI + ClinGen",
      x = "Gene Score (SFARI)",
      y = "log2FC (Ribo-seq)",
      color = "HI Score (ClinGen)"
    ) +
    theme_prism(base_size = 12) +
    scale_color_brewer(palette = "Dark2")
}

# Puedes seguir agregando análisis o visualizaciones con esta tabla combinada:
# ribo_sfari_clingen contiene columnas:
# gene, ribo_baseMean, ribo_log2FC, ribo_padj, localization, gene_score, syndromic, category (SFARI),
# hi_score, ts_score, online_report (ClinGen)

