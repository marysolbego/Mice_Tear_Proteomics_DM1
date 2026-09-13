
# Proteomics Differential Expression Analysis (DEA)

#Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(matrixStats)
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  library(ggrepel)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(readxl)
  library(ggvenn)
  library(patchwork)
  library(vegan)
  library(ggpubr)
})


#Define sample sheet file
metadata_file <- 'Raw_data/sample_sheet.csv'

#Define Raw data files
raw_files <- list.files(path = 'Raw_data', pattern = "(?i)_total_proteins\\.csv$", full.names = TRUE)


OUTPUT_DIR <- file.path(getwd(), "results")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

# Statistical Parameters 
pvalue_threshold <- 0.05
ratio_threshold <- log2(1.5)
peptide_count_threshold <- 1
unique_peptide_threshold <- 1
n_best <- 5

# Load Sample Sheet
the_samples <- read_csv(metadata_file, show_col_types = FALSE) %>% 
  mutate(muestra = make_clean_names(muestra))

calc_condition_stats <- function(quant_df, samples_df, cond_name) {
  ids <- samples_df %>% filter(condition == cond_name) %>% pull(muestra)
  if(length(ids) == 0 || all(ids == "")) stop(paste("Condition", cond_name, "not found in metadata."))
  
  regex <- paste(ids, collapse = "|")
  pos_samples <- str_detect(colnames(quant_df), regex) %>% which()
  if(length(pos_samples) > length(ids)) pos_samples <- pos_samples[1:length(ids)]
  
  df <- quant_df %>% select(all_of(pos_samples)) %>%
    mutate(across(everything(), ~ as.numeric(str_replace(as.character(.), ",", "."))))
  
  mat <- as.matrix(df)
  stats_tibble <- tibble(
    average = rowMeans(mat, na.rm = TRUE),
    log10_avg = log10(rowMeans(mat, na.rm = TRUE) + 1),
    sd = rowSds(mat, na.rm = TRUE),
    cv = rowSds(mat, na.rm = TRUE) / rowMeans(mat, na.rm = TRUE),
    count = rowSums(mat > 0, na.rm = TRUE)
  )
  
  stats_tibble <- stats_tibble %>% rename_with(~ paste0(., ".", cond_name))
  return(stats_tibble)
}

format_pub_label <- function(label) {
  label <- str_replace_all(label, "(?i)Control", "Non-diabetic")
  label <- str_replace_all(label, "(?i)STZ", "Diabetic")
  label <- str_replace_all(label, "(?i)_meses", " months")
  label <- str_replace_all(label, "(?i)_mes", " month")
  label <- str_replace_all(label, "_", " ") 
  return(trimws(label))
}

for (current_file in raw_files) {
  file_name <- basename(current_file)
  
  clean_name <- str_remove_all(file_name, "(?i)_total_proteins\\.csv")
  name_parts <- str_split(clean_name, "_vs_")[[1]]
  cond_1 <- name_parts[1] 
  cond_2 <- name_parts[2] 
  
  raw_data <- read_csv(current_file, skip = 2, show_col_types = FALSE) %>% clean_names() 
  raw_data <- raw_data[, 1:(ncol(raw_data)-6)]
  
  quantified <- raw_data %>%
    filter(unique_peptides != 0,
           max_fold_change != Inf,
           !str_detect(toupper(accession), "REVERSE"))
  
  stats_cond_1 <- calc_condition_stats(quantified, the_samples, cond_1)
  stats_cond_2 <- calc_condition_stats(quantified, the_samples, cond_2)
  
  base_data <- quantified %>%
    select(accession, peptide_count, unique_peptides, confidence_score, anova_p, q_value, description) %>%
    mutate(
      gene_name = str_extract(description, "GN=\\S+"),
      gene_name = str_remove(gene_name, "GN="),
      gene_name = ifelse(is.na(gene_name), accession, gene_name)
    )
  
  processed_data <- bind_cols(base_data, stats_cond_1, stats_cond_2)
  
  avg_c1_col <- paste0("average.", cond_1)
  avg_c2_col <- paste0("average.", cond_2)
  count_c1_col <- paste0("count.", cond_1)
  count_c2_col <- paste0("count.", cond_2)
  
  final_data <- processed_data %>%
    mutate(
      nlog10_anova_p = -log10(anova_p),
      ratio = !!sym(avg_c2_col) / !!sym(avg_c1_col),
      log2_ratio = log2(ratio),
      significance = case_when(
        log2_ratio > ratio_threshold & anova_p < pvalue_threshold ~ "UP",
        log2_ratio < -ratio_threshold & anova_p < pvalue_threshold ~ "DOWN",
        TRUE ~ "Not significant"
      )
    )
  
  final_data$significance <- factor(final_data$significance, levels = c("UP", "DOWN", "Not significant"))
  
  filtered_data <- final_data %>%
    filter(
      !!sym(count_c1_col) >= 1, 
      !!sym(count_c2_col) >= 1,
      peptide_count >= peptide_count_threshold,
      unique_peptides >= unique_peptide_threshold
    )
  
  peptides_up <- filtered_data %>% filter(significance == "UP")
  peptides_down <- filtered_data %>% filter(significance == "DOWN")
  
  n_total  <- nrow(filtered_data)
  n_up     <- nrow(peptides_up)
  n_down   <- nrow(peptides_down)
  
  base_output_name <- str_remove(file_name, "(?i)\\.csv$")
  outputfile1 <- file.path(OUTPUT_DIR, paste0(base_output_name, ".volcano.png"))
  outputfile2 <- file.path(OUTPUT_DIR, paste0(base_output_name, ".named.png"))
  outputfile3 <- file.path(OUTPUT_DIR, paste0(base_output_name, ".Publication_Tables.xlsx"))
  csv_final_path <- file.path(OUTPUT_DIR, paste0(base_output_name, "_Supp_All.csv"))
  csv_top_path <- file.path(OUTPUT_DIR, paste0(base_output_name, "_Table1_Top5.csv"))
  
  #Volcano Plot
  y_max <- quantile(filtered_data$nlog10_anova_p, 0.99, na.rm = TRUE)
  
  volcano1 <- ggplot(filtered_data, aes(x = log2_ratio, y = nlog10_anova_p, color = significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = c(-ratio_threshold, ratio_threshold), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", linewidth = 0.4, color = "black") +
    scale_color_manual(values = c("UP" = "#D55E00", "DOWN" = "#0072B2", "Not significant" = "#CCCCCC")) +
    scale_y_continuous(limits = c(0, y_max * 1.15), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(x = expression(log[2]~Fold~Change), y = expression(-log[10](italic(p))), color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = c(0.88, 0.80), legend.background = element_blank(), legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(30, 30, 20, 30)
    ) +
    annotate("text", x = min(filtered_data$log2_ratio, na.rm = TRUE), y = y_max * 1.12, hjust = 0, label = paste0("Total = ", n_total), size = 4, fontface = "bold") +
    annotate("text", x = min(filtered_data$log2_ratio, na.rm = TRUE), y = y_max * 1.05, hjust = 0, label = paste0("UP = ", n_up, "   |   DOWN = ", n_down), size = 3.6)
  
  ggsave(filename = outputfile1, plot = volcano1, width = 10, height = 7, dpi = 600)
  
  top_up_plot <- peptides_up %>% arrange(desc(log2_ratio)) %>% slice_head(n = n_best)
  top_down_plot <- peptides_down %>% arrange(log2_ratio) %>% slice_head(n = n_best)
  
  volcano2 <- volcano1 +
    geom_label_repel(data = top_up_plot, aes(label = gene_name), size = 3, box.padding = 0.3, point.padding = 0.2, show.legend = FALSE) +
    geom_label_repel(data = top_down_plot, aes(label = gene_name), size = 3, box.padding = 0.3, point.padding = 0.2, show.legend = FALSE)
  
  ggsave(filename = outputfile2, plot = volcano2, width = 10, height = 7, dpi = 600)
  
  # Tables
  up_pub <- peptides_up %>%
    mutate(
      `Protein Description` = str_trim(str_remove(description, " OS=.*")),
      `Log2 Fold Change`    = round(log2_ratio, 2),
      `Linear Fold Change`  = round(ifelse(ratio >= 1, ratio, -1/ratio), 2),
      `p-value`             = signif(anova_p, 3)
    ) %>%
    select(Accession = accession, `Gene Name` = gene_name, `Protein Description`, 
           `Log2 Fold Change`, `Linear Fold Change`, `p-value`) %>%
    arrange(desc(`Log2 Fold Change`))
  
  down_pub <- peptides_down %>%
    mutate(
      `Protein Description` = str_trim(str_remove(description, " OS=.*")),
      `Log2 Fold Change`    = round(log2_ratio, 2),
      `Linear Fold Change`  = round(ifelse(ratio >= 1, ratio, -1/ratio), 2),
      `p-value`             = signif(anova_p, 3)
    ) %>%
    select(Accession = accession, `Gene Name` = gene_name, `Protein Description`, 
           `Log2 Fold Change`, `Linear Fold Change`, `p-value`) %>%
    arrange(`Log2 Fold Change`)
  
  top5_up <- up_pub %>% slice_head(n = n_best) %>% mutate(Regulation = "UP-REGULATED")
  top5_down <- down_pub %>% slice_head(n = n_best) %>% mutate(Regulation = "DOWN-REGULATED")
  
  table1_main <- bind_rows(top5_up, top5_down) %>%
    select(Regulation, Accession, `Gene Name`, `Protein Description`, 
           `Fold Change` = `Linear Fold Change`, `p-value`)
  
  wb <- createWorkbook()
  addWorksheet(wb, "Table_1_Top_Hits")
  writeData(wb, "Table_1_Top_Hits", table1_main)
  addWorksheet(wb, "Supp_All_UP")
  writeData(wb, "Supp_All_UP", up_pub)
  addWorksheet(wb, "Supp_All_DOWN")
  writeData(wb, "Supp_All_DOWN", down_pub)
  
  saveWorkbook(wb, outputfile3, overwrite = TRUE)
  write.csv(table1_main, file = csv_top_path, row.names = FALSE)
  
  all_proteins_pub <- bind_rows(
    up_pub %>% mutate(Regulation = "UP"), 
    down_pub %>% mutate(Regulation = "DOWN")
  )
  write.csv(all_proteins_pub, file = csv_final_path, row.names = FALSE)
}

#Heatmap Visualization 

excel_files <- list.files(path = OUTPUT_DIR, pattern = "(?i)\\.Publication_Tables\\.xlsx$", full.names = TRUE)

for (current_excel in excel_files) {
  base_name <- basename(current_excel) %>% str_remove("(?i)\\.Publication_Tables\\.xlsx$")
  
  clean_name <- str_remove_all(base_name, "(?i)_total_proteins")
  name_parts <- str_split(clean_name, "_vs_")[[1]]
  
  if (length(name_parts) < 2) next
  
  cond_1 <- name_parts[1]
  cond_2 <- name_parts[2]
  
  lbl_cond1 <- format_pub_label(cond_1)
  lbl_cond2 <- format_pub_label(cond_2)
  
  raw_file_pattern <- paste0("(?i)", cond_1, ".*_vs_.*", cond_2, ".*\\.csv$")
  matching_raw_files <- grep(raw_file_pattern, raw_files, value = TRUE)
  
  if (length(matching_raw_files) == 0) next
  raw_file <- matching_raw_files[1] 
  
  sheet_up <- if("Supp_All_UP" %in% excel_sheets(current_excel)) "Supp_All_UP" else "UP_Proteins"
  sheet_down <- if("Supp_All_DOWN" %in% excel_sheets(current_excel)) "Supp_All_DOWN" else "DOWN_Proteins"
  
  up_data <- suppressMessages(read_excel(current_excel, sheet = sheet_up))
  down_data <- suppressMessages(read_excel(current_excel, sheet = sheet_down))
  
  up_accessions <- up_data$Accession
  down_accessions <- down_data$Accession
  all_hits <- c(up_accessions, down_accessions)
  
  if (length(all_hits) < 2) next
  
  dict_genes <- bind_rows(up_data, down_data) %>%
    select(Accession, Gene_Name = `Gene Name`)
  
  data_total <- read_csv(raw_file, skip = 2, show_col_types = FALSE) %>% clean_names()
  col_accession <- grep("accession", colnames(data_total), ignore.case = TRUE)[1]
  
  regex_cond1 <- paste0("(?i)", str_replace_all(cond_1, "_", ".*"))
  regex_cond2 <- paste0("(?i)", str_replace_all(cond_2, "_", ".*"))
  
  cols_cond1 <- grep(regex_cond1, colnames(data_total))
  cols_cond2 <- grep(regex_cond2, colnames(data_total))
  
  if (length(cols_cond1) > 3) cols_cond1 <- cols_cond1[1:3]
  if (length(cols_cond2) > 3) cols_cond2 <- cols_cond2[1:3]
  
  if (length(cols_cond1) == 0 || length(cols_cond2) == 0) next
  
  expression_matrix <- data_total %>%
    select(Accession = all_of(col_accession), all_of(cols_cond1), all_of(cols_cond2)) %>%
    filter(Accession %in% all_hits) %>%
    filter(Accession != "" & !is.na(Accession)) %>%
    distinct(Accession, .keep_all = TRUE) %>% 
    column_to_rownames("Accession") %>%
    mutate(across(everything(), ~as.numeric(gsub(",", ".", .)))) %>%
    drop_na()
  
  log_matrix <- log2(expression_matrix + 1)
  
  regulation_vector <- factor(
    ifelse(rownames(log_matrix) %in% up_accessions, "UP", "DOWN"), 
    levels = c("UP", "DOWN")
  )
  
  zscore_matrix <- t(scale(t(log_matrix)))
  nombres_genes <- dict_genes$Gene_Name[match(rownames(log_matrix), dict_genes$Accession)]
  nombres_genes <- ifelse(is.na(nombres_genes) | nombres_genes == "", rownames(log_matrix), nombres_genes)
  rownames(zscore_matrix) <- nombres_genes
  
  display_colnames <- c(rep(lbl_cond1, length(cols_cond1)), rep(lbl_cond2, length(cols_cond2)))
  display_colnames <- paste0(display_colnames, " N", c(1:length(cols_cond1), 1:length(cols_cond2)))
  colnames(zscore_matrix) <- display_colnames
  
  ha <- HeatmapAnnotation(
    Group = factor(c(rep(lbl_cond1, length(cols_cond1)), rep(lbl_cond2, length(cols_cond2))), levels = c(lbl_cond1, lbl_cond2)),
    col = list(Group = setNames(c("#56B4E9", "#CC79A7"), c(lbl_cond1, lbl_cond2)))
  )
  
  ra <- rowAnnotation(
    Regulation = regulation_vector,
    col = list(Regulation = c("UP" = "#E69F00", "DOWN" = "#0072B2"))
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#D55E00"))
  
  ht <- Heatmap(
    zscore_matrix,
    name = "Z-score",
    col = col_fun,
    top_annotation = ha,
    left_annotation = ra,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    row_split = regulation_vector,
    row_title = NULL,                                               
    column_title = paste("Significant Hits Tear Proteome (", lbl_cond1, "vs", lbl_cond2, ")"), 
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = gpar(fontsize = 7),                              
    column_names_gp = gpar(fontsize = 10)                           
  )
  
  output_pdf <- file.path(OUTPUT_DIR, paste0(base_name, "_Heatmap.pdf"))
  output_png <- file.path(OUTPUT_DIR, paste0(base_name, "_Heatmap.png"))
  output_tiff <- file.path(OUTPUT_DIR, paste0(base_name, "_Heatmap.tiff"))
  
  pdf(output_pdf, width = 8, height = 12)
  draw(ht, merge_legend = TRUE)
  dev.off()
  
  png(output_png, width = 8, height = 12, units = "in", res = 600)
  draw(ht, merge_legend = TRUE)
  dev.off()
  
  tiff(output_tiff, width = 8, height = 12, units = "in", res = 600, compression = "lzw")
  draw(ht, merge_legend = TRUE)
  dev.off()
}

#Venn Diagrams & Biological Signatures Extraction


file_ctrl <- file.path(OUTPUT_DIR, "Control_1_mes_vs_Control_3_meses_Total_proteins.Publication_Tables.xlsx")
file_stz  <- file.path(OUTPUT_DIR, "STZ_1_mes_vs_STZ_3_meses_Total_proteins.Publication_Tables.xlsx")

if (!file.exists(file_ctrl)) stop("Control file not found in results folder.")
if (!file.exists(file_stz)) stop("STZ file not found in results folder.")

#Protein lists
ctrl_up_df   <- read_excel(file_ctrl, sheet = "Supp_All_UP")
ctrl_down_df <- read_excel(file_ctrl, sheet = "Supp_All_DOWN")
stz_up_df    <- read_excel(file_stz, sheet = "Supp_All_UP")
stz_down_df  <- read_excel(file_stz, sheet = "Supp_All_DOWN")

#Extract Accessions
ctrl_up   <- ctrl_up_df$Accession
ctrl_down <- ctrl_down_df$Accession
stz_up    <- stz_up_df$Accession
stz_down  <- stz_down_df$Accession

list_up <- list(
  "Nondiabetic" = ctrl_up,
  "Diabetic" = stz_up
)

list_down <- list(
  "Nondiabetic" = ctrl_down,
  "Diabetic" = stz_down
)

cb_palette <- c("#56B4E9", "#E69F00")

p_up <- ggvenn(
  list_up, 
  fill_color = cb_palette,
  stroke_size = 0.5,        
  set_name_size = 5,        
  text_size = 6,            
  show_percentage = FALSE   
)

p_down <- ggvenn(
  list_down, 
  fill_color = cb_palette,
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 6,
  show_percentage = FALSE
)

final_plot <- p_up + p_down

ruta_pdf  <- file.path(OUTPUT_DIR, "Venn_Diagrams.pdf")
ruta_png  <- file.path(OUTPUT_DIR, "Venn_Diagrams.png")
ruta_tiff <- file.path(OUTPUT_DIR, "Venn_Diagrams.tiff")

ggsave(ruta_pdf,  final_plot, width = 12, height = 6, dpi = 300)
ggsave(ruta_png,  final_plot, width = 12, height = 6, dpi = 300)
ggsave(ruta_tiff, final_plot, width = 12, height = 6, dpi = 300, compression = "lzw")


exclusive_diabetic_up   <- setdiff(stz_up, ctrl_up)
exclusive_diabetic_down <- setdiff(stz_down, ctrl_down)

exclusive_basal_up   <- setdiff(ctrl_up, stz_up)
exclusive_basal_down <- setdiff(ctrl_down, stz_down)

intersection_aging_up   <- intersect(ctrl_up, stz_up)
intersection_aging_down <- intersect(ctrl_down, stz_down)

stz_up_exclusive_df   <- stz_up_df %>% filter(Accession %in% exclusive_diabetic_up)
stz_down_exclusive_df <- stz_down_df %>% filter(Accession %in% exclusive_diabetic_down)

ctrl_up_exclusive_df   <- ctrl_up_df %>% filter(Accession %in% exclusive_basal_up)
ctrl_down_exclusive_df <- ctrl_down_df %>% filter(Accession %in% exclusive_basal_down)

aging_up_df   <- ctrl_up_df %>% filter(Accession %in% intersection_aging_up)
aging_down_df <- ctrl_down_df %>% filter(Accession %in% intersection_aging_down)

wb <- createWorkbook()

addWorksheet(wb, "UP_Diabetic_Exclusive")
writeData(wb, "UP_Diabetic_Exclusive", stz_up_exclusive_df)
addWorksheet(wb, "DOWN_Diabetic_Exclusive")
writeData(wb, "DOWN_Diabetic_Exclusive", stz_down_exclusive_df)

addWorksheet(wb, "UP_Basal_Exclusive")
writeData(wb, "UP_Basal_Exclusive", ctrl_up_exclusive_df)
addWorksheet(wb, "DOWN_Basal_Exclusive")
writeData(wb, "DOWN_Basal_Exclusive", ctrl_down_exclusive_df)

addWorksheet(wb, "UP_Aging_Intersection")
writeData(wb, "UP_Aging_Intersection", aging_up_df)
addWorksheet(wb, "DOWN_Aging_Intersection")
writeData(wb, "DOWN_Aging_Intersection", aging_down_df)

ruta_excel_venn <- file.path(OUTPUT_DIR, "Venn_Filtered_Proteins_Lists.xlsx")
saveWorkbook(wb, ruta_excel_venn, overwrite = TRUE)

#Divergent Lollipop Plot for Diabetic Signature

archivo_venn <- file.path(OUTPUT_DIR, "Venn_Filtered_Proteins_Lists.xlsx")
if (!file.exists(archivo_venn)) stop("Venn_Filtered_Proteins_Lists.xlsx not found.")

stz_up <- read_excel(archivo_venn, sheet = "UP_Diabetic_Exclusive") %>% 
  mutate(Direction = "Up-regulated")

stz_down <- read_excel(archivo_venn, sheet = "DOWN_Diabetic_Exclusive") %>% 
  mutate(Direction = "Down-regulated")

lollipop_data <- bind_rows(stz_up, stz_down) %>%
  select(`Gene Name`, `Log2 Fold Change`, Direction) %>%
  arrange(`Log2 Fold Change`) %>%
  mutate(`Gene Name` = factor(`Gene Name`, levels = unique(`Gene Name`)))

plot_stz <- ggplot(lollipop_data, aes(x = `Gene Name`, y = `Log2 Fold Change`, color = Direction)) +
  geom_segment(aes(x = `Gene Name`, xend = `Gene Name`, y = 0, yend = `Log2 Fold Change`), 
               color = "grey80", linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Down-regulated" = "#56B4E9", "Up-regulated" = "#E69F00")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9, color = "grey20"),
    axis.text.x = element_text(size = 10, color = "grey20"),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x = "Protein (Gene Name)",
    y = expression(log[2]~Fold~Change)
  )

ruta_lollipop_pdf  <- file.path(OUTPUT_DIR, "Lollipop_Diabetic_Signature.pdf")
ruta_lollipop_png  <- file.path(OUTPUT_DIR, "Lollipop_Diabetic_Signature.png")
ruta_lollipop_tiff <- file.path(OUTPUT_DIR, "Lollipop_Diabetic_Signature.tiff")

ggsave(ruta_lollipop_pdf,  plot_stz, width = 7, height = 9, dpi = 300)
ggsave(ruta_lollipop_png,  plot_stz, width = 7, height = 9, dpi = 300)
ggsave(ruta_lollipop_tiff, plot_stz, width = 7, height = 9, dpi = 300, compression = "lzw")


#Principal Component Analysis (PCA)

for (current_file in raw_files) {
  file_name <- basename(current_file)
  base_name <- str_remove(file_name, "(?i)\\.csv$")
  
  clean_name <- str_remove_all(base_name, "(?i)_total_proteins")
  name_parts <- str_split(clean_name, "_vs_")[[1]]
  
  if(length(name_parts) < 2) next
  
  cond_1 <- name_parts[1] 
  cond_2 <- name_parts[2] 
  
  lbl_cond1 <- format_pub_label(cond_1)
  lbl_cond2 <- format_pub_label(cond_2)
  
  ids_c1 <- the_samples %>% filter(condition == cond_1) %>% pull(muestra)
  ids_c2 <- the_samples %>% filter(condition == cond_2) %>% pull(muestra)
  
  if(length(ids_c1) == 0 || length(ids_c2) == 0) next
  
  raw_data <- read_csv(current_file, skip = 2, show_col_types = FALSE) %>% clean_names()
  
  regex_c1 <- paste(ids_c1, collapse = "|")
  regex_c2 <- paste(ids_c2, collapse = "|")
  
  pos_c1 <- which(str_detect(colnames(raw_data), regex_c1))
  pos_c2 <- which(str_detect(colnames(raw_data), regex_c2))
  
  if(length(pos_c1) > length(ids_c1)) pos_c1 <- pos_c1[1:length(ids_c1)]
  if(length(pos_c2) > length(ids_c2)) pos_c2 <- pos_c2[1:length(ids_c2)]
  
  if(length(pos_c1) == 0 || length(pos_c2) == 0) next
  
  expr_data <- raw_data %>%
    select(accession, all_of(pos_c1), all_of(pos_c2)) %>%
    filter(!is.na(accession) & accession != "") %>%
    distinct(accession, .keep_all = TRUE) %>%
    column_to_rownames("accession") %>%
    mutate(across(everything(), ~ as.numeric(str_replace(as.character(.), ",", "."))))
  
  expr_matrix <- as.matrix(expr_data)
  expr_matrix[is.na(expr_matrix)] <- 0
  log_matrix <- log2(expr_matrix + 1)
  
  variances <- apply(log_matrix, 1, var)
  log_matrix <- log_matrix[variances > 0, ]
  
  pca_input <- t(log_matrix)
  pca_result <- prcomp(pca_input, scale. = TRUE)
  
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
  pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")
  
  pca_df <- as.data.frame(pca_result$x) %>%
    rownames_to_column("SampleID") %>%
    mutate(
      Condition = case_when(
        str_detect(SampleID, regex_c1) ~ lbl_cond1,
        str_detect(SampleID, regex_c2) ~ lbl_cond2,
        TRUE ~ "Unknown"
      )
    ) %>%
    group_by(Condition) %>%
    mutate(RepLabel = paste0("N", row_number())) %>%
    ungroup()
  
  pca_df$Condition <- factor(pca_df$Condition, levels = c(lbl_cond1, lbl_cond2))
  
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, fill = Condition)) +
    stat_ellipse(geom = "polygon", alpha = 0.1, linetype = "dashed", show.legend = FALSE, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(size = 4, alpha = 0.9, shape = 21, color = "black", stroke = 0.5) +
    geom_text_repel(aes(label = RepLabel), size = 4, color = "black", show.legend = FALSE, box.padding = 0.6, point.padding = 0.5) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
    scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    theme_classic(base_size = 14) +
    labs(x = pc1_label, y = pc2_label, color = NULL, fill = NULL) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
  
  output_pdf  <- file.path(OUTPUT_DIR, paste0(base_name, "_PCA.pdf"))
  output_png  <- file.path(OUTPUT_DIR, paste0(base_name, "_PCA.png"))
  output_tiff <- file.path(OUTPUT_DIR, paste0(base_name, "_PCA.tiff"))
  
  ggsave(output_pdf,  plot = pca_plot, width = 7, height = 6, dpi = 300)
  ggsave(output_png,  plot = pca_plot, width = 7, height = 6, dpi = 300)
  ggsave(output_tiff, plot = pca_plot, width = 7, height = 6, dpi = 600, compression = "lzw")
}


#Combined Longitudinal PCA & PERMANOVA (Exact 12 Samples)


cond1 <- "Control_1_mes"
cond2 <- "Control_3_meses"
cond3 <- "STZ_1_mes"
cond4 <- "STZ_3_meses"

pattern_ctrl <- paste0("(?i)(", cond1, ".*", cond2, "|", cond2, ".*", cond1, ").*\\.csv$")
file_ctrl <- raw_files[grepl(pattern_ctrl, basename(raw_files), perl = TRUE)]

pattern_stz <- paste0("(?i)(", cond3, ".*", cond4, "|", cond4, ".*", cond3, ").*\\.csv$")
file_stz <- raw_files[grepl(pattern_stz, basename(raw_files), perl = TRUE)]


  
  extract_6_samples <- function(file_path, c1_name, c2_name) {
    df <- read_csv(file_path, skip = 2, show_col_types = FALSE) %>% clean_names()
    ids_c1 <- the_samples %>% filter(condition == c1_name) %>% pull(muestra)
    ids_c2 <- the_samples %>% filter(condition == c2_name) %>% pull(muestra)
    regex_c1 <- paste(ids_c1, collapse = "|"); pos_c1 <- which(str_detect(colnames(df), regex_c1))
    if(length(pos_c1) > length(ids_c1)) pos_c1 <- pos_c1[1:length(ids_c1)]
    regex_c2 <- paste(ids_c2, collapse = "|"); pos_c2 <- which(str_detect(colnames(df), regex_c2))
    if(length(pos_c2) > length(ids_c2)) pos_c2 <- pos_c2[1:length(ids_c2)]
    
    sub_df <- df %>% select(accession, all_of(pos_c1), all_of(pos_c2)) %>%
      filter(!is.na(accession) & accession != "") %>% distinct(accession, .keep_all = TRUE)
    colnames(sub_df) <- c("accession", paste0(c1_name, "_N", 1:length(ids_c1)), paste0(c2_name, "_N", 1:length(ids_c2)))
    return(sub_df)
  }
  
  df_ctrl <- extract_6_samples(file_ctrl[1], cond1, cond2)
  df_stz  <- extract_6_samples(file_stz[1], cond3, cond4)
  
  expr_data_all <- full_join(df_ctrl, df_stz, by = "accession") %>%
    column_to_rownames("accession") %>%
    mutate(across(everything(), ~ as.numeric(str_replace(as.character(.), ",", "."))))
  
  expr_matrix_all <- as.matrix(expr_data_all)
  expr_matrix_all[is.na(expr_matrix_all)] <- 0
  log_matrix_all <- log2(expr_matrix_all + 1)
  
  variances_all <- apply(log_matrix_all, 1, var)
  log_matrix_all <- log_matrix_all[variances_all > 0, ]
  
  pca_input_all <- t(log_matrix_all)
  pca_result_all <- prcomp(pca_input_all, scale. = TRUE)
  var_explained_all <- summary(pca_result_all)$importance[2, 1:2] * 100
  pc1_label_all <- paste0("PC1 (", round(var_explained_all[1], 1), "%)")
  pc2_label_all <- paste0("PC2 (", round(var_explained_all[2], 1), "%)")
  
  pca_df_all <- as.data.frame(pca_result_all$x) %>%
    rownames_to_column("SampleCol") %>%
    mutate(
      Disease = ifelse(str_detect(SampleCol, "(?i)Control"), "Non-diabetic", "Diabetic"),
      Time = ifelse(str_detect(SampleCol, "(?i)1_mes_N"), "1 month", "3 months"),
      RepLabel = str_extract(SampleCol, "N[1-3]$"),
      TrajectoryGroup = paste(Disease, RepLabel)
    )
  pca_df_all$Disease <- factor(pca_df_all$Disease, levels = c("Non-diabetic", "Diabetic"))
  pca_df_all$Time <- factor(pca_df_all$Time, levels = c("1 month", "3 months"))
  
  set.seed(123) 
  dist_matrix <- vegdist(pca_input_all, method = "euclidean")
  permanova_res <- adonis2(dist_matrix ~ Disease * Time, data = pca_df_all, permutations = 999)
  
  p_val_disease <- permanova_res[["Pr(>F)"]][1] 
  
  combined_pca_plot <- ggplot(pca_df_all, aes(x = PC1, y = PC2, color = Disease)) +
    stat_ellipse(aes(group = Disease), linetype = "dashed", alpha = 0.8, show.legend = FALSE, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_path(aes(group = TrajectoryGroup), arrow = arrow(length = unit(0.08, "inches"), type = "closed"), 
              color = "grey60", linewidth = 0.6, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(shape = Time), size = 4.5, alpha = 0.9, stroke = 1) +
    geom_text_repel(aes(label = RepLabel), size = 4, color = "black", show.legend = FALSE, box.padding = 0.6, point.padding = 0.5) +
    scale_color_manual(values = c("Non-diabetic" = "#56B4E9", "Diabetic" = "#E69F00")) +
    scale_shape_manual(values = c("1 month" = 16, "3 months" = 15)) + 
    theme_classic(base_size = 14) +
    labs(x = pc1_label_all, y = pc2_label_all, color = "Group", shape = "Time") +
    theme(legend.position = "right", legend.text = element_text(size = 12),
          axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"))
  
  output_pdf_comb  <- file.path(OUTPUT_DIR, "Combined_Longitudinal_PCA.pdf")
  output_png_comb  <- file.path(OUTPUT_DIR, "Combined_Longitudinal_PCA.png")
  ggsave(output_pdf_comb, plot = combined_pca_plot, width = 8, height = 6, dpi = 300)
  ggsave(output_png_comb, plot = combined_pca_plot, width = 8, height = 6, dpi = 300)
  
  raw_annotation_df <- read_csv(file_ctrl[1], skip = 2, show_col_types = FALSE) %>% clean_names() %>%
    select(accession, description) %>%
    mutate(gene_name = str_remove(str_extract(description, "GN=\\S+"), "GN="),
           gene_name = ifelse(is.na(gene_name), accession, gene_name),
           protein_desc = str_trim(str_remove(description, " OS=.*"))) %>%
    select(Accession = accession, Gene_Name = gene_name, Description = protein_desc)
  
  pca_loadings <- as.data.frame(pca_result_all$rotation) %>% rownames_to_column("Accession") %>%
    select(Accession, PC1, PC2) %>% filter(!str_detect(toupper(Accession), "REVERSE")) %>%
    left_join(raw_annotation_df, by = "Accession") %>%
    mutate(Vector_Magnitude = sqrt(PC1^2 + PC2^2)) %>%
    arrange(desc(Vector_Magnitude))
  
  loadings_excel_path <- file.path(OUTPUT_DIR, "Table_PCA_Trajectory_Signatures.xlsx")
  wb_pca <- createWorkbook(); addWorksheet(wb_pca, "PCA_Signatures")
  writeData(wb_pca, "PCA_Signatures", pca_loadings); saveWorkbook(wb_pca, loadings_excel_path, overwrite = TRUE)

