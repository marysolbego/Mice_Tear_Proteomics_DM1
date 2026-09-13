#PPI Network & Enrichment analysis 

library(tidyverse)
library(igraph)
library(ggraph)
library(readxl)
library(patchwork)
library(openxlsx)


file_edges <-'STRING_data/D1m_ND1m/string_interactions.tsv'

file_excel<-'results/Control_1_mes_vs_STZ_1_mes_Total_proteins.Publication_Tables.xlsx'

file_interpro<-'STRING_data/D1m_ND1m/enrichment.InterPro.tsv'

file_go <- 'STRING_data/D1m_ND1m/enrichment.Process.tsv'

#Supplementary Table All Enriched Terms

file_enrichment <- 'STRING_data/D1m_ND1m/enrichment.all.tsv'


OUTPUT_DIR <- file.path(getwd(), "results_D1m_ND1m")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)}

#Build network

nodes_df <- read_excel(file_excel, sheet = "Supp_All_DOWN") %>%
  mutate(FC_Value = as.numeric(str_replace(as.character(`Log2 Fold Change`), ",", "."))) %>%
  select(name = `Gene Name`, FC_Value) %>%
  mutate(name = str_trim(name)) %>%
# NAME CORRECTION TO MATCH STRING DATABASE 
  mutate(name = case_when(
    name == "Tf" ~ "Trf",
    name == "H1-2" ~ "H1f2",
    name == "FRRS1" ~ "Frrs1", 
    TRUE ~ name
  )) %>%
  filter(!is.na(FC_Value) & name != "" & !is.na(name)) %>% 
  distinct(name, .keep_all = TRUE) # ISOFORM  FILTER

edges_raw <- read_tsv(file_edges, show_col_types = FALSE)
edges_df <- edges_raw %>% select(from = 1, to = 2) %>% mutate(from = str_trim(from), to = str_trim(to))
edges_filtered <- edges_df %>% filter(from %in% nodes_df$name & to %in% nodes_df$name)

g <- graph_from_data_frame(d = edges_filtered, vertices = nodes_df, directed = FALSE)
max_abs <- max(abs(V(g)$FC_Value), na.rm = TRUE)

set.seed(42)
p_net <- ggraph(g, layout = "kk") + 
  geom_edge_link(color = "grey75", edge_width = 1, alpha = 0.8) +
  geom_node_point(aes(fill = FC_Value), shape = 21, size = 7, color = "grey40", stroke = 0.6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "italic", size = 3, color = "black",
                 point.padding = 0.2, bg.color = "white", bg.r = 0.15) +
  scale_fill_gradient2(
    low = "#56B4E9", mid = "#FFFFFF", high = "#E69F00", midpoint = 0,
    limits = c(-max_abs, 0), 
    name = expression(log[2]~FC)
  ) +
  theme_void() +                    
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom", 
    legend.key.width = unit(1, "cm"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )


std_cols <- c("term_id", "term_description", "observed_genes", "bg_genes", 
              "strength", "signal", "fdr", "matching_ids", "matching_labels")

df_interpro_raw <- read_tsv(file_interpro, show_col_types = FALSE)
colnames(df_interpro_raw) <- std_cols

df_go_raw <- read_tsv(file_go, show_col_types = FALSE)
colnames(df_go_raw) <- std_cols


global_min_signal <- min(c(df_interpro_raw$signal, df_go_raw$signal), na.rm = TRUE)
global_max_signal <- max(c(df_interpro_raw$signal, df_go_raw$signal), na.rm = TRUE)

min_x_limit <- floor(global_min_signal * 10) / 10 
max_x_limit <- ceiling(global_max_signal * 10) / 10 


min_genes <- min(c(df_interpro_raw$observed_genes, df_go_raw$observed_genes), na.rm = TRUE)
max_genes <- max(c(df_interpro_raw$observed_genes, df_go_raw$observed_genes), na.rm = TRUE)

min_fdr <- min(c(df_interpro_raw$fdr, df_go_raw$fdr), na.rm = TRUE)
max_fdr <- max(c(df_interpro_raw$fdr, df_go_raw$fdr), na.rm = TRUE)

#Plot enrichment

plot_enrichment <- function(df_input, title_text, top_n = 5, x_min_val, x_max_val, g_min, g_max, fdr_min, fdr_max) {
  
  df_clean <- df_input %>%
    arrange(fdr) %>%
    head(top_n) %>%
    mutate(
      term_description = str_to_sentence(term_description), 
      term_description = str_wrap(term_description, width = 40), 
      term_description = factor(term_description, levels = rev(term_description)) 
    )
  
  ggplot(df_clean, aes(x = signal, y = term_description)) +

    geom_segment(aes(x = x_min_val, xend = signal, y = term_description, yend = term_description), color = "#e0e0e0", linewidth = 1.5) +
    geom_point(aes(size = observed_genes, color = fdr), alpha = 0.9) +

    scale_color_gradientn(
      colors = c("#a1db58", "#4ea3a1", "#1d4477"),
      trans = "log10",
      name = "FDR",
      limits = c(fdr_min, fdr_max)
    ) +
    scale_size_continuous(
      name = "Gene count", 
      range = c(4, 10),
      limits = c(g_min, g_max)
    ) +

    scale_x_continuous(limits = c(x_min_val, x_max_val), breaks = seq(x_min_val, x_max_val, length.out = 5)) +
    labs(title = title_text, x = "Signal", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_line(color = "grey75", linetype = "dashed", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5, color = "black"),
      legend.background = element_rect(fill = "#f2f2f2", color = "grey60", linewidth = 0.5),
      legend.key = element_rect(fill = "#f2f2f2", color = NA),
      legend.text = element_text(color = "black", size = 9),
      legend.title = element_text(color = "black", face = "bold", size = 10)
    )
}

p_interpro <- plot_enrichment(df_interpro_raw, "InterPro Functional Domains", 5, min_x_limit, max_x_limit, min_genes, max_genes, min_fdr, max_fdr)
p_go       <- plot_enrichment(df_go_raw,       "Gene Ontology - Biological Process", 5, min_x_limit, max_x_limit, min_genes, max_genes, min_fdr, max_fdr)


# Save plots

#p_net
ggsave(file.path(OUTPUT_DIR, "Network_Plot_1Month.pdf"), plot = p_net, width = 8, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Network_Plot_1Month.png"), plot = p_net, width = 8, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Network_Plot_1Month.tiff"), plot = p_net, width = 8, height = 8, dpi = 300, compression = "lzw")

#p_interpro
ggsave(file.path(OUTPUT_DIR, "Enrichment_InterPro_1Month.pdf"), plot = p_interpro, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Enrichment_InterPro_1Month.png"), plot = p_interpro, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Enrichment_InterPro_1Month.tiff"), plot = p_interpro, width = 7, height = 5, dpi = 300, compression = "lzw")

#p_go
ggsave(file.path(OUTPUT_DIR, "Enrichment_GO_1Month.pdf"), plot = p_go, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Enrichment_GO_1Month.png"), plot = p_go, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Enrichment_GO_1Month.tiff"), plot = p_go, width = 7, height = 5, dpi = 300, compression = "lzw")



#Supplementary Table All Enriched Terms


raw_data <- read_tsv(file_enrichment, show_col_types = FALSE)


colnames(raw_data)[1] <- "Category"
colnames(raw_data) <- str_to_title(str_replace_all(colnames(raw_data), "[_\\.]", " "))

# Identify names STRING used for FDR, Description, and Strength
fdr_col <- grep("(?i)fdr|false discovery", colnames(raw_data), value = TRUE)[1]
desc_col <- grep("(?i)description", colnames(raw_data), value = TRUE)[1]
str_col <- grep("(?i)strength", colnames(raw_data), value = TRUE)[1]


processed_data <- raw_data %>%
  mutate(
    !!sym(fdr_col) := as.numeric(!!sym(fdr_col)),
    !!sym(desc_col) := str_to_sentence(!!sym(desc_col))
  )

if (!is.na(str_col)) {
  processed_data <- processed_data %>%
    mutate(!!sym(str_col) := round(as.numeric(!!sym(str_col)), 2))
}

# Order by Category and significance
processed_data <- processed_data %>% arrange(Category, !!sym(fdr_col))



wb <- createWorkbook()

header_style <- createStyle(
  fontName = "Segoe UI", fontSize = 11, fontColour = "#FFFFFF",
  fgFill = "#1F4E79", halign = "left", valign = "center", textDecoration = "bold",
  border = "TopBottomLeftRight", borderColour = "#1F4E79"
)

data_style <- createStyle(
  fontName = "Segoe UI", fontSize = 10, halign = "left", valign = "center"
)

scientific_style <- createStyle(
  fontName = "Segoe UI", fontSize = 10, numFmt = "0.00E+00", halign = "left"
)


categories <- unique(processed_data$Category)

for (cat_name in categories) {
  

  cat_df <- processed_data %>% 
    filter(Category == cat_name) %>%
    select(-Category) 
  
  sheet_name <- str_replace_all(cat_name, " ", "_")
  sheet_name <- str_remove_all(sheet_name, "\\(|\\)")
  sheet_name <- substr(sheet_name, 1, 30) 
  
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = cat_df, startRow = 1, startCol = 1)

  addStyle(wb, sheet = sheet_name, style = header_style, rows = 1, cols = 1:ncol(cat_df), gridExpand = TRUE)
  addStyle(wb, sheet = sheet_name, style = data_style, rows = 2:(nrow(cat_df)+1), cols = 1:ncol(cat_df), gridExpand = TRUE)
  

  fdr_index <- which(colnames(cat_df) == fdr_col)
  if (length(fdr_index) > 0) {
    addStyle(wb, sheet = sheet_name, style = scientific_style, rows = 2:(nrow(cat_df)+1), cols = fdr_index, gridExpand = TRUE)
  }

  setRowHeights(wb, sheet = sheet_name, rows = 1, heights = 26)
  setRowHeights(wb, sheet = sheet_name, rows = 2:(nrow(cat_df)+1), heights = 20)
  showGridLines(wb, sheet = sheet_name, show = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = 1:ncol(cat_df), widths = "auto")
}

#Save
output_excel_path <- file.path(OUTPUT_DIR, "Supplementary_Table_S2_All_Enriched_Terms.xlsx")
saveWorkbook(wb, output_excel_path, overwrite = TRUE)

