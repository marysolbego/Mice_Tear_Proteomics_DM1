#Protein-Protein Interaction Network Diabetic signature

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(readxl)


file_edges <-'STRING_data/Diabetic_signature/string_interactions.tsv'

file_nodes <- 'results/Venn_Filtered_Proteins_Lists.xlsx'

OUTPUT_DIR <- file.path(getwd(), "results_diabetic_signature")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)}

# Load Data 
edges_df <- read_tsv(file_edges, show_col_types = FALSE)

nodes_up <- read_excel(file_nodes, sheet = "UP_Diabetic_Exclusive")
nodes_down <- read_excel(file_nodes, sheet = "DOWN_Diabetic_Exclusive")
nodes_df <- bind_rows(nodes_up, nodes_down)

# Clean numeric values
nodes_df <- nodes_df %>%
  mutate(FC_Diabetic = as.numeric(str_replace(as.character(`Log2 Fold Change`), ",", "."))) %>%
  filter(!is.na(FC_Diabetic))

# Format Nodes and Filter Edges 
colnames(edges_df)[1:2] <- c("from", "to")
string_names_present <- unique(c(edges_df$from, edges_df$to))


nodes_df_prepared <- nodes_df %>%
  mutate(name = case_when(
    `Gene Name` %in% string_names_present ~ `Gene Name`,
    Accession %in% string_names_present ~ Accession,
    TRUE ~ `Gene Name` 
  )) %>%
  relocate(name)

valid_node_names <- nodes_df_prepared$name


edges_filtered <- edges_df %>%
  filter(from %in% valid_node_names & to %in% valid_node_names)

#Build Network Structure 

g <- graph_from_data_frame(d = edges_filtered, vertices = nodes_df_prepared, directed = FALSE)

max_abs <- max(abs(V(g)$FC_Diabetic), na.rm = TRUE)

#Network visualization
set.seed(42)
p <- ggraph(g, layout = "kk") + 
  geom_edge_link(color = "grey75", edge_width = 1, alpha = 0.8) +
  geom_node_point(aes(fill = FC_Diabetic), shape = 21, size = 8, color = "grey40", stroke = 0.6) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "italic", size = 3.5, 
                 point.padding = 0.2, bg.color = "white", bg.r = 0.15) +
  scale_fill_gradient2(
    low = "#56B4E9", mid = "#FFFFFF", high = "#E69F00", midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = expression(log[2]~Fold~Change),
    guide = guide_colorbar(
      title.position = "top", 
      title.hjust = 0.5,
      barwidth = unit(7, "cm"),       
      barheight = unit(0.5, "cm"), 
      frame.colour = "black",         
      ticks.colour = "black"
    )
  ) +
  theme_void() +                      
  theme(
    legend.position = "bottom",       
    legend.margin = margin(t = 20),   
    plot.margin = margin(15, 15, 15, 15)
  )

# Save 
ruta_net_pdf  <- file.path(OUTPUT_DIR, "Network_Diabetic_Signature.pdf")
ruta_net_png  <- file.path(OUTPUT_DIR, "Network_Diabetic_Signature.png")
ruta_net_tiff <- file.path(OUTPUT_DIR, "Network_Diabetic_Signature.tiff")

ggsave(ruta_net_pdf,  plot = p, width = 8, height = 9, dpi = 300)
ggsave(ruta_net_png,  plot = p, width = 8, height = 9, dpi = 300)
ggsave(ruta_net_tiff, plot = p, width = 8, height = 9, dpi = 300, compression = "lzw")
