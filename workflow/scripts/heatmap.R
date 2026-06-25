#!/usr/bin/env Rscript

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Il pacchetto 'ggplot2' è richiesto. Installalo con install.packages('ggplot2').")
}

library(ggplot2)

# --- INPUT ---
args <- commandArgs(trailingOnly = TRUE)
input_csv <- args[1]
out_dir <- args[2]

corr_matrix <- read.csv(input_csv, row.names = 1, check.names = FALSE)
corr_matrix <- as.matrix(corr_matrix)
corr_matrix <- round(corr_matrix, 2)

# --- FUNZIONI DI SUPPORTO ---

# 1. Riordino la matrice con il clustering gerarchico e salvo l'ordine
reorder_corr_matrix <- function(mat){
  dd <- as.dist((1 - mat) / 2)
  hc <- hclust(dd)
  ordered_mat <- mat[hc$order, hc$order]
  # Restituisco sia la matrice riordinata sia il vettore con i nomi in ordine gerarchico
  return(list(matrix = ordered_mat, order = rownames(ordered_mat)))
}

# 2. Prendo solo il triangolo inferiore (metto NA nel superiore)
get_upper_tri <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}

# 3. Applico le elaborazioni
reorder_result <- reorder_corr_matrix(corr_matrix)
corr_matrix_ordered <- reorder_result$matrix
ordered_names <- reorder_result$order

lower_tri_matrix <- get_upper_tri(corr_matrix_ordered)

# 4. Melt dei dati con correzione dei livelli per ggplot2
melt_lower_tri <- function(mat, names_order) {
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  df <- data.frame(
    Var1 = rep(row_names, times = length(col_names)),
    Var2 = rep(col_names, each = length(row_names)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  # Rimuovo le righe corrispondenti al triangolo superiore (i valori NA)
  df <- df[!is.na(df$value), ]
  
  # --- CORREZIONE CRITICA ---
  # Forziamo l'ordine gerarchico come fattori, bloccando l'ordine alfabetico di default.
  # Usiamo rev() su Var1 per far sì che l'asse Y venga disegnato dall'alto verso il basso.
  df$Var1 <- factor(df$Var1, levels = rev(names_order))
  df$Var2 <- factor(df$Var2, levels = names_order)
  
  return(df)
}

melted_matrix <- melt_lower_tri(lower_tri_matrix, ordered_names)


# --- PLOT HEATMAP ---
ggheatmap <- ggplot(melted_matrix, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = 'bold', size = 10)
  ) +
  coord_fixed() # Forza le celle a rimanere perfettamente quadrate


# Aggiungo i valori numerici e il titolo
ggheatmap <- ggheatmap +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5, fontface = "bold") +
  ggtitle("Cross-Method Spearman Correlation Matrix") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 20)))


# --- SALVATAGGIO ---
output_path <- file.path(out_dir, "Spearman_correlation_R_corrected.png")
ggsave(output_path, plot = ggheatmap, width = 8, height = 7, dpi = 300, bg = "white")

print(paste("[SUCCESS] Grafico a triangolo generato con successo:", output_path))