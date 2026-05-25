plot_ascat_sc_segments <- function(seg_file, target_id, out_dir) {
  
  # 1. Importo  file ASCAT.sc
  df_seg <- fread(seg_file) %>% as.data.frame()
  
  # 2. Pulizia e Preparazione
  chr_levels <- c(as.character(1:22), "X", "Y")
  
  df_plot_prep <- df_seg %>%
    mutate(
      chr = factor(chromosome, levels = chr_levels),
      # Usiamo la colonna continua per il plot
      CN = round(total_copy_number_logr, 3), 
      ID = target_id
    ) %>%
    filter(!is.na(chr)) %>%
    arrange(chr, start)
  
  # 3. Calcolo Coordinate Cumulative (Tua logica originale)
  chrs_info <- df_plot_prep %>%
    group_by(chr) %>% 
    summarise(chr_len = max(end), .groups = "drop") %>%
    mutate(tot = cumsum(as.numeric(lag(chr_len, default = 0))))
  
  df_plot <- df_plot_prep %>%
    left_join(chrs_info, by = "chr") %>%
    mutate(start_cum = start + tot,
           end_cum = end + tot)
  
  # 4. Logica Colori (Tua soglia: 1.8 e 2.2)
  df_plot <- df_plot %>%
    mutate(status = case_when(
      CN > 2.2 ~ "Gain",
      CN < 1.8 ~ "Loss",
      TRUE ~ "Neutral"
    ))
  
  # 5. Generazione Plot
  plot_out <- ggplot(df_plot) +
    # Disegniamo i segmenti usando le coordinate cumulative
    geom_segment(aes(x = start_cum, xend = end_cum, y = CN, yend = CN, color = status), 
                 linewidth = 1.5, lend = "butt") +
    # Estetica identica al tuo script originale
    scale_color_manual(values = c("Gain" = "#E41A1C", "Loss" = "#377EB8", "Neutral" = "grey60")) +
    geom_hline(yintercept = 2, linetype = "solid", color = "black", alpha = 0.3) +
    geom_hline(yintercept = c(1.8, 2.2), linetype = "dashed", color = "black", alpha = 0.5) +
    facet_grid(. ~ chr, scales = "free_x", space = "free_x") + 
    theme_bw() +
    labs(title = paste("CNA Profile - Patient:", target_id),
         subtitle = "Continuous Signal (total_copy_number_logr)",
         x = "Chromosomes", 
         y = "Copy Number",
         color = "CNA Status") +
    scale_y_continuous(breaks = seq(0, 10, 2), limits = c(-0.2, 8.5)) +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.1, "lines"), 
      strip.background = element_rect(fill = "white")
    )
  
  # 6. Salvataggio
  file_out <- paste0("CNA_Profile_Continuous_", target_id)
  ggsave(file.path(out_dir, paste0(file_out, ".png")), plot = plot_out, width = 14, height = 7, dpi = 300)
  
  return(plot_out)
}