# !/usr/bin/r

# input_file --> load_input_file
#last update --> 02/02+3/2026

# CARICA LIBRERIE
library(data.table) # fread per legere il dataser
library(tidyverse) # manipolazione dato (marge, filro, visual)
library(openxlsx) # importare e leggere exel

# ----MAIN----
# set enviroment
wd <- setwd('C:/Users/rossi/ASCAT_UMA/')

# import dataset and check 
df <- read.xlsx('db_uma/metadata_filt.xlsx')
str(df) # info df, numero righe, num colonne, tipo di dato (dtype python)
#colnames(df) per nome colonne 

df_filt <- df %>%
  filter(purity >= 0.9 & purity != 1, 
         tot_reads >= 4000000, 
         MAD <= 0.15)


