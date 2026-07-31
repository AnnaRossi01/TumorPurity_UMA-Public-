import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. LETTURA DATI
exel = pd.read_excel('metadata.xlsx')
db_filt = exel[['ID','MAD','tot_reads','purity']].copy()

# ==========================================
# 2. ELABORAZIONE DATI (Creazione dei gruppi)
# ==========================================

# Filtro MAD
binsMAD = [-np.inf, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, np.inf]
labelsMAD = ['< 0.15', '0.15-0.16', '0.16-0.17', '0.17-0.18', '0.18-0.19', '0.19-0.20', '> 0.20']
db_filt['MAD_Group'] = pd.cut(db_filt['MAD'], bins=binsMAD, labels=labelsMAD)
freq_MAD = db_filt['MAD_Group'].value_counts().reindex(labelsMAD).fillna(0)

# Filtro PURITY
binsPUR = [-np.inf, 0.50, 0.60, 0.70, 0.80, 0.90, 0.99, np.inf]
labelsPUR = ["<0.50", "0.50-0.60", "0.60-0.70", "0.70-0.80", "0.80-0.90", "0.90-0.99", "1"]
db_filt['Purity_Group'] = pd.cut(db_filt['purity'], bins=binsPUR, labels=labelsPUR)
freq_PUR = db_filt['Purity_Group'].value_counts().reindex(labelsPUR).fillna(0)

# Filtro READS (Completato!)
binsREADS = [-np.inf, 4_000_000, np.inf]
labelsREADS = ['<= 4.000.000', '> 4.000.000']
db_filt['Reads_Group'] = pd.cut(db_filt['tot_reads'], bins=binsREADS, labels=labelsREADS)
freq_READS = db_filt['Reads_Group'].value_counts().reindex(labelsREADS).fillna(0)


# ==========================================
# 3. VISUALIZZAZIONE (Un'unica immagine)
# ==========================================

# Creo una "tela" (fig) larga 18 pollici e alta 6, divisa in 1 riga e 3 colonne (axes)
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))

# -- GRAFICO 1: MAD (assegnato ad axes[0]) --
freq_MAD.plot(kind='bar', color='#4C72B0', edgecolor='black', ax=axes[0])
axes[0].bar_label(axes[0].containers[0], padding=3)
axes[0].set_title('Frequency of samples for MAD value')
axes[0].set_xlabel('MAD intervals')
axes[0].set_ylabel('Number of samples')
axes[0].tick_params(axis='x', rotation=45)

# -- GRAFICO 2: PURITY (assegnato ad axes[1]) --
freq_PUR.plot(kind='bar', color='#4C72B0', edgecolor='black', ax=axes[1])
axes[1].bar_label(axes[1].containers[0], padding=3)
axes[1].set_title('Frequency of samples for Purity value')
axes[1].set_xlabel('Purity intervals')
axes[1].tick_params(axis='x', rotation=45)

# -- GRAFICO 3: TOTAL READS (assegnato ad axes[2]) --
freq_READS.plot(kind='bar', color='#4C72B0', edgecolor='black', ax=axes[2])
axes[2].bar_label(axes[2].containers[0], padding=3)
axes[2].set_title('Frequency of samples for Total Reads')
axes[2].set_xlabel('Total Reads intervals')
# Qui metto rotazione 0 perché le etichette sono solo due e leggibilissime in orizzontale
axes[2].tick_params(axis='x', rotation=0) 

# Ottimizzo gli spazi in modo che le scritte non si sovrappongano
plt.tight_layout()

# Salvo l'immagine unica
plt.savefig('report_qualita_campioni.png')


# TEST ULTIVARIATO

# ==========================================
# 1. DEFINISCO I PARAMETRI DEL TEST CORRENTE
# ==========================================
TEST_MAD_MAX = 0.15          # Voglio un MAD basso (es. minore o uguale a 0.15)
TEST_PURITY_MIN = 0.80       # Voglio un'alta purezza (es. maggiore o uguale a 0.80)
TEST_READS_MIN = 4_000_000   # Voglio abbastanza reads (es. maggiore o uguale a 4M)

# ==========================================
# 2. APPLICO LA MASCHERA MULTIPLA
# ==========================================

maschera_coorte = (
    (db_filt['MAD'] <= TEST_MAD_MAX) & 
    (db_filt['purity'] >= TEST_PURITY_MIN) & 
    (db_filt['purity'] != 1) &
    (db_filt['tot_reads'] >= TEST_READS_MIN)
)

# Creo il DataFrame finale contenente SOLO i campioni perfetti
coorte_studio = db_filt[maschera_coorte].copy()

# ==========================================
# 4. VALUTO IL RISULTATO
# ==========================================
campioni_iniziali = len(db_filt)
campioni_rimasti = len(coorte_studio)
tasso_ritenzione = (campioni_rimasti / campioni_iniziali) * 100

print(f"--- RISULTATI TEST COORTE ---")
print(f"Campioni totali iniziali: {campioni_iniziali}")
print(f"Campioni sopravvissuti ai filtri: {campioni_rimasti}")
print(f"Tasso di ritenzione: {tasso_ritenzione:.2f}%\n")

# Mostro i primi 5 ID scelti
print("Primi campioni inclusi:")
print(coorte_studio.head())