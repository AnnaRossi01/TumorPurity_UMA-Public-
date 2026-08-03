import pandas as pd
import itertools
import numpy as np

# 1. Carico i dati veri
exel = pd.read_excel('metadata.xlsx')
db_filt = exel[['ID','MAD','tot_reads','purity', 'RUN']].copy()

# FILTRO DALLA RUN 40 IN POI
db_filt = db_filt[db_filt['RUN'] >= 40].copy()

# 2. Generzione dei Range
# MAD: da 0.15 a 0.20 (incluso, quindi mi fermo a 0.21), salti di 0.01
test_mad = np.round(np.arange(0.07, 0.21, 0.01), 2).tolist()

# PURITY: da 0.50 a 1 (escluso, quindi mi fermo a 1.00), salti di 0.10
test_purity = np.round(np.arange(0.50, 1.00, 0.10), 2).tolist()

# READS: due scenari logici ("minore" o "maggiore/uguale" a 4 milioni)
test_reads = ['< 4_000_000', '>= 4_000_000']

print(f"MAD da testare: {test_mad}")
print(f"Purity da testare: {test_purity}")
print(f"Reads da testare: {test_reads}\n")

# 3. Creo tutte le combinazioni possibili
combinazioni = list(itertools.product(test_mad, test_purity, test_reads))

risultati = []

# ==========================================
# 4. CICLO DI TEST SUL DATABASE
# ==========================================
for mad_val, pur_val, read_cond in combinazioni:
    
    # Gestisco l'interruttore logico per le reads
    if read_cond == '< 4_000_000':
        mask_reads = (db_filt['tot_reads'] < 4_000_000)
    else:
        mask_reads = (db_filt['tot_reads'] >= 4_000_000)
        
    # Applico la maschera completa
    maschera = (
        (db_filt['MAD'] <= mad_val) & 
        (db_filt['purity'] >= pur_val) &
        (db_filt['purity'] < 1.0) & 
        mask_reads
    )
    
    # Conto quanti campioni sopravvivono
    campioni_salvati = maschera.sum()
    
    # Salvo il risultato
    risultati.append({
        'MAD_Max': mad_val,
        'Purity_Min': pur_val,
        'Condizione_Reads': read_cond,
        'Campioni_Rimasti': campioni_salvati
    })

# 5. Converto e stampo la tabella riassuntiva finale
df_risultati = pd.DataFrame(risultati)

# Ordino dal migliore al peggiore in termini di campioni salvati
df_risultati_ordinati = df_risultati.sort_values(by='Campioni_Rimasti', ascending=False)

print("--- CLASSIFICA COMBINAZIONI (Prime 15) ---")
print(df_risultati_ordinati.head(15).to_string(index=False))

with open('combination.txt', 'w') as combinazioni: 
    combinazioni.write(df_risultati_ordinati.to_string(index=False))

combinazioni.close()


import seaborn as sns
import matplotlib.pyplot as plt

# ==========================================
# VISUALIZZAZIONE GRID SEARCH (Heatmap)
# ==========================================

# 1. Preparo la "tela": 1 riga, 2 colonne (per le due condizioni di Reads)
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))

# --- PANNELLO 1: Reads < 4.000.000 ---
# Filtro il dataframe solo per questa condizione
df_minori = df_risultati[df_risultati['Condizione_Reads'] == '< 4_000_000']
# "Pivot" trasforma la tabella lunga in una griglia 2D (Asse Y, Asse X, Valori)
pivot_minori = df_minori.pivot(index='Purity_Min', columns='MAD_Max', values='Campioni_Rimasti')

# Creo la mappa di calore. annot=True stampa il numero esatto dentro al quadrato
sns.heatmap(pivot_minori, annot=True, cmap='Blues', fmt='g', ax=axes[0],
            cbar_kws={'label': 'Numero di Campioni'})
axes[0].set_title('Coorte con Reads < 4.000.000', fontsize=14)
# Inverto l'asse Y per avere la Purity più alta (1.00) in cima al grafico
axes[0].invert_yaxis()


# --- PANNELLO 2: Reads >= 4.000.000 ---
# Ripeto l'operazione per la seconda condizione
df_maggiori = df_risultati[df_risultati['Condizione_Reads'] == '>= 4_000_000']
pivot_maggiori = df_maggiori.pivot(index='Purity_Min', columns='MAD_Max', values='Campioni_Rimasti')

sns.heatmap(pivot_maggiori, annot=True, cmap='Blues', fmt='g', ax=axes[1],
            cbar_kws={'label': 'Numero di Campioni'})
axes[1].set_title('Coorte con Reads >= 4.000.000', fontsize=14)
axes[1].invert_yaxis()

# Ottimizzo il layout e salvo
plt.tight_layout()
plt.savefig('grid_search_heatmap_RUN40.png')