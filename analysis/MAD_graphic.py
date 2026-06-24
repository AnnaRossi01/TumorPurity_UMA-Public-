import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# --- GESTIONE DINAMICA DELLE DIRECTORY ---
# 1. Trova la cartella esatta in cui si trova questo script (analysis/)
script_dir = Path(__file__).resolve().parent
# 2. Definisce la root directory del progetto (un livello sopra analysis/)
root_dir = script_dir.parent

# 3. Costruisce i percorsi in modo inattaccabile
csv_path = script_dir / "purecn_local_test_QC.csv"
excel_path = root_dir / "flow_cytometry_db.xlsx"

# 1. Caricamento dati
data_files = pd.read_csv(csv_path)
exel = pd.read_excel(excel_path, sheet_name="db_purity_no1")

# 2. Definizione funzione per il grafico
def plot_mad_distribution(df, title_suffix):
    # Unione con i dati di riferimento
    merged = pd.merge(df, exel, left_on='Sample_ID', right_on='ID', how='inner')

    plt.figure(figsize=(8, 5))

    # Grafico con seaborn (Migliorato l'uso della palette mappando i valori esatti)
    sns.histplot(data=merged, x='MAD', hue='Status', kde=True,
                 palette={'Failed': '#e63a01', 'Success': '#2ca843'}, 
                 element="step", alpha=0.3)

    # Linee delle mediane
    for status, color in zip(['Failed', 'Success'], ["#e63a01", "#2ca843"]):
        # Aggiunto un controllo di sicurezza nel caso un gruppo sia vuoto
        subset = merged[merged['Status'] == status]
        if not subset.empty:
            median = subset['MAD'].median()
            plt.axvline(median, color=color, linestyle='--', label=f'Median {status}: {median:.2f}')

    plt.title(f'MAD Distribution: {title_suffix}')
    plt.xlabel('MAD (Mean Absolute Difference)')
    plt.ylabel('Count')
    plt.legend()
    plt.tight_layout()
    
    # Salva il grafico dinamicamente nella stessa cartella dello script (analysis/)
    output_img = script_dir / f'MAD_Distribution_{title_suffix}.png'
    plt.savefig(output_img, dpi=300)
    plt.close() # Buona pratica: chiudere la figura per liberare memoria

# 3. Esecuzione per ogni file
if __name__ == "__main__":
    plot_mad_distribution(data_files, "Overall")