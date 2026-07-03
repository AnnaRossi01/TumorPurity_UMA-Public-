import pandas as pd
import subprocess
import os
from pathlib import Path

# ==========================================
# 1. CONFIGURAZIONI INIZIALI
# ==========================================
excel_metadata = '../metadata_GoldStandard_filtrato.xlsx'      # Il tuo file Excel
output_csv = 'CNVkit_classification_results.csv' # Il file finale con le classificazioni
cartella_seg_output = 'risultati_seg' # Cartella dove salveremo i .seg convertiti

log2_threshold = 0.2  # Soglia di seg.mean per considerare il segmento alterato (amplificato o deleto)
cnb_threshold = 0.15  # Soglia di genoma alterato per classificare "Altered" (15%)

# Creiamo la cartella di output per i .seg se non esiste
Path(cartella_seg_output).mkdir(parents=True, exist_ok=True)

# ==========================================
# 2. LETTURA DEL DATABASE
# ==========================================
print("Loading metadata...")
df_meta = pd.read_excel(excel_metadata)

risultati_finali = []
campioni_processati = 0
campioni_mancanti = 0

# ==========================================
# 3. MOTORE DI ELABORAZIONE
# ==========================================
for index, row in df_meta.iterrows():
    sample_id = str(row['ID']).strip()
    
    # Formattazione della RUN (gestisce sia se scritta come numero intero 40, sia come stringa)
    run_num = str(row['RUN']).strip()
    
    # Costruzione dinamica del path come da tua indicazione
    cns_path = Path(f"/mnt/drive02/UMA/run_{run_num}/CNVkit/remDup_recal_{sample_id}.cns")
    seg_path = Path(f"{cartella_seg_output}/{sample_id}.seg")
    
    # Controllo esistenza file
    if not cns_path.exists():
        print(f"[WARNING] File not found, skipping: {cns_path}")
        campioni_mancanti += 1
        continue

    # --- STEP 1: CONVERSIONE DA CNS A SEG TRAMITE CNVKIT ---
    try:
        # Costruiamo il comando per il terminale: cnvkit.py export seg input.cns -o output.seg
        cmd = ["cnvkit.py", "export", "seg", str(cns_path), "-o", str(seg_path)]
        
        # Eseguiamo il comando in modo invisibile
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed export of CNVkit for {sample_id}. Error: {e.stderr}")
        continue

    # --- STEP 2: ANALISI E CLASSIFICAZIONE DEL .SEG ---
    try:
        # Carico il .seg generato (separato da tab)
        df_seg = pd.read_csv(seg_path, sep='\t')
        
        # Estraggo le metriche
        # 'num.mark' indica il numero di sonde/bin che compongono quel segmento
        total_markers = df_seg['num.mark'].sum()
        
        # Filtro per trovare solo i segmenti alterati (sopra 0.2 o sotto -0.2)
        segmenti_alterati = df_seg[df_seg['seg.mean'].abs() > log2_threshold]
        altered_markers = segmenti_alterati['num.mark'].sum()
        
        # Calcolo matematico del Burden
        if total_markers > 0:
            cn_burden = altered_markers / total_markers
        else:
            cn_burden = 0.0
            
        # Classificazione 
        status = 'Altered' if cn_burden >= cnb_threshold else 'Flat'
        
        # Salvo i dati
        risultati_finali.append({
            'ID': sample_id,
            'RUN': run_num,
            'Total_Markers': total_markers,
            'Altered_Markers': altered_markers,
            'CN_Burden_Fraction': round(cn_burden, 4),
            'Genomic_Status': status
        })
        
        campioni_processati += 1
        print(f"[OK] {sample_id} processato: {status} (Burden: {cn_burden:.1%})")

    except Exception as e:
        print(f"[ERROR] Problem reading .seg file for {sample_id}: {e}")

# ==========================================
# 4. SALVATAGGIO DEI RISULTATI
# ==========================================
if risultati_finali:
    df_risultati = pd.DataFrame(risultati_finali)
    df_risultati.to_csv(output_csv, index=False)
    print("\n" + "="*40)
    print("ANALYSIS COMPLETED SUCCESSFULLY")
    print(f"Samples processed and classified: {campioni_processati}")
    print(f"Samples not found (missing file): {campioni_mancanti}")
    print(f"Summary file saved to: {output_csv}")
    print("="*40)
    
    # Riassunto statistico rapido della coorte
    conteggi = df_risultati['Genomic_Status'].value_counts()
    print("\nResoconto Coorte:")
    print(conteggi.to_string())
else:
    print("\n[CRITICAL] No data processed. Check the paths and the Excel file.")