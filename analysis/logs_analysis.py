import os
import glob
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse_purecn_log(filepath, failed_samples_report=None):
    """
    Legge un log di PureCN ed estrae metriche chimiche,
    di copertura e lo status di fallimento.
    """
    sample_id = os.path.basename(filepath).replace('_run.log', '').replace('.log', '')
    failed_samples_report = failed_samples_report or set()
    
    # Dizionario base per il campione
    data = {
        'Sample_ID': sample_id,
        'Status': 'Success', # Assumiamo successo fino a prova contraria
        'Tumor_Coverage': None,
        'Normal_Coverage': None,
        'Noise_Ratio': None,
        'AT_GC_Dropout': None,
        "High_GC_Bias": None,
        "On_Target_intervals": None,
        "Off_Target_intervals": None,
        'Total_Variants': None,
        'Rejected_LowBQ': None,
        'Segments': None,
        'Median_Segment_Size_Mb': None,
        "MAPD_Overall": None,
    }
    
    try:
        with open(filepath, 'r') as file:
            log_text = file.read()
            
            # --- STATUS DEL CAMPIONE ---
            if failed_samples_report:
                # Se esiste un report esterno, usalo come fonte autorevole.
                data['Status'] = 'Failed' if sample_id in failed_samples_report else 'Success'
            elif "Execution halted" in log_text or "PureCN failed due to QC/Noise limits" in log_text:
                data['Status'] = 'Failed'
                
            # --- ESTRAZIONE METRICHE VIA REGEX ---
            # 1. Coverage
            cov_match = re.search(r'Mean target coverages: (\d+)X \(tumor\) (\d+)X \(normal\)', log_text)
            if cov_match:
                data['Tumor_Coverage'] = float(cov_match.group(1))
                data['Normal_Coverage'] = float(cov_match.group(2))
            
            # 2. Intervals 
            intervals_match = re.search(r'Using (\d+) intervals \((\d+) on-target, (\d+) off-target\)', log_text)
            if intervals_match:
                data['On_Target_intervals'] = int(intervals_match.group(2))
                data['Off_Target_intervals'] = int(intervals_match.group(3))
            
            # 3. dropout & GC Bias 
            droupout_match = re.search(r'AT/GC dropout: ([\d\.]+) \(tumor\)', log_text)
            if droupout_match:
                data['AT_GC_Dropout'] = float(droupout_match.group(1))
            if "High GC-bias in normal and/or tunor coverage" in log_text:
                data['High_GC_Bias'] = True

            # 4. Noise Ratio (Il parametro critico per i fallimenti)
            noise_match = re.search(r'Tumor/normal noise ratio: ([\d\.]+)', log_text)
            if noise_match:
                data['Noise_Ratio'] = float(noise_match.group(1))
                
            # 5. MAPD (per capire se il campione è rumoroso o meno)
            mapd_match = re.search(r'Mean standard deviation of log-ratios: [\d\.]+ \(MAPD: ([\d\.]+)\)', log_text)
            if mapd_match:
                data['MAPD_Overall'] = float(mapd_match.group(1))

            # 6. Varianti eliminate per scarsa qualità (Il problema del BQ < 25)
            total_var = re.search(r'Found (\d+) variants in VCF file', log_text)
            if total_var:
                data['Total_Variants'] = int(total_var.group(1))
                
            low_bq = re.search(r'Removing (\d+) low quality variants with non-offset BQ < 25', log_text)
            if low_bq:
                data['Rejected_LowBQ'] = int(low_bq.group(1))
            
            var_used = re.search(r'Using (\d+) variants\.', log_text)
            if var_used:
                data['Variants_Used'] = int(var_used.group(1))

            # 7. Segmentazione (per capire se il tool ha avuto problemi a segmentare il genoma)
            segm_numb = re.search(r'Found (\d+) segments with median size of ([\d\.]+)Mb', log_text)
            if segm_numb:
                data['Segments'] = int(segm_numb.group(1))
                data['Median_Segment_Size_Mb'] = float(segm_numb.group(2))
                

    except Exception as e:
        print(f"Impossibile leggere il file {filepath}: {e}")
        
    return data


def load_failed_samples_from_report(report_path):
    """Legge il file quality_filter_report.txt e restituisce i campioni esclusi."""
    failed_samples = set()
    in_excluded_section = False
    try:
        with open(report_path, 'r') as report:
            for line in report:
                stripped = line.strip()
                if stripped.startswith('EXCLUDED SAMPLES'):
                    in_excluded_section = True
                    continue
                if stripped.startswith('INCLUDED SAMPLES'):
                    break
                if in_excluded_section and stripped.startswith('- '):
                    sample = stripped[2:]
                    if sample:
                        failed_samples.add(sample)
    except FileNotFoundError:
        print(f"Report non trovato: {report_path}. Continuo senza status esterni.")
    return failed_samples


def generate_qc_plots(df):
    # Usa una metrica disponibile per tutti i campioni, in modo che i falliti non vengano esclusi dal grafico.
    if 'MAPD_Overall' in df.columns and df['MAPD_Overall'].notna().any():
        noise_metric = 'MAPD_Overall'
    elif 'Noise_Ratio' in df.columns and df['Noise_Ratio'].notna().any():
        noise_metric = 'Noise_Ratio'
    else:
        noise_metric = 'Tumor_Coverage'

    fig = plt.figure(figsize=(15, 10))
    palette_status = {"Success": "green", "Failed": "red"}

    # 1. Noise Ratio/MAPD vs Status
    plt.subplot(2, 2, 1)

    sns.boxplot(x='Status', y=noise_metric, hue='Status', data=df, palette=palette_status, width=0.5, fliersize=0, legend=False)
    sns.stripplot(x='Status', y=noise_metric, data=df, color='black', alpha=0.6, jitter=True)
    plt.title(f'A. {noise_metric} by Status', fontweight='bold')
    plt.ylabel(noise_metric)
    plt.xlabel('Status')
    plt.legend([],[], frameon=False)

    # # Traccia empiricamente la massima tolleranza osservata nei campioni validi
    #if 'Success' in df['Status'].values:
    #    max_tol = df[df['Status']=='Success']['Noise_Ratio'].max()
    #    plt.axhline(y=max_tol, color='grey', linestyle='--', label=f'Max Tolerance (~{max_tol:.1f})')
    #    plt.legend()

    # 2. Iper segmentation (Segments e Median Segment Size)
    plt.subplot(2, 2, 2)
    sns.scatterplot(x='Segments', y='Median_Segment_Size_Mb', hue='Status', data=df, palette=palette_status, s=100)
    plt.title('B. PSCBS Hyper-segmentation Dynamics', fontweight='bold')
    plt.xlabel('Number of Copy Number Segments (Log Scale)')
    plt.ylabel('Median Segment Size (Mb)')
    plt.xscale('log') # Log scale necessaria perché i segmenti esplodono a migliaia nel fallimento


    # 3. Coverage vs quality loss
    plt.subplot(2, 2, 3)
    sns.scatterplot(x='Tumor_Coverage', y='Pct_Variants_Lost_BQ', hue='Status', data=df, palette=palette_status, edgecolor='black', s=100)
    plt.title('C. Tumor Coverage vs Quality Loss', fontweight='bold')
    plt.xlabel('Mean Tumor Coverage (X)')       
    plt.ylabel('Percentage of Variants Lost due to Low BQ (< 25)')



    #4. Intervals on/off target
    plt.subplot(2, 2, 4)
    sns.scatterplot(x='On_Target_intervals', y='Off_Target_intervals', hue='Status', data=df, palette=palette_status, edgecolor='black', s=100)
    plt.title('D. On-target vs Off-target Intervals', fontweight='bold')
    plt.xlabel('Number of On-target Intervals')
    plt.ylabel('Number of Off-target Intervals')

    

    plt.tight_layout()
    plt.savefig("./analysis/PureCN_QC_Analytics.png", dpi=300)
    print("Grafici salvati in PureCN_QC_Analytics.png")



def main():
    # Percorso relativo alla wd (TESI_PROJECT)
    log_directory = "logs/purecn" 
    
    # Cerca tutti i file .log nella cartella 'logs'
    log_files = glob.glob(os.path.join(log_directory, "*run.log"))
    
    if len(log_files) == 0:
        print(f"Nessun log trovato in '{log_directory}'")
        return
        
    print(f" Inizio il parsing di {len(log_files)} file di log trovati in locale...\n")

    report_path = os.path.join("results", "purity_comparison", "quality_filter_report.txt")
    failed_samples_report = load_failed_samples_from_report(report_path)
    if failed_samples_report:
        print(f"Trovati {len(failed_samples_report)} campioni falliti nel report: {sorted(failed_samples_report)}\n")
    else:
        print("Nessun report di failed samples trovato o report vuoto. Continuerò con la sola analisi dai log.\n")

    results = [parse_purecn_log(f, failed_samples_report=failed_samples_report) for f in log_files]
    df = pd.DataFrame(results)

    # Calcolo della percentuale di varianti rimosse
    # Usiamo un controllo di sicurezza per evitare divisioni per zero
    if 'Total_Variants' in df.columns and 'Rejected_LowBQ' in df.columns:
        df['Pct_Variants_Lost_BQ'] = (df['Rejected_LowBQ'] / df['Total_Variants']) * 100

    # Salvataggio dell'output in locale
    output_csv = "purecn_local_test_QC.csv"
    outputDir = "./analysis/"
    df.to_csv(outputDir + output_csv, index=False)

    # Stampa un riassunto elegante a terminale
    display_cols = ['Sample_ID', 'Status', 'Tumor_Coverage', 'Noise_Ratio', 'Pct_Variants_Lost_BQ', 'Segments', 'Median_Segment_Size_Mb']
    print(df.to_string(index=False, columns=[c for c in display_cols if c in df.columns]))
    print(f"\n Dati estratti e salvati con successo in: {outputDir + output_csv}")

    generate_qc_plots(df)

    print("Numbers of failed samples:", df[df['Status'] == 'Failed'].shape[0])

if __name__ == "__main__":
    main()
  