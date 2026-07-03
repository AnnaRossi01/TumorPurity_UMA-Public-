# !usr/bin/python3

# UPDATE-26-05-2026
#    
# AIM -- Compare the tumor purity from flow citometry, ASCAT.sc, ichorCNA and PureCN results for the same samples

# INPUT -- 1. Flow citometry results (tumor purity) for the samples --> flow_cytometry_db.xlsx 
#          2. ASCAT.sc results (tumor purity) for the same samples --> summary_project_{sample_id}.txt 
#          3. ichorCNA results (tumor purity) for the same samples --> {sample_id}.params.txt
#          4. PureCN results (tumor purity) for the same samples --> {sample_id}.csv

# OUTPUT -- A table (tsv) with ID as rows and two columns (Flow citometry tumor purity and ASCAT.sc tumor purity)


# ---- LIBRARIES ----
import argparse
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys 
from visualization import TumorPurityVisualization



class TumorPurityComparison:
    """ Class to compare tumor purity results from flow citometry and ASCAT.sc for the same samples.
    input_dir: directory where the input files are located (flow cytometry results and ASCAT.sc results)
    output_path: file path where the merged results will be saved
    """

    def __init__(self, manifest_path, flow_db_path, ascat_dir, purecn_dir, output_path, cnb_data_path=None):
        " Initialize the paths of snakefile"
        self.manifest_path = Path(manifest_path)
        self.flow_db_path = Path(flow_db_path)
        self.ascat_base_dir = Path(ascat_dir)
        self.purecn_base_dir = Path(purecn_dir)
        self.output_path = Path(output_path)
        self.cnb_data_path = Path(cnb_data_path) if cnb_data_path else None # NUOVO
        self.output_dir = self.output_path.parent
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.manifest_df = pd.read_csv(self.manifest_path)
        self.manifest_df['sample_id'] = self.manifest_df['sample_id'].astype(str).str.strip()
        #self.manifest_path['sample_id'] = self.manifest_df['sample_id'].astype(str).str.strip() # ensure that the sample_id column is treated as string

        # ---- INIZIO BLOCCO FIX & DEBUG SUI 6 CAMPIONI ----
        duplicates = self.manifest_df[self.manifest_df.duplicated(subset=['sample_id'])]
        if not duplicates.empty:
            print(f"\n[ATTENZIONE] Trovati {len(duplicates)} campioni duplicati nel manifest.csv:")
            print(duplicates['sample_id'].tolist())
            print("Rimozione dei duplicati per prevenire la 'Merge Explosion'...\n")
            self.manifest_df = self.manifest_df.drop_duplicates(subset=['sample_id'], keep='first')
        # ---- FINE BLOCCO FIX & DEBUG ----

        # Initialize quality_flag if not present (backward compatibility)
        if 'quality_flag' not in self.manifest_df.columns:
            self.manifest_df['quality_flag'] = 'PASS'

        self.df_flow = None
        self.df_tools = {} # dictionary for all the tools i will use (ASCAT.sc, PureCN...)
        self.df_merged = None
        self.failed_samples = []  # Track failed samples for reporting


    def load_flow_citometry_results(self):
        """Reads the flow citometry results from an Excel file and returns a DataFrame. Select only the purity value for each patient
        # possibile problema --> l.exel di input con dati di flow dentro cartella di ASCAT_UMA, in futuro magari va modiifcato il percorso per dividere meglio i dati 
        """
        if not self.flow_db_path.exists():
            raise FileNotFoundError(f"ERROR: no excel file found at {self.flow_db_path.resolve()}")

        print(f"\n#### --- Loading Flow Cytometry --- ####\n")
        df = pd.read_excel(self.flow_db_path)       # because I have more than one sheet in the excel, I need to specify which one I want to load (in this case, the one with purity results)
        
        # 2. Filtro usando .isin() per tenere solo le righe il cui 'ID' è nel manifest
        df_filtered = df[df['ID'].isin(self.manifest_df['sample_id'])].copy()
        
        # Rimuovo eventuali ID duplicati tenendo la prima misurazione trovata
        df_filtered = df_filtered.drop_duplicates(subset=['ID'])
        
        # 3. Seleziono solo le colonne che mi servono e rinomino
        self.df_flow = df_filtered[['ID', 'purity']].rename(columns={'purity': 'purity_flow'})
        print(f"Flow Cytometry Data loaded: {len(self.df_flow)} samples\n")

    def load_ASCAT_results(self): 
        "Reads ASCAT summary file (txt) for each sample and obtain the purity results"
        
        ascat_list = []
        print(f"#### --- Loading ASCAT.sc results --- ####\n")

        # iterate over the sample IDs in the flow citometry DataFrame and look for the corresponding ASCAT.sc results in the input directory. 
        # If found, extract the purity value and store it in a list of dictionaries, which will be converted to a DataFrame at the end.
        for sample_id in self.manifest_df['sample_id']: 
            file_ASCAT = self.ascat_base_dir / f"{sample_id}" / f"summary_project_{sample_id}.txt"

            if file_ASCAT.exists(): 
                purity_val = pd.read_csv(file_ASCAT, sep="\t")['purity'].iloc[0]
                ascat_list.append({"ID": sample_id, 'purity_ascat': purity_val})

        self.df_tools['ASCAT.sc'] = pd.DataFrame(ascat_list)
        print(f"ASCAT.sc Data loaded: {len(self.df_tools['ASCAT.sc'])} samples\n")  



    def load_ichor_results(self):
        """Reads ichorCNA summary file (txt) for each sample and obtain the purity results"""
        ichorCNA_list = []
        print(f"#### --- Loading ichorCNA results --- ####\n")

        for _,row in self.manifest_df.iterrows():
            sample_id = row['sample_id']    
            ichor_path = str(row['ichor_res'])  
            

            if ichor_path != "NA" and Path(ichor_path).exists():
                tumor_fraction = pd.read_csv(ichor_path, sep="\t", nrows=1)['Tumor Fraction'].iloc[0]
                ichorCNA_list.append({"ID": sample_id, 'purity_ichorCNA': tumor_fraction})
            else:
                ichorCNA_list.append({"ID": sample_id, 'purity_ichorCNA': pd.NA})
        
        self.df_tools['ichorCNA'] = pd.DataFrame(ichorCNA_list)
        print(f"ichorCNA Data loaded: {len(self.df_tools['ichorCNA'])} samples\n")


    def load_PureCN_results(self):
            """Reads PureCN summary file (csv) for each sample and obtain the purity results"""
            PureCN_list = []
            print(f"#### --- Loading PureCN results --- ####\n")

            for sample_id in self.manifest_df['sample_id']:
                file_PureCN = self.purecn_base_dir / sample_id / f"{sample_id}.csv"

                if file_PureCN.exists():
                    try:
                        df_purecn = pd.read_csv(file_PureCN)
                        
                        # Trova la colonna della purezza ignorando maiuscole/minuscole
                        purity_col = [col for col in df_purecn.columns if col.lower() == 'purity']
                        
                        if purity_col:
                            purity_val = df_purecn[purity_col[0]].iloc[0]
                        else:
                            print(f"  [ERROR] Colonna purity non trovata in {sample_id}.csv")
                            purity_val = pd.NA
                            
                    except Exception as e:
                        print(f"  [ERROR] Impossibile leggere o parsare {file_PureCN}: {e}")
                        purity_val = pd.NA

                    # Check if PureCN failed (NA/NaN value)
                    if pd.isna(purity_val):
                        print(f"  [FILTERED] PureCN failed for {sample_id}: insufficient somatic variants (likely PSCBS crash)")
                        self.failed_samples.append(sample_id)
                        # Mark manifest entry as FAIL
                        self.manifest_df.loc[self.manifest_df['sample_id'] == sample_id, 'quality_flag'] = 'FAIL_LOW_VARIANTS'

                    PureCN_list.append({"ID": sample_id, 'purity_PureCN': purity_val})
                else:
                    PureCN_list.append({"ID": sample_id, 'purity_PureCN': pd.NA})
                    print(f"  [WARNING] PureCN result not found for: {sample_id}, added NA\n")

            # Manteniamo la protezione anti-crash di Pandas
            self.df_tools['PureCN'] = pd.DataFrame(PureCN_list, columns=["ID", "purity_PureCN"])
            print(f"PureCN Data loaded: {len(self.df_tools['PureCN'])} samples (Failed: {len(self.failed_samples)})\n")


    def merge_and_save(self):
        """Merges the flow cytometry and all tool results on the ID column.
        All samples are kept in the final dataset to avoid survivorship bias."""
        self.df_merged = self.df_flow.copy()

        for tool_name, df_tool in self.df_tools.items():
            if df_tool is not None and not df_tool.empty:
                self.df_merged = pd.merge(self.df_merged, df_tool, on='ID', how='left')
        
        if self.cnb_data_path and self.cnb_data_path.exists():
            print(f"#### --- Integrating Genomic Status (CNB) --- ####")
            df_cnb = pd.read_csv(self.cnb_data_path)
            cnb_cols = ['ID']
            if 'Genomic_Status' in df_cnb.columns:
                cnb_cols.append('Genomic_Status')
            if 'CN_Burden_Fraction' in df_cnb.columns:
                cnb_cols.append('CN_Burden_Fraction')

            if len(cnb_cols) > 1:
                self.df_merged = pd.merge(self.df_merged, df_cnb[cnb_cols], on='ID', how='left')
                print(f"Genomic Status and CNB burden successfully integrated.\n")
            else:
                print(f"[WARNING] No usable CNB columns found in {self.cnb_data_path.name}\n")

        print(f"\n#### --- Quality Summary --- ####")
        print(f"Total samples analyzed in the cohort: {len(self.df_merged)}")
        print(f"Samples exhibiting non-convergence (NA) in PureCN: {len(self.failed_samples)}\n")

        if len(self.failed_samples) > 0:
            print(f"[NOTE] The {len(self.failed_samples)} samples lacking PureCN convergence are kept in the dataset to accurately evaluate the robustness of ASCAT.sc and IchorCNA.\n")

        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # FIX: Impostato il separatore su virgola per coerenza con l'estensione .csv
        self.df_merged.to_csv(self.output_path, sep=',', index=False)
        print(f"#### --- Results merged and saved to: {self.output_path} --- ####\n")

    def plot_burden_titration(self, output_dir=None):
        """Compute MAE across increasing CNB thresholds and save a combined plot and table."""
        if self.df_merged is None:
            raise ValueError("No merged dataset available. Run merge_and_save() first.")

        plot_dir = Path(output_dir) if output_dir else self.output_dir / "plots"
        plot_dir.mkdir(exist_ok=True, parents=True)

        if 'CN_Burden_Fraction' not in self.df_merged.columns:
            print("[WARNING] CN_Burden_Fraction column not found; skipping burden titration analysis.")
            return None

        df_clean = self.df_merged[self.df_merged['purity_flow'] < 1.0].copy()
        df_clean['CN_Burden_Fraction'] = pd.to_numeric(df_clean['CN_Burden_Fraction'], errors='coerce')
        print(f"\n#### --- Burden Titration Analysis --- ####")
        print(f"Samples retained after removing Flow=1.0: {len(df_clean)}")

        thresholds = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
        results = []
        tools = ['purity_ascat', 'purity_ichorCNA', 'purity_PureCN']

        for threshold in thresholds:
            df_subset = df_clean[df_clean['CN_Burden_Fraction'] >= threshold]
            row_data = {'Threshold': threshold, 'N_Samples': len(df_subset)}

            for tool in tools:
                if tool in df_subset.columns and len(df_subset) > 0:
                    valid_mask = df_subset[tool].notna() & df_subset['purity_flow'].notna()
                    if valid_mask.sum() > 0:
                        mae = np.mean(np.abs(df_subset.loc[valid_mask, tool] - df_subset.loc[valid_mask, 'purity_flow']))
                        row_data[tool] = mae * 100.0
                    else:
                        row_data[tool] = np.nan
                else:
                    row_data[tool] = np.nan

            results.append(row_data)

        df_results = pd.DataFrame(results)
        df_results.to_csv(plot_dir / "burden_titration_results.csv", index=False)

        fig, ax1 = plt.subplots(figsize=(10, 6))
        colors = {'purity_ascat': '#2ca02c', 'purity_ichorCNA': '#ff7f0e', 'purity_PureCN': '#d62728'}
        labels = {'purity_ascat': 'ASCAT.sc', 'purity_ichorCNA': 'IchorCNA', 'purity_PureCN': 'PureCN'}

        for tool in tools:
            ax1.plot(df_results['Threshold'] * 100, df_results[tool], marker='o', linewidth=2, color=colors[tool], label=labels[tool])

        ax1.set_xlabel('Soglia Minima di Genoma Alterato (%)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Mean Absolute Error (MAE %)', fontsize=12, fontweight='bold')
        ax1.grid(True, linestyle='--', alpha=0.5)

        ax2 = ax1.twinx()
        ax2.bar(df_results['Threshold'] * 100, df_results['N_Samples'], width=2.5, alpha=0.2, color='gray', label='N° Campioni nel Subset')
        ax2.set_ylabel('Numero di Campioni (N)', fontsize=12, color='gray', fontweight='bold')
        ax2.tick_params(axis='y', labelcolor='gray')

        for i, txt in enumerate(df_results['N_Samples']):
            ax2.annotate(txt, (df_results['Threshold'].iloc[i] * 100, df_results['N_Samples'].iloc[i]), ha='center', va='bottom', color='gray', fontsize=9)

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')

        plt.title("Dinamica dell'Errore (MAE) vs. Soglia di Instabilità Genomica", fontsize=14, fontweight='bold')
        plt.tight_layout()
        output_file = plot_dir / "MAE_Titration_vs_Burden.png"
        plt.savefig(output_file, dpi=300)
        plt.close()

        print(f"Burden titration plot saved to: {output_file}")
        print(df_results.to_string(index=False))
        return df_results

    def save_quality_report(self):
        """Save a quality filter report to accompany the final results"""
        report_path = self.output_dir / "quality_filter_report.txt"
        with open(report_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("TUMOR PURITY COMPARISON - QUALITY FILTER REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"Total samples processed: {len(self.manifest_df)}\n")
            f.write(f"Samples excluded (PureCN failure): {len(self.failed_samples)}\n")
            f.write(f"Samples in final cohort: {len(self.df_merged)}\n\n")

            if len(self.failed_samples) > 0:
                f.write("EXCLUDED SAMPLES (reason: insufficient somatic variants after PureCN filtering):\n")
                f.write("-" * 80 + "\n")
                for sample_id in self.failed_samples:
                    f.write(f"  - {sample_id}\n")
                f.write("\n")

            f.write("INCLUDED SAMPLES:\n")
            f.write("-" * 80 + "\n")
            for sample_id in self.df_merged['ID'].tolist():
                f.write(f"  - {sample_id}\n")
            f.write("\n")

            f.write("=" * 80 + "\n")
            f.write("NOTE: Excluded samples had purity values marked as 'NA' in PureCN output,\n")
            f.write("indicating insufficient high-quality variants for reliable segmentation (PSCBS crash).\n")
            f.write("These samples are also excluded from ASCAT.sc and ichorCNA analysis\n")
            f.write("to ensure statistical coherence of the comparison.\n")
            f.write("=" * 80 + "\n")

        print(f"Quality filter report saved to: {report_path}\n")

    def _run_group_analysis(self, df_group, group_name, output_dir, visualization_dir):
        """Create a stratified CSV and run visualization/metrics export for one genomic group."""
        group_file = output_dir / f"{self.output_path.stem}_{group_name.lower()}.csv"
        df_group.to_csv(group_file, index=False)

        vis_dir = visualization_dir / f"{group_name}_Genomes"
        vis_dir.mkdir(parents=True, exist_ok=True)

        print(f"Generating plots for {group_name} genomes (n={len(df_group)})...")
        vis = TumorPurityVisualization(group_file, vis_dir)
        vis.load_data()
        vis.generate_all_plots()
        metrics_df = vis.save_metrics_summary(vis_dir / f"{group_name.lower()}_metrics_summary.csv")

        return group_file, vis_dir, metrics_df

    def run_stratified_analysis(self, df_master, output_dir, visualization_dir):
        """Create separate Flat/Altered outputs and compare their metrics."""
        if 'Genomic_Status' not in df_master.columns:
            return None

        output_dir.mkdir(parents=True, exist_ok=True)
        visualization_dir.mkdir(parents=True, exist_ok=True)

        group_metrics = {}
        for group_name in ['Flat', 'Altered']:
            df_group = df_master[df_master['Genomic_Status'] == group_name].copy()
            _, _, metrics_df = self._run_group_analysis(df_group, group_name, output_dir, visualization_dir)
            metrics_df.insert(0, 'Genomic_Status', group_name)
            group_metrics[group_name] = metrics_df

        summary_df = pd.concat(group_metrics.values(), ignore_index=True)
        summary_path = output_dir / "stratified_metrics_summary.csv"
        summary_df.to_csv(summary_path, index=False)

        flat_df = group_metrics['Flat']
        altered_df = group_metrics['Altered']
        comparison_df = pd.merge(
            flat_df[['Tool', 'n_valid_samples', 'pearson_r', 'spearman_r', 'r_squared', 'mae', 'rmse']],
            altered_df[['Tool', 'n_valid_samples', 'pearson_r', 'spearman_r', 'r_squared', 'mae', 'rmse']],
            on='Tool',
            suffixes=('_Flat', '_Altered')
        )
        comparison_df['delta_pearson_r'] = comparison_df['pearson_r_Altered'] - comparison_df['pearson_r_Flat']
        comparison_df['delta_spearman_r'] = comparison_df['spearman_r_Altered'] - comparison_df['spearman_r_Flat']
        comparison_df['delta_r_squared'] = comparison_df['r_squared_Altered'] - comparison_df['r_squared_Flat']
        comparison_df['delta_mae'] = comparison_df['mae_Altered'] - comparison_df['mae_Flat']
        comparison_df['delta_rmse'] = comparison_df['rmse_Altered'] - comparison_df['rmse_Flat']

        comparison_path = output_dir / "stratified_metrics_difference.csv"
        comparison_df.to_csv(comparison_path, index=False)

        print(f"\nStratified metrics saved to: {summary_path}")
        print(f"Group comparison saved to: {comparison_path}")
        return summary_path, comparison_path


def main():
    parser = argparse.ArgumentParser(description="OOP Purity Comparison Parser")
    parser.add_argument("--manifest", required=True, help="Path to manifest.csv")
    parser.add_argument("--flow_db", required=True, help="Path to Excel flow cytometry DB")
    parser.add_argument("--ascat_dir", required=True, help="Base directory for ASCAT results")
    parser.add_argument("--purecn_dir", required=True, help="Base directory for PureCN results")
    parser.add_argument("--output", required=True, help="Final CSV output path")
    parser.add_argument("--cnb_data", default=None, help="Path to CNVkit classification results (Optional)") 
    parser.add_argument("--visualization_dir", default=None,
                        help="Directory where plotting outputs will be saved (default: <output>/plots)")
    args = parser.parse_args()

    output_file = Path(args.output)          # Es: cartella_xxx/purity.csv
    output_dir = output_file.parent          # Es: cartella_xxx/
    output_dir.mkdir(parents=True, exist_ok=True)
    visualization_dir = Path(args.visualization_dir) if args.visualization_dir else output_dir / "plots"

    tumor_comparison = TumorPurityComparison(
        manifest_path=args.manifest,
        flow_db_path=args.flow_db,
        ascat_dir=args.ascat_dir,
        purecn_dir=args.purecn_dir,
        output_path=output_file,
        cnb_data_path=args.cnb_data
    )
    
    tumor_comparison.load_flow_citometry_results()
    tumor_comparison.load_ASCAT_results()
    tumor_comparison.load_ichor_results()
    tumor_comparison.load_PureCN_results()
    tumor_comparison.merge_and_save()
    tumor_comparison.plot_burden_titration(output_dir=output_dir / "plots")
    tumor_comparison.save_quality_report()

    # Esecuzione del modulo di visualizzazione
    try:
        df_master = tumor_comparison.df_merged

        if 'Genomic_Status' in df_master.columns:
            print("-" * 100)
            print("#### --- STRATIFIED VISUALIZATION INITIALIZED --- ####")
            tumor_comparison.run_stratified_analysis(df_master, output_dir, visualization_dir)
            print(f"\n #### --- All stratified visualizations completed --- ####\n")
        else:
            vis_path = Path(args.output)
            visualization = TumorPurityVisualization(vis_path, visualization_dir)
            print("-" * 100)
            visualization.load_data()
            visualization.generate_all_plots()
            print(f" #### --- Standard visualization completed, results in: {visualization_dir} --- ####\n")

    except Exception as e:
        print(f"[ERROR] Visualization module failed: {e}")

if __name__ == "__main__":
    main()
