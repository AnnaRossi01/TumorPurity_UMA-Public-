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

    def __init__(self, manifest_path, flow_db_path, ascat_dir, purecn_dir, output_path):
        """Initialize the workflow paths and the main data structures."""
        self.manifest_path = Path(manifest_path)
        self.flow_db_path = Path(flow_db_path)
        self.ascat_base_dir = Path(ascat_dir)
        self.purecn_base_dir = Path(purecn_dir)
        self.output_path = Path(output_path)
        self.output_dir = self.output_path.parent
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.manifest_df = pd.read_csv(self.manifest_path)
        self.manifest_df['sample_id'] = self.manifest_df['sample_id'].astype(str).str.strip()
        #self.manifest_path['sample_id'] = self.manifest_df['sample_id'].astype(str).str.strip() # ensure that the sample_id column is treated as string

        # ---- DUPLICATE SAMPLE HANDLING ----
        duplicates = self.manifest_df[self.manifest_df.duplicated(subset=['sample_id'])]
        if not duplicates.empty:
            print(f"\n[WARNING] Found {len(duplicates)} duplicate samples in the manifest.csv:")
            print(duplicates['sample_id'].tolist())
            print("Removing duplicates to prevent a merge explosion...\n")
            self.manifest_df = self.manifest_df.drop_duplicates(subset=['sample_id'], keep='first')

        # Initialize quality_flag if not present (backward compatibility)
        if 'quality_flag' not in self.manifest_df.columns:
            self.manifest_df['quality_flag'] = 'PASS'

        self.df_flow = None
        self.df_tools = {}  # Dictionary for all tools (ASCAT.sc, PureCN, ichorCNA)
        self.df_merged = None
        self.failed_samples = []  # Track failed samples for reporting
        self.ichorcna_cutoff_summary = None


    def load_flow_citometry_results(self):
        """Load flow-cytometry purity values from the Excel input file."""
        if not self.flow_db_path.exists():
            raise FileNotFoundError(f"ERROR: no excel file found at {self.flow_db_path.resolve()}")

        print(f"\n#### --- Loading Flow Cytometry --- ####\n")

        excel_file = pd.ExcelFile(self.flow_db_path)
        sheet_name = excel_file.sheet_names[0] if excel_file.sheet_names else None
        df = pd.read_excel(self.flow_db_path, sheet_name=sheet_name)
        df.columns = [str(col).strip() for col in df.columns]

        id_col = next((col for col in df.columns if col.lower() == 'id'), None)
        purity_col = next((col for col in df.columns if col.lower() == 'purity'), None)

        if id_col is None or purity_col is None:
            raise ValueError(f"Expected columns 'ID' and 'purity' in {self.flow_db_path}, found: {list(df.columns)}")

        df_filtered = df[df[id_col].astype(str).str.strip().isin(self.manifest_df['sample_id'].astype(str).str.strip())].copy()
        df_filtered = df_filtered.drop_duplicates(subset=[id_col])

        df_filtered[id_col] = df_filtered[id_col].astype(str).str.strip()
        self.df_flow = df_filtered[[id_col, purity_col]].rename(columns={purity_col: 'purity_flow'})
        self.df_flow.columns = ['ID', 'purity_flow']
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
                            print(f"  [ERROR] The purity column was not found in {sample_id}.csv")
                            purity_val = pd.NA
                            
                    except Exception as e:
                        print(f"  [ERROR] Failed to read or parse {file_PureCN}: {e}")
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
                    print(f"  [WARNING] PureCN result not found for: {sample_id}; added NA\n")

            # Keep Pandas robust to missing values and parsing issues
            self.df_tools['PureCN'] = pd.DataFrame(PureCN_list, columns=["ID", "purity_PureCN"])
            print(f"PureCN Data loaded: {len(self.df_tools['PureCN'])} samples (Failed: {len(self.failed_samples)})\n")


    def merge_and_save(self):
        """Merge the flow and tool results on the sample ID and export the final table."""
        self.df_merged = self.df_flow.copy()

        for _, df_tool in self.df_tools.items():
            if df_tool is not None and not df_tool.empty:
                self.df_merged = pd.merge(self.df_merged, df_tool, on='ID', how='left')

        for col in ['purity_flow', 'purity_ascat', 'purity_ichorCNA', 'purity_PureCN']:
            if col in self.df_merged.columns:
                self.df_merged[col] = pd.to_numeric(self.df_merged[col], errors='coerce')

        self.compute_ichorcna_cutoff()

        print(f"\n#### --- Quality Summary --- ####")
        print(f"Total samples analyzed in the cohort: {len(self.df_merged)}")
        print(f"Samples exhibiting non-convergence (NA) in PureCN: {len(self.failed_samples)}\n")

        if len(self.failed_samples) > 0:
            print(f"[NOTE] The {len(self.failed_samples)} samples lacking PureCN convergence are kept in the dataset to evaluate robustness across methods.\n")

        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.df_merged.to_csv(self.output_path, sep=',', index=False)
        print(f"#### --- Results merged and saved to: {self.output_path} --- ####\n")

    def compute_ichorcna_cutoff(self):
        """Apply the ichorCNA cutoff and report MAD and correlation before/after filtering."""
        if self.df_merged is None:
            raise ValueError("No merged dataset available. Run merge_and_save() first.")

        if 'purity_ichorCNA' not in self.df_merged.columns:
            self.ichorcna_cutoff_summary = { 
                'ichorCNA_cutoff': np.nan,
                'mad_before_cutoff': np.nan,
                'mad_after_cutoff': np.nan,
                'spearman_before_cutoff': np.nan,
                'spearman_after_cutoff': np.nan,
            }
            return self.ichorcna_cutoff_summary

        ichor_values = self.df_merged['purity_ichorCNA'].dropna()
        if ichor_values.empty:
            self.ichorcna_cutoff_summary = {
                'ichorCNA_cutoff': np.nan,
                'mad_before_cutoff': np.nan,
                'mad_after_cutoff': np.nan,
                'spearman_before_cutoff': np.nan,
                'spearman_after_cutoff': np.nan,
            }
            return self.ichorcna_cutoff_summary


        cutoff = 0.2383

        self.df_merged['ichorCNA_cutoff'] = cutoff
        self.df_merged['ichorCNA_cleaned'] = self.df_merged['purity_ichorCNA'].where(
            self.df_merged['purity_ichorCNA'] > cutoff, np.nan
        )
        self.df_merged['ichorCNA_cleaned'] = pd.to_numeric(self.df_merged['ichorCNA_cleaned'], errors='coerce')

        before_pairs = self.df_merged.dropna(subset=['purity_flow', 'purity_ichorCNA'])
        after_pairs = self.df_merged.dropna(subset=['purity_flow', 'ichorCNA_cleaned'])

        n_above_cutoff = int((ichor_values > cutoff).sum())
        n_below_cutoff = int((ichor_values <= cutoff).sum())

        mad_before = float(np.median(np.abs(before_pairs['purity_ichorCNA'] - np.median(before_pairs['purity_ichorCNA']))))
        mad_after = float(np.median(np.abs(after_pairs['ichorCNA_cleaned'] - np.median(after_pairs['ichorCNA_cleaned']))))

        spearman_before = float(before_pairs['purity_flow'].corr(before_pairs['purity_ichorCNA'], method='spearman')) if len(before_pairs) > 1 else np.nan
        spearman_after = float(after_pairs['purity_flow'].corr(after_pairs['ichorCNA_cleaned'], method='spearman')) if len(after_pairs) > 1 else np.nan

        summary = {
            'ichorCNA_cutoff': cutoff,
            'n_ichorCNA_values': int(len(ichor_values)),
            'n_ichorCNA_above_cutoff': n_above_cutoff,
            'n_ichorCNA_below_or_equal_cutoff': n_below_cutoff,
            'mad_before_cutoff': mad_before,
            'mad_after_cutoff': mad_after,
            'spearman_before_cutoff': spearman_before,
            'spearman_after_cutoff': spearman_after,
            'pearson_before_cutoff': float(before_pairs['purity_flow'].corr(before_pairs['purity_ichorCNA'], method='pearson')) if len(before_pairs) > 1 else np.nan,
            'pearson_after_cutoff': float(after_pairs['purity_flow'].corr(after_pairs['ichorCNA_cleaned'], method='pearson')) if len(after_pairs) > 1 else np.nan,
        }

        self.ichorcna_cutoff_summary = summary
        summary_path = self.output_dir / 'ichorcna_cutoff_summary.csv'
        pd.DataFrame([summary]).to_csv(summary_path, index=False)

        print(f"IchorCNA cutoff: {cutoff:.4f}")
        print(f"MAD before cutoff: {mad_before:.4f}")
        print(f"MAD after cutoff:  {mad_after:.4f}")
        print(f"Spearman before cutoff: {spearman_before:.4f}")
        print(f"Spearman after cutoff:  {spearman_after:.4f}")
        print(f"Cutoff summary saved to: {summary_path}")
        return summary

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



def main():
    parser = argparse.ArgumentParser(description="OOP Purity Comparison Parser")
    parser.add_argument("--manifest", required=True, help="Path to manifest.csv")
    parser.add_argument("--flow_db", required=True, help="Path to Excel flow cytometry DB")
    parser.add_argument("--ascat_dir", required=True, help="Base directory for ASCAT results")
    parser.add_argument("--purecn_dir", required=True, help="Base directory for PureCN results")
    parser.add_argument("--output", required=True, help="Final CSV output path")
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
    )
    
    tumor_comparison.load_flow_citometry_results()
    tumor_comparison.load_ASCAT_results()
    tumor_comparison.load_ichor_results()
    tumor_comparison.load_PureCN_results()
    tumor_comparison.merge_and_save()
    tumor_comparison.save_quality_report()

    try:
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
