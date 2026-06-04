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
    output_dir: directory where the merged results will be saved
    """

    def __init__(self, manifest_path, flow_db_path, ascat_dir, purecn_dir, output_dir):
        " Initialize the paths of snakefile"
        self.manifest_path = Path(manifest_path)
        self.flow_db_path = Path(flow_db_path)
        self.ascat_base_dir = Path(ascat_dir)
        self.purecn_base_dir = Path(purecn_dir)
        self.output_dir = Path(output_dir)

        self.manifest_df = pd.read_csv(self.manifest_path)
        self.manifest_df['sample_id'] = self.manifest_df['sample_id'].astype(str).str.strip()
        #self.manifest_path['sample_id'] = self.manifest_df['sample_id'].astype(str).str.strip() # ensure that the sample_id column is treated as string

        self.df_flow = None
        self.df_tools = {} # dictionary for all the tools i will use (ASCAT.sc, PureCN...)
        self.df_merged = None


    def load_flow_citometry_results(self):
        """Reads the flow citometry results from an Excel file and returns a DataFrame. Select only the purity value for each patient
        # possibile problema --> l'exel di input con dati di flow dentro cartella di ASCAT_UMA, in futuro magari va modiifcato il percorso per dividere meglio i dati 
        """
        if not self.flow_db_path.exists():
            raise FileNotFoundError(f"ERROR: no excel file found at {self.flow_db_path.resolve()}")

        print(f"\n#### --- Loading Flow Cytometry --- ####\n")
        df = pd.read_excel(self.flow_db_path, sheet_name="db_purity_no7500")       # because I have more than one sheet in the excel, I need to specify which one I want to load (in this case, the one with purity results)
        self.df_flow = df[['ID', 'purity']].rename(columns={'purity': 'purity_flow'})
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
            ichor_path = str(row['ichor_res'])  # Assuming the manifest has a column 'ichor_res' with the path to the ichorCNA results
            

            if ichor_path != "NA" and Path(ichor_path).exists():
                tumor_fraction = pd.read_csv(ichor_path, sep="\t", nrows=1)['Tumor Fraction'].iloc[0]
                ichorCNA_list.append({"ID": sample_id, 'purity_ichorCNA': tumor_fraction})
            else:
                ichorCNA_list.append({"ID": sample_id, 'purity_ichorCNA': pd.NA})
        
        self.df_tools['ichorCNA'] = pd.DataFrame(ichorCNA_list)
        print(f"ichorCNA Data loaded: {len(self.df_tools['ichorCNA'])} samples\n")


    def load_PureCN_results(self):
        """Reads PureCN summary file (csv) for each sample and obtain the purity results - results into server"""
        PureCN_list = []
        print(f"#### --- Loading PureCN results --- ####\n")

        for sample_id in self.manifest_df['sample_id']:
            file_PureCN = self.purecn_base_dir / sample_id / f"{sample_id}.csv"   
                 
            if file_PureCN.exists(): 
                purity_val = pd.read_csv(file_PureCN)['Purity'].iloc[0]
                PureCN_list.append({"ID": sample_id, 'purity_PureCN': purity_val})
            else:
                PureCN_list.append({"ID": sample_id, 'purity_PureCN': pd.NA})
                print(f"  [WARNING] PureCN result not found for: {sample_id}, added NA\n")

        # Manteniamo la protezione anti-crash di Pandas
        self.df_tools['PureCN'] = pd.DataFrame(PureCN_list, columns=["ID", "purity_PureCN"])
        print(f"PureCN Data loaded: {len(self.df_tools['PureCN'])} samples\n")


    def merge_and_save(self): 
        """Merges the flow citometry and all tool results on the ID column"""
        self.df_merged = self.df_flow.copy()
        
        for tool_name, df_tool in self.df_tools.items():
            if df_tool is not None and not df_tool.empty:
                self.df_merged = pd.merge(self.df_merged, df_tool, on='ID', how='left')
        
        self.output_dir.parent.mkdir(exist_ok=True, parents=True)
        self.df_merged.to_csv(self.output_dir, sep='\t', index=False)
        print(f"#### --- Results merged and saved to: {self.output_dir} --- ####\n")
    
def main():
    parser = argparse.ArgumentParser(description="OOP Purity Comparison Parser")
    parser.add_argument("--manifest", required=True, help="Path to manifest.csv")
    parser.add_argument("--flow_db", required=True, help="Path to Excel flow cytometry DB")
    parser.add_argument("--ascat_dir", required=True, help="Base directory for ASCAT results")
    parser.add_argument("--purecn_dir", required=True, help="Base directory for PureCN results")
    parser.add_argument("--output", required=True, help="Final TSV output path")
    args = parser.parse_args()

    tumor_comparison = TumorPurityComparison(
        manifest_path=args.manifest,
        flow_db_path=args.flow_db,
        ascat_dir=args.ascat_dir,
        purecn_dir=args.purecn_dir,
        output_dir=args.output
    )
    
    tumor_comparison.load_flow_citometry_results()
    tumor_comparison.load_ASCAT_results()
    tumor_comparison.load_ichor_results()
    tumor_comparison.load_PureCN_results()
    tumor_comparison.merge_and_save()

    # Esecuzione del modulo di visualizzazione
    try:
        vis_path = Path(args.output)
        out_dir = vis_path.parent
        visualization = TumorPurityVisualization(vis_path, out_dir)

        print("-" * 100)
        visualization.load_data()
        visualization.generate_all_plots()
        print(f" #### --- Visualization completed, results saved in: {out_dir} --- ####\n")
    except Exception as e:
        print(f"[ERROR] Visualization module failed: {e}")

if __name__ == "__main__":
    main()
