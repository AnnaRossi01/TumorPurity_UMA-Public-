# !usr/bin/python3

# UPDATE-18-05-2026
#    
# AIM -- Compare the tumor purity from flow citometry, ASCAT.sc, ichorCNA and PureCN results for the same samples

# INPUT -- 1. Flow citometry results (tumor purity) for the samples --> flow_cytometry_db.xlsx 
#          2. ASCAT.sc results (tumor purity) for the same samples --> summary_project_{sample_id}.txt 
#          3. ichorCNA results (tumor purity) for the same samples --> {sample_id}.params.txt
#          4. PureCN results (tumor purity) for the same samples --> {sample_id}.csv

# OUTPUT -- A table (tsv) with ID as rows and two columns (Flow citometry tumor purity and ASCAT.sc tumor purity)


# ---- LIBRARIES ----
import pandas as pd
from pathlib import Path        # for handling file paths
import numpy as np
import matplotlib.pyplot as plt
import sys 
from visualization import TumorPurityVisualization


class TumorPurityComparison:
    """ Class to compare tumor purity results from flow citometry and ASCAT.sc for the same samples.
    input_dir: directory where the input files are located (flow cytometry results and ASCAT.sc results)
    output_dir: directory where the merged results will be saved
    """

    def __init__(self, input_dir, output_dir):
        " Initialize the paths and data containers"
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.df_flow = None
        self.df_tools = {} # dictionary for all the tools i will use (ASCAT.sc, PureCN...)
        self.df_merged = None


    def load_flow_citometry_results(self,file_flow, sheet="db_purity_0.90"):
        """Reads the flow citometry results from an Excel file and returns a DataFrame. Select only the purity value for each patient
        # possibile problema --> l'exel di input con dati di flow dentro cartella di ASCAT_UMA, in futuro magari va modiifcato il percorso per dividere meglio i dati 
        """

        # Construct the full path to the file
        path = self.input_dir / "db_uma" / file_flow
        
        if not path.exists():
            raise FileNotFoundError(f"Flow data not found at: {path}")

        print(f"\n#### --- Loading Flow Cytometry --- ####\n")
        df = pd.read_excel(path, sheet_name="db_purity_0.90")       # because I have more than one sheet in the excel, I need to specify which one I want to load (in this case, the one with purity results)
        self.df_flow = df[['ID', 'purity']].rename(columns={'purity': 'purity_flow'})
        print(f"Flow Cytometry Data loaded: {len(self.df_flow)} samples\n")


    def load_ASCAT_results(self): 
        "Reads ASCAT summary file (txt) for each sample and obtain the purity results"
        
        ascat_list = []
        print(f"#### --- Loading ASCAT.sc results --- ####\n")

        # iterate over the sample IDs in the flow citometry DataFrame and look for the corresponding ASCAT.sc results in the input directory. 
        # If found, extract the purity value and store it in a list of dictionaries, which will be converted to a DataFrame at the end.
        for sample_id in self.df_flow['ID']: 
            file_ASCAT = self.input_dir / f"output_{sample_id}" / f"summary_project_{sample_id}.txt"

            if file_ASCAT.exists(): 
                purity_val = pd.read_csv(file_ASCAT, sep="\t")['purity'].iloc[0]
                ascat_list.append({"ID": sample_id, 'purity_ascat': purity_val})

        self.df_tools['ASCAT.sc'] = pd.DataFrame(ascat_list)
        print(f"ASCAT.sc Data loaded: {len(self.df_tools['ASCAT.sc'])} samples\n")  



    def load_ichor_results(self):
        """Reads ichorCNA summary file (txt) for each sample and obtain the purity results"""
        ichorCNA_list = []
        print(f"#### --- Loading ichorCNA results --- ####\n")

        for sample_id in self.df_flow['ID']:
            base_dir = Path("/mnt/c/Users/rossi/comparison_tools")
            file_ichorCNA = base_dir / "ichorCNA" / f"{sample_id}.params.txt"

            if file_ichorCNA.exists(): 
                tumor_fraction = pd.read_csv(file_ichorCNA, sep="\t", nrows=1)['Tumor Fraction'].iloc[0]
                ichorCNA_list.append({"ID": sample_id, 'purity_ichorCNA': tumor_fraction})

        self.df_tools['ichorCNA'] = pd.DataFrame(ichorCNA_list)
        print(f"ichorCNA Data loaded: {len(self.df_tools['ichorCNA'])} samples\n")

    def load_PureCN_results(self):
        """Reads PureCN summary file (csv) for each sample and obtain the purity results - results into server"""
        PureCN_list = []
        print(f"#### --- Loading PureCN results --- ####\n")


        # ASSICURATI DI METTERE IL TUO UTENTE REALE QUI
        mounted_drive_base = Path("/home/anna_rossi/Server_drive02")
        server_results_base = mounted_drive_base / "PureCN" / "output" / "purecn_results"

        for sample_id in self.df_flow['ID']:
            sample_id_str = str(sample_id).strip() # .strip() rimuove eventuali spazi invisibili
            
            file_PureCN = server_results_base / sample_id_str / f"{sample_id_str}.csv"

            # ---------------------------------------------------------
            # DEBUGGING
            # ---------------------------------------------------------
            #print(f"[DEBUG PureCN] ID Campione analizzato: '{sample_id_str}'")
            #print(f"[DEBUG PureCN] Path assoluto generato: {file_PureCN}")
            #print(f"[DEBUG PureCN] Il file esiste fisicamente? {file_PureCN.exists()}\n")
            # ---------------------------------------------------------

            if file_PureCN.exists(): 
                purity_val = pd.read_csv(file_PureCN)['Purity'].iloc[0]
                PureCN_list.append({"ID": sample_id_str, 'purity_PureCN': purity_val})
            else:
                print(f"  [WARNING] PureCN result not found for: {sample_id_str}\n")

        # Manteniamo la protezione anti-crash di Pandas
        self.df_tools['PureCN'] = pd.DataFrame(PureCN_list, columns=["ID", "purity_PureCN"])
        print(f"PureCN Data loaded: {len(self.df_tools['PureCN'])} samples\n")



    def merge_and_save(self, filename="comparison_purity.tsv"): 
        """Merges the flow citometry and all tool results on the ID column"""
        
        self.df_merged = self.df_flow.copy()
        
        # Merge all tools dynamically - the pd.merge function can merge only two df at the time, so I have to use a iterative method to merge all the tools in the dictionary
        for tool_name, df_tool in self.df_tools.items():
            self.df_merged = pd.merge(self.df_merged, df_tool, on="ID", how="left")

        self.output_dir.mkdir(exist_ok=True, parents=True)  #create the output_dir if doesn't exist
        out_path = self.output_dir / filename

        self.df_merged.to_csv(out_path, sep='\t', index=False)
        print(f"#### --- Results merged and saved to: {out_path} --- ####\n")



if __name__ == "__main__":
    INPUT_DIR = Path("/home/anna_rossi/ASCAT_UMA")
    OUTPUT_DIR = Path("/home/anna_rossi/comparison_tools/output")
    output_filename = "tumor_purity_comparison.tsv"

    tumor_comparison = TumorPurityComparison(INPUT_DIR, OUTPUT_DIR)
    tumor_comparison.load_flow_citometry_results("flow_cytometry_db.xlsx")
    tumor_comparison.load_ASCAT_results()
    tumor_comparison.load_ichor_results()
    tumor_comparison.load_PureCN_results()
    tumor_comparison.merge_and_save(filename=output_filename)
    print(tumor_comparison.df_merged.to_string(index=False))

    # Visualization
    vis_path = OUTPUT_DIR / output_filename
    visualization = TumorPurityVisualization(vis_path, OUTPUT_DIR)

    print("-"*100)
    visualization.load_data()
    visualization.generate_all_plots()
    print(f" #### --- Visualization completed, results saved in: {OUTPUT_DIR} --- ####\n")


