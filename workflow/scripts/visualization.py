# usr/bin/python3

# UPDATE -- 23-04-2026

# AIM -- Visualize the comparison of tumor purity from different methods (flow citometry, ASCAT.sc, ichorCNA, PureCN, etc.)
# INPUT -- A table (csv) with ID as rows and columns for Flow citometry tumor purity and tools tumor purity
# OUTPUT -- Linear regression, bland-altman plot, correlation coefficient, scatter plots, and multi-tool Friedman test

# ---- LIBRARIES ----
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import seaborn as sns
from scipy import stats     # algoritmi e strumenti matematici per l'analisi statistica
from pathlib import Path


class TumorPurityVisualization:
    """ Class to visualize the comparison of tumor purity from different methods (flow citometry and multiple tools)"""

    def __init__(self, data_file, output_dir):
        self.data_file = Path(data_file)
        self.output_dir = Path(output_dir)
        self.df = None
        self.output_dir.mkdir(exist_ok=True, parents=True)  #create the output_dir if doesn't exist

    def load_data(self):
        self.df = pd.read_csv(self.data_file, sep=',')
        self.df.columns = self.df.columns.str.strip()    

    def get_tool_columns(self):
        """Get all tool columns (exclude ID and purity_flow)"""
        return [col for col in self.df.columns if col.startswith('purity_') and col != 'purity_flow']
    
    def generate_all_plots(self):
        """Generate plots for each tool compared to flow cytometry, plus a multi-tool summary"""
        tool_columns = self.get_tool_columns()
        for tool_col in tool_columns:
            tool_name = tool_col.replace('purity_', '')
            print(f"\n#### --- Generating plots for {tool_name} vs Flow Cytometry --- ####\n")
            self.scatter_plot(tool_col, tool_name)
            self.linear_regression_correlation(tool_col, tool_name)
            self.bland_altman_plot(tool_col, tool_name)
            
        print(f"\n#### --- Generating Multiple Tools Comparison Plot & Friedman Test --- ####\n")
        self.multiple_tools_comparison_plot()

    def scatter_plot(self, tool_col, tool_name):
        """ Explorative visualization: scatter plot with identity line to compare the two methods """
        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        # rimuovere i NaN prima dei test statistici
        mask = ~np.isnan(self.df[tool_col])
        x = self.df['purity_flow'][mask]
        y = self.df[tool_col][mask]

        analyses = [
            ('Default(Linregress)', f"scatter_plot_{tool_name}_default.png"),
            ('Pearson Correlation', f"scatter_plot_{tool_name}_pearson.png"),        
            ('Spearman Correlation', f"scatter_plot_{tool_name}_spearman.png")       
        ]

        for title_suffix, filename in analyses:
            plt.figure(figsize=(15, 10))
            sns.scatterplot(x='purity_flow', y=tool_col, data=self.df, s=100, color='blue', label='Samples', alpha=0.7)
            plt.plot([0, 1], [0, 1], 'r--')  # Line of perfect agreement

            # Calcolo statistiche specifiche sui dati filtrati (x, y)
            if "Default" in title_suffix:
                slope, intercept, r, p, _ = stats.linregress(x, y)
                ann = f'Eq: $y = {slope:.2f}x + {intercept:.2f}$\n$R^2$: {r**2:.2f}\np: {p:.2e}'
            elif "Pearson" in title_suffix:
                r, p = stats.pearsonr(x, y)
                ann = f'Pearson r: {r:.2f}\np: {p:.2e}'
            else:
                r, p = stats.spearmanr(x, y)
                ann = f'Spearman r: {r:.2f}\np: {p:.2e}'

            # Annotazione e Settings
            plt.gca().text(0.02, 0.95, ann, transform=plt.gca().transAxes, 
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
            plt.title(f"Tumor Purity Comparison ({tool_name.upper()} vs Flow Cytometry) - {title_suffix}")
            plt.xlabel('Flow Cytometry Tumor Purity')
            plt.ylabel(f'{tool_name.upper()} Tumor Purity')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.grid(True, linestyle=':', alpha=0.6)
            plt.legend(loc='lower right')
            plt.tight_layout()
            
            plt.savefig(self.output_dir / filename, dpi=300)
            plt.close()

    def linear_regression_correlation(self, tool_col, tool_name):
        """Linear regression analysis to quantify the relationship between the two methods"""
        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        mask = ~np.isnan(self.df[tool_col])
        x = self.df['purity_flow'][mask]
        y = self.df[tool_col][mask]

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        
        with open(self.output_dir / f"linear_regression_results_{tool_name}.txt", "w") as f:
             f.write(f"LINEAR REGRESSION ANALYSIS ({tool_name.upper()} vs Flow Cytometry):\n")
             f.write(f"Valid Samples (no NA): {len(x)}\n")
             f.write(f"Equation: y = {slope:.2f}x + {intercept:.2f}\n")
             f.write(f"R-squared: {r_value**2:.2f}\n")
             f.write(f"P-value: {p_value:.3e}\n")
             f.write(f"Standard Error: {std_err:.3f}\n")
             f.write(f"Mean Absolute Error: {np.mean(np.abs(x - y)):.2%}\n")

    def bland_altman_plot(self, tool_col, tool_name):
        """Bland-Altman plot to assess the agreement between the two methods"""
        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        mask = ~np.isnan(self.df[tool_col])
        x = self.df['purity_flow'][mask]
        y = self.df[tool_col][mask]

        mean = (x + y) / 2
        diff = x - y
        md = np.mean(diff)
        sd = np.std(diff, axis=0)

        plt.figure(figsize=(8, 6))
        plt.scatter(mean, diff, color='blue', s=100, alpha=0.7)
        plt.axhline(md, color='red', linestyle='--')
        plt.axhline(md + 1.96*sd, color='gray', linestyle='--', c='orange', label=r'+1.96SD')
        plt.axhline(md - 1.96*sd, color='gray', linestyle='--', c='orange', label=r'-1.96SD')
        plt.legend(['Samples', 'Mean Difference', 'Limits of Agreement'])
        plt.xlabel(f'Mean Tumor Purity (Flow Cytometry & {tool_name.upper()})')
        plt.ylabel(f'Difference in Tumor Purity (Flow Cytometry - {tool_name.upper()})')
        plt.title(f'Bland-Altman Plot for Tumor Purity Comparison ({tool_name.upper()} vs Flow Cytometry)')
        plt.grid()
        plt.tight_layout()

        plt.savefig(self.output_dir / f"bland_altman_plot_{tool_name}.png", dpi=300)
        plt.close()

    def multiple_tools_comparison_plot(self):
        """Generates an advanced sorted plot of all tools against the Flow Cytometry baseline"""
        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        cols_to_test = ['purity_flow'] + self.get_tool_columns()
        df_complete = self.df.dropna(subset=cols_to_test).copy()
        
        if len(df_complete) < 2:
            print("[WARNING] Not enough complete samples.")
            return

        # 1. Ordinamento per Flow Cytometry (Cruciale per la forma a "S" o in salita)
        df_complete = df_complete.sort_values(by='purity_flow').reset_index(drop=True)
        
        # 2. Setup della figura: molto larga per far respirare i 113 ID
        plt.figure(figsize=(22, 8))
        
        # 3. Disegnamo la Flow Cytometry come LINEA CONTINUA CENTRALE
        plt.plot(df_complete['ID'], df_complete['purity_flow'], 
                 color='black', linewidth=2.5, zorder=3, 
                 label='Flow Cytometry (Gold Standard)')
        
        # 4. Disegnamo i Tools con colori DIVERSI e punti grandi per evidenziare le differenze
        # Usiamo zorder=4 per farli comparire "sopra" la linea nera
        plt.scatter(df_complete['ID'], df_complete['purity_ascat'], 
                    color='#27ae60', s=70, alpha=0.8, label='ASCAT.sc', zorder=4)
        
        plt.scatter(df_complete['ID'], df_complete['purity_PureCN'], 
                    color='#e74c3c', s=70, alpha=0.8, label='PureCN', zorder=4)
        
        plt.scatter(df_complete['ID'], df_complete['purity_ichorCNA'], 
                    color='#f39c12', s=80, alpha=0.8, label='IchorCNA', zorder=4)

        # 5. Aggiungiamo le linee verticali di separazione per ogni campione (come chiedevi)
        # Questo crea un effetto "binario" che guida l'occhio verticalmente per ogni ID
        for i in range(len(df_complete)):
            plt.axvline(x=i, color='gray', linestyle=':', alpha=0.3, zorder=1)

        # 6. Estetica e formattazione
        plt.title('Tumor Purity Estimates: Computational Tools vs Flow Cytometry', fontsize=18, pad=20)
        plt.xlabel('Sample ID (Sorted by Ascending Purity)', fontsize=14, labelpad=15)
        plt.ylabel('Tumor Purity Fraction', fontsize=14)
        
        # Rotazione degli ID a 90 gradi con font piccolo per farli entrare tutti
        plt.xticks(rotation=90, fontsize=8) 
        plt.ylim(0, 1.05)
        
        # Legenda esterna pulita
        plt.legend(title='Method', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=12, title_fontsize=14)
        
        # Layout stretto per non tagliare le etichette
        plt.tight_layout()
        plt.savefig(self.output_dir / "multi_tool_comparison_advanced.png", dpi=300, bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    import argparse

    # Configurazione dinamica per permettere l'esecuzione dello script standalone
    parser = argparse.ArgumentParser(description="Generazione grafici Purity Comparison")
    parser.add_argument("--input", required=True, help="Percorso del file CSV unito (es. purity_comparison.csv)")
    parser.add_argument("--out_dir", required=True, help="Cartella di destinazione dei grafici")
    args = parser.parse_args()

    try:
        visualization = TumorPurityVisualization(args.input, args.out_dir)
        visualization.load_data()
        visualization.generate_all_plots()
        print(f" #### --- Visualization completed, results saved in: {args.out_dir} --- ####\n")
    except Exception as e:
        print(f"[ERROR] Generazione grafici fallita: {e}")