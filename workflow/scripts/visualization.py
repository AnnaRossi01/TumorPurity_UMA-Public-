# usr/bin/python3

# UPDATE -- 23-04-2026

# AIM -- Visualize the comparison of tumor purity from different methods (flow citometry, ASCAT.sc, ichorCNA, PureCN, etc.)

# INPUT -- A table (tsv) with ID as rows and two columns (Flow citometry tumor purity, ASCAT.sc tumor purity and other tools tumor purity)

# OUTPUT -- Linear regression, bland-altman plot, t-test, correlation coefficient, and a scatter plot comparing the two methods


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
        self.df = pd.read_csv(self.data_file, sep='\t')
    
    def get_tool_columns(self):
        """Get all tool columns (exclude ID and purity_flow)"""
        return [col for col in self.df.columns if col.startswith('purity_') and col != 'purity_flow']
    
    def generate_all_plots(self):
        """Generate plots for each tool compared to flow cytometry"""
        tool_columns = self.get_tool_columns()
        for tool_col in tool_columns:
            tool_name = tool_col.replace('purity_', '')
            print(f"\n#### --- Generating plots for {tool_name} vs Flow Cytometry --- ####\n")
            self.scatter_plot(tool_col, tool_name)
            self.linear_regression_correlation(tool_col, tool_name)
            self.bland_altman_plot(tool_col, tool_name)

    def scatter_plot(self, tool_col, tool_name):
        """ Explorative visualization: scatter plot with identity line to compare the two methods
            -both on the diagonal line --> perfect agreement between the two methods
            -points above the line --> tool estimates higher purity than flow citometry
            -points below the line --> flow citometry estimates higher purity than tool
            """

        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        # statistical analysis
        x = self.df['purity_flow']
        y = self.df[tool_col]

        # list analysis 
        analyses = [
            ('Default(Linregress)', f"scatter_plot_{tool_name}_default.png"),
            ('Pearson Correlation', f"scatter_plot_{tool_name}_pearson.png"),        
            ('Spearman Correlation', f"scatter_plot_{tool_name}_spearman.png")       
        ]

        for title_suffix, filename in analyses:
                plt.figure(figsize=(8, 6))
                sns.scatterplot(x='purity_flow', y=tool_col, data=self.df, s=100, color='blue', label='Samples', alpha=0.7)
                plt.plot([0, 1], [0, 1], 'r--')  # Line of perfect agreement

                # Calcolo statistiche specifiche
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
        """Linear regression analysis to quantify the relationship between the two methods, not visual but with statistical output (slope, intercept, R-squared, p-value, standard error)
        -slope close to 1 and intercept close to 0 --> strong agreement between the two methods
        -R-squared close to 1 --> strong correlation between the two methods, say how much of the variance in tool purity can be explained by flow citometry purity  
        -p-value < 0.05 --> statistically significant relationship between the two methods
        -standard error --> precision of the slope estimate (lower is better)
        """

        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        x = self.df['purity_flow']
        y = self.df[tool_col]

        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        
        with open(self.output_dir / f"linear_regression_results_{tool_name}.txt", "w") as f:
             f.write(f"LINEAR REGRESSION ANALYSIS ({tool_name.upper()} vs Flow Cytometry):\n")
             f.write(f"Equation: y = {slope:.2f}x + {intercept:.2f}\n")
             f.write(f"R-squared: {r_value**2:.2f}\n")
             f.write(f"P-value: {p_value:.3e}\n")
             f.write(f"Standard Error: {std_err:.3f}\n")
             f.write(f"Mean Absolute Error: {np.mean(np.abs(x - y)):.2%}\n")
        f.close()


    def bland_altman_plot(self, tool_col, tool_name):
        """Bland-Altman plot to assess the agreement between the two methods, visualizing the difference between the two methods against their average
        - red line: if it's near zero --> good agreement between the two methods
        - limits of agreement (mean difference ± 1.96*SD): if most points are within these limits --> good agreement between the two methods
        """
        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")
        
        x = self.df['purity_flow']
        y = self.df[tool_col]

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

        # saving figure
        plt.savefig(self.output_dir / f"bland_altman_plot_{tool_name}.png", dpi=300)
        plt.close()

        

if __name__ == "__main__":
    data_path = "tumor_purity_comparison/tumor_purity_comparison.tsv"
    out_dir = "tumor_purity_comparison"

    visualization = TumorPurityVisualization(data_path, out_dir)
    visualization.load_data()
    visualization.generate_all_plots()
    print(f" #### --- Visualization completed, results saved in: {out_dir} --- ####\n")
    