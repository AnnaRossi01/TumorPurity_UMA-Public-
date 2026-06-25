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
import subprocess
import sys


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

    def _assert_data_loaded(self):
        if self.df is None:
            raise ValueError("Data not loaded. Please call load_data() before plotting.")

    def _pretty_tool_name(self, tool_col):
        tool_name = tool_col.replace('purity_', '')
        normalized = tool_name.lower()
        if normalized == 'ascat':
            return 'ASCAT.sc'
        if normalized == 'purecn':
            return 'PureCN'
        if normalized == 'ichorcna':
            return 'IchorCNA'
        return tool_name

    def _tool_color(self, tool_col):
        return {
            'purity_ascat': '#27ae60',
            'purity_PureCN': '#e74c3c',
            'purity_ichorCNA': '#f39c12'
        }.get(tool_col, '#3498db')

    def _filtered_xy(self, tool_col):
        mask = ~np.isnan(self.df['purity_flow']) & ~np.isnan(self.df[tool_col])
        return self.df.loc[mask, 'purity_flow'], self.df.loc[mask, tool_col]

    def _tool_stats(self, tool_col):
        x, y = self._filtered_xy(tool_col)
        linreg = stats.linregress(x, y)
        return {
            'x': x,
            'y': y,
            'n': len(x),
            'pearson_r': stats.pearsonr(x, y)[0],
            'pearson_p': stats.pearsonr(x, y)[1],
            'spearman_r': stats.spearmanr(x, y)[0],
            'spearman_p': stats.spearmanr(x, y)[1],
            'slope': linreg.slope,
            'intercept': linreg.intercept,
            'r_value': linreg.rvalue,
            'p_value': linreg.pvalue,
            'std_err': linreg.stderr,
            'rmse': np.sqrt(np.mean((x - y) ** 2)),
            'mae': np.mean(np.abs(x - y)),
            'mean_diff': np.mean(x - y),
            'std_diff': np.std(x - y, axis=0),
        }

    def _barplot_for_metric(self, metric_name, title, ylabel, filename, fmt='.3f'):
        df_metrics = []
        for tool_col in self.get_tool_columns():
            stats = self._tool_stats(tool_col)
            df_metrics.append({
                'Tool': self._pretty_tool_name(tool_col),
                'Value': stats[metric_name],
                'Color': self._tool_color(tool_col)
            })

        df_metrics = pd.DataFrame(df_metrics).sort_values('Value').reset_index(drop=True)

        plt.figure(figsize=(10, 6))
        bars = sns.barplot(x='Tool', y='Value', data=df_metrics,
                           palette=df_metrics['Color'].tolist(),
                           edgecolor='black', linewidth=1.2)

        for p in bars.patches:
            bars.annotate(format(p.get_height(), fmt),
                          (p.get_x() + p.get_width() / 2., p.get_height()),
                          ha='center', va='center',
                          xytext=(0, 9), textcoords='offset points',
                          fontsize=12, fontweight='bold')

        plt.title(title, fontsize=16, pad=20)
        plt.xlabel('Computational Method', fontsize=14)
        plt.ylabel(ylabel, fontsize=14)
        plt.ylim(0, df_metrics['Value'].max() * 1.2)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        sns.despine()
        plt.tight_layout()
        plt.savefig(self.output_dir / filename, dpi=300, bbox_inches='tight')
        plt.close()

    def generate_all_plots(self):
        """Generate plots for each tool compared to flow cytometry, plus a multi-tool summary"""
        self._assert_data_loaded()
        tool_columns = self.get_tool_columns()
        for tool_col in tool_columns:
            tool_name = self._pretty_tool_name(tool_col)

            print(f"\n#### --- Generating plots for {tool_name} vs Flow Cytometry --- ####\n")
            self.scatter_plot(tool_col, tool_name)
            self.linear_regression_correlation(tool_col, tool_name)
            self.bland_altman_plot(tool_col, tool_name)
            
        print(f"\n#### --- Generating Multiple Tools Comparison Plot & Friedman Test --- ####\n")
        self.multiple_tools_comparison_plot()
        self.RMSE_plot()
        self.MAE_plot()
        self.comparative_boxplot()
        self.correlation_matrix_heatmap()
        self.metrics_bubble_heatmap()

    def scatter_plot(self, tool_col, tool_name):
        """ Explorative visualization: scatter plot with identity line to compare the two methods """
        self._assert_data_loaded()
        stats = self._tool_stats(tool_col)
        df_plot = self.df.loc[~np.isnan(self.df[tool_col]) & ~np.isnan(self.df['purity_flow']), :]

        analyses = [
            ('Default (Linregress)', f"scatter_plot_{tool_name}_default.png"),
            ('Pearson Correlation', f"scatter_plot_{tool_name}_pearson.png"),
            ('Spearman Correlation', f"scatter_plot_{tool_name}_spearman.png")
        ]

        for title_suffix, filename in analyses:
            plt.figure(figsize=(15, 10))
            sns.scatterplot(x='purity_flow', y=tool_col, data=df_plot,
                            s=100, color='blue', label='Samples', alpha=0.7)
            plt.plot([0, 1], [0, 1], 'r--')

            if 'Default' in title_suffix:
                ann = (f'Eq: $y = {stats["slope"]:.2f}x + {stats["intercept"]:.2f}$\n'
                       f'$R^2$: {stats["r_value"]**2:.2f}\np: {stats["p_value"]:.2e}')
            elif 'Pearson' in title_suffix:
                ann = f'Pearson r: {stats["pearson_r"]:.2f}\np: {stats["pearson_p"]:.2e}'
            else:
                ann = f'Spearman r: {stats["spearman_r"]:.2f}\np: {stats["spearman_p"]:.2e}'

            plt.gca().text(0.02, 0.95, ann, transform=plt.gca().transAxes,
                           verticalalignment='top', fontsize=14, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.6', facecolor='white', alpha=0.9,
                                     edgecolor='black', linewidth=1.5))
            plt.title(f"Tumor Purity Comparison ({tool_name} vs Flow Cytometry) - {title_suffix}")
            plt.xlabel('Flow Cytometry Tumor Purity')
            plt.ylabel(f'{tool_name} Tumor Purity')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.grid(True, linestyle=':', alpha=0.6)
            plt.legend(loc='lower right', prop={'size': 18, 'weight': 'bold', 'family': 'serif'},
                       frameon=True, fancybox=True, shadow=True, framealpha=0.95)
            plt.tight_layout()
            plt.savefig(self.output_dir / filename, dpi=300)
            plt.close()

    def linear_regression_correlation(self, tool_col, tool_name):
        """Linear regression analysis to quantify the relationship between the two methods"""
        self._assert_data_loaded()
        stats = self._tool_stats(tool_col)

        with open(self.output_dir / f"linear_regression_results_{tool_name}.txt", "w") as f:
            f.write(f"LINEAR REGRESSION ANALYSIS ({tool_name} vs Flow Cytometry):\n")
            f.write(f"Valid Samples (no NA): {stats['n']}\n")
            f.write(f"Equation: y = {stats['slope']:.2f}x + {stats['intercept']:.2f}\n")
            f.write(f"R-squared: {stats['r_value']**2:.2f}\n")
            f.write(f"P-value: {stats['p_value']:.3e}\n")
            f.write(f"Standard Error: {stats['std_err']:.3f}\n")
            f.write(f"Mean Absolute Error: {stats['mae']:.2%}\n")

    def bland_altman_plot(self, tool_col, tool_name):
        """Bland-Altman plot to assess the agreement between the two methods"""
        self._assert_data_loaded()
        stats = self._tool_stats(tool_col)

        mean = (stats['x'] + stats['y']) / 2
        diff = stats['x'] - stats['y']

        plt.figure(figsize=(8, 6))
        plt.scatter(mean, diff, color='blue', s=100, alpha=0.7)
        plt.axhline(stats['mean_diff'], color='red', linestyle='--')
        plt.axhline(stats['mean_diff'] + 1.96 * stats['std_diff'], color='orange', linestyle='--', label=r'+1.96SD')
        plt.axhline(stats['mean_diff'] - 1.96 * stats['std_diff'], color='orange', linestyle='--', label=r'-1.96SD')
        plt.legend(['Samples', 'Mean Difference', 'Limits of Agreement'])
        plt.xlabel(f'Mean Tumor Purity (Flow Cytometry & {tool_name})')
        plt.ylabel(f'Difference in Tumor Purity (Flow Cytometry - {tool_name})')
        plt.title(f'Bland-Altman Plot for Tumor Purity Comparison ({tool_name} vs Flow Cytometry)')
        plt.grid()
        plt.tight_layout()
        plt.savefig(self.output_dir / f"bland_altman_plot_{tool_name}.png", dpi=300)
        plt.close()

    def multiple_tools_comparison_plot(self):
        """Generates an advanced sorted plot of all tools against the Flow Cytometry baseline"""
        self._assert_data_loaded()
        cols_to_test = ['purity_flow'] + self.get_tool_columns()
        df_complete = self.df.dropna(subset=cols_to_test).copy()
        if len(df_complete) < 2:
            print("[WARNING] Not enough complete samples.")
            return

        df_complete = df_complete.sort_values(by='purity_flow').reset_index(drop=True)
        x = np.arange(len(df_complete))

        plt.figure(figsize=(22, 8))
        plt.plot(x, df_complete['purity_flow'], color='black', linewidth=2.5, zorder=3,
                 label='Flow Cytometry (Gold Standard)')

        for tool_col in self.get_tool_columns():
            plt.scatter(x, df_complete[tool_col], color=self._tool_color(tool_col), s=70,
                        alpha=0.8, label=self._pretty_tool_name(tool_col), zorder=4)

        for sample_index in x:
            plt.axvline(x=sample_index, color='gray', linestyle=':', alpha=0.3, zorder=1)

        plt.title('Tumor Purity Estimates: Computational Tools vs Flow Cytometry', fontsize=18, pad=20)
        plt.xlabel('Sample ID (Sorted by Ascending Purity)', fontsize=14, labelpad=15)
        plt.ylabel('Tumor Purity Fraction', fontsize=14)
        plt.xticks(x, df_complete['ID'], rotation=90, fontsize=8)
        plt.ylim(0, 1.05)
        plt.legend(title='Method', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=12, title_fontsize=14)
        plt.tight_layout()
        plt.savefig(self.output_dir / "multi_tool_comparison_advanced.png", dpi=300, bbox_inches='tight')
        plt.close()

    def RMSE_plot(self):
        """Generates a plot of RMSE (Root Mean Square Error) for all the tools compared to Flow Cytometry in a single bar chart"""
        self._assert_data_loaded()
        print("\n### --- Generating RMSE Comparison Bar Chart --- ####\n")
        self._barplot_for_metric('rmse', 'Root Mean Square Error (RMSE) Comparison Across Tools',
                                 'RMSE (Lower is Better)', 'RMSE_comparison_barchart.png')

    def MAE_plot(self):
        """Calcola e plotta il MAE per tutti i tools rispetto alla Flow Cytometry"""
        self._assert_data_loaded()
        print("\n#### --- Generating MAE Comparison Bar Chart --- ####\n")
        self._barplot_for_metric('mae', 'Mean Absolute Error (MAE) Comparison Across Tools',
                                 'MAE', 'MAE_comparison_barchart.png')

    def metrics_bubble_heatmap(self):
        """Generates a comprehensive dot-plot matrix comparing tools across all metrics."""
        if self.df is None:
            raise ValueError("Data not loaded.")

        print("\n#### --- Generating Comprehensive Metrics Bubble Chart --- ####\n")

        tools = self.get_tool_columns()
        metrics_dict = {'Tool': [], 'Metric': [], 'Value': [], 'Goodness': []}

        for tool_col in tools:
            mask = ~np.isnan(self.df['purity_flow']) & ~np.isnan(self.df[tool_col])
            y_true = self.df['purity_flow'][mask]
            y_pred = self.df[tool_col][mask]

            tool_name = tool_col.replace('purity_', '')
            if tool_name.lower() == 'ascat': tool_name = 'ASCAT.sc'
            elif tool_name.lower() == 'purecn': tool_name = 'PureCN'
            elif tool_name.lower() == 'ichorcna': tool_name = 'IchorCNA'

            # 1. Pearson r
            p_r, _ = stats.pearsonr(y_true, y_pred)
            metrics_dict['Tool'].append(tool_name)
            metrics_dict['Metric'].append('Pearson r')
            metrics_dict['Value'].append(p_r)
            metrics_dict['Goodness'].append(max(0, p_r)) # Normalizzato 0-1

            # 2. Spearman r
            s_r, _ = stats.spearmanr(y_true, y_pred)
            metrics_dict['Tool'].append(tool_name)
            metrics_dict['Metric'].append('Spearman r')
            metrics_dict['Value'].append(s_r)
            metrics_dict['Goodness'].append(max(0, s_r))

            # 3. R-squared
            _, _, r_value, _, _ = stats.linregress(y_true, y_pred)
            r2 = r_value**2
            metrics_dict['Tool'].append(tool_name)
            metrics_dict['Metric'].append('R-squared')
            metrics_dict['Value'].append(r2)
            metrics_dict['Goodness'].append(max(0, r2))

            # 4. MAE (L'errore viene invertito: 0 errore = 1 Goodness. Se errore >= 40% -> 0 Goodness)
            mae = np.mean(np.abs(y_true - y_pred))
            metrics_dict['Tool'].append(tool_name)
            metrics_dict['Metric'].append('MAE')
            metrics_dict['Value'].append(mae)
            metrics_dict['Goodness'].append(max(0, 1 - (mae * 2.5))) 

            # 5. RMSE
            rmse = np.sqrt(np.mean((y_true - y_pred)**2))
            metrics_dict['Tool'].append(tool_name)
            metrics_dict['Metric'].append('RMSE')
            metrics_dict['Value'].append(rmse)
            metrics_dict['Goodness'].append(max(0, 1 - (rmse * 2.5)))

        df_m = pd.DataFrame(metrics_dict)

        # --- PLOTTING ---
        fig, ax = plt.subplots(figsize=(12, 6))

        # Mappatura assi
        y_labels = df_m['Tool'].unique()
        x_labels = df_m['Metric'].unique()

        y_map = {tool: i for i, tool in enumerate(y_labels)}
        x_map = {metric: i for i, metric in enumerate(x_labels)}

        df_m['x'] = df_m['Metric'].map(x_map)
        df_m['y'] = df_m['Tool'].map(y_map)

        # Creazione dello scatter plot (Le bolle)
        # Dimensione: Goodness + 0.1 (per non far sparire del tutto i valori a zero)
        scatter = ax.scatter(x=df_m['x'], y=df_m['y'],
                             s=(df_m['Goodness'] + 0.1) * 2500, 
                             c=df_m['Goodness'], cmap='RdYlGn', 
                             alpha=0.85, edgecolors='black', linewidth=1.5)

        # Annotazioni di testo sotto ogni bolla
        for i, row in df_m.iterrows():
            if row['Metric'] in ['MAE', 'RMSE']:
                val_text = f"{row['Value']:.1%}" # In percentuale per gli errori
            else:
                val_text = f"{row['Value']:.2f}" # Frazione per correlazioni

            ax.text(row['x'], row['y'] - 0.28, val_text,
                    ha='center', va='top', fontsize=11, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none'))

        # Formattazione estetica accademica
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, fontsize=12, fontweight='bold')
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=12, fontweight='bold')

        ax.set_xlim(-0.5, len(x_labels) - 0.5)
        ax.set_ylim(-0.6, len(y_labels) - 0.2)
        
        # Griglia leggera per guidare l'occhio
        ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.4)

        plt.title('Comprehensive Tool Performance Matrix', fontsize=18, pad=25, fontweight='bold')
        sns.despine(left=True, bottom=True)
        plt.tight_layout()

        # Salvataggio
        output_file = self.output_dir / "Metrics_Bubble_Heatmap.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()


    def comparative_boxplot(self):
        """Generates a comparative boxplot of purity distributions across all methods."""
        if self.df is None:
            raise ValueError("Data not loaded.")
            
        print("\n#### --- Generating Comparative Boxplot --- ####\n")

        # Seleziono solo le colonne di interesse
        cols_to_keep = ['purity_flow'] + self.get_tool_columns()
        df_filtered = self.df[cols_to_keep].copy()
        
        # Rinomino le colonne per renderle eleganti nel grafico
        rename_dict = {'purity_flow': 'Flow Cytometry'}
        for col in self.get_tool_columns():
            tool_name = col.replace('purity_', '')
            # Normalizzo i nomi
            if tool_name.lower() == 'ascat': tool_name = 'ASCAT.sc'
            elif tool_name.lower() == 'purecn': tool_name = 'PureCN'
            elif tool_name.lower() == 'ichorcna': tool_name = 'IchorCNA'
            rename_dict[col] = tool_name
            
        df_filtered.rename(columns=rename_dict, inplace=True)

        # Trasformo il dataframe dal formato "largo" al formato "lungo" (necessario per seaborn)
        df_melted = df_filtered.melt(var_name='Method', value_name='Tumor Purity').dropna()

        # Dizionario colori (Flow Cytometry in grigio per fare da baseline)
        palette = {
            'Flow Cytometry': '#95a5a6', 
            'ASCAT.sc': '#27ae60', 
            'PureCN': '#e74c3c', 
            'IchorCNA': '#f39c12'
        }

        plt.figure(figsize=(12, 7))
        
        # 1. Disegno i Boxplot
        sns.boxplot(x='Method', y='Tumor Purity', data=df_melted, 
                    palette=palette, width=0.5, boxprops=dict(alpha=0.7), 
                    showfliers=False) # Nascondo gli outlier del boxplot perché li disegno con lo stripplot
        
        # 2. Aggiungo i singoli campioni (puntini) sopra le scatole
        sns.stripplot(x='Method', y='Tumor Purity', data=df_melted, 
                      color='black', alpha=0.4, size=4, jitter=True)

        # Estetica
        plt.title('Distribution of Tumor Purity Estimates Across Methods', fontsize=18, pad=20, fontweight='bold')
        plt.xlabel('Analytical Method', fontsize=14, fontweight='bold')
        plt.ylabel('Tumor Fraction', fontsize=14, fontweight='bold')
        plt.ylim(-0.05, 1.05)
        
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        sns.despine()

        plt.tight_layout()
        output_file = self.output_dir / "Comparative_Boxplot.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    
    def correlation_matrix_heatmap(self):
        """Generates the correlation heatmap by calling the STHDA R script."""
        if self.df is None:
            raise ValueError("Data not loaded.")

        print("\n#### --- Delegating Correlation Matrix to R (ggplot2) --- ####\n")

        # 1. Preparazione dei dati e rinomina colonne
        cols_to_keep = ['purity_flow'] + self.get_tool_columns()
        df_corr_data = self.df[cols_to_keep].copy()
        
        rename_dict = {'purity_flow': 'Flow Cytometry'}
        for col in self.get_tool_columns():
            tool_name = col.replace('purity_', '')
            if tool_name.lower() == 'ascat': tool_name = 'ASCAT.sc'
            elif tool_name.lower() == 'purecn': tool_name = 'PureCN'
            elif tool_name.lower() == 'ichorcna': tool_name = 'IchorCNA'
            rename_dict[col] = tool_name
            
        df_corr_data.rename(columns=rename_dict, inplace=True)

        # 2. CALCOLO DELLA MATRICE IN PYTHON
        corr_matrix = df_corr_data.corr(method='spearman')

        # 3. Salvataggio della matrice in un CSV temporaneo
        matrix_csv_path = self.output_dir / "spearman_matrix_for_R.csv"
        corr_matrix.to_csv(matrix_csv_path)

        r_script_path = Path(__file__).resolve().parent / "heatmap.R"
        
        try:
            result = subprocess.run(
                ["Rscript", str(r_script_path), str(matrix_csv_path), str(self.output_dir)],
                check=True,
                capture_output=True,
                text=True
            )
            print("[SUCCESS] Heatmap generated using ggplot2")
            if result.stdout:
                print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] R script failed with return code {e.returncode}")
            print(e.stderr or e.stdout)
            sys.exit(1)
        except FileNotFoundError:
            print("[ERROR] Rscript not found.")
            sys.exit(1)

    
    
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