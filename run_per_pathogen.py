import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

# Define paths
current_dir = Path(__file__).parent
results_csv = current_dir / "validation" / "output" / "validation_results_scored.csv"
output_per_pathogen_csv = current_dir / "validation" / "output" / "per_pathogen_summary.csv"

df = pd.read_csv(results_csv)
pathogens = sorted(df['Pathogen'].unique())
pathogen_summary = []

for pathogen in pathogens:
    sub = df[df['Pathogen'] == pathogen]
    pos = sub[sub['Label'] == 1]
    neg = sub[sub['Label'] == 0]
    
    auc = roc_auc_score(sub['Label'], sub['ACS'])
    stat, pval = mannwhitneyu(pos['ACS'], neg['ACS'], alternative='greater')
    
    pathogen_summary.append({
        'Pathogen': pathogen,
        'Positives (n)': len(pos),
        'Controls (n)': len(neg),
        'Pos Mean ACS': round(pos['ACS'].mean(), 2),
        'Control Mean ACS': round(neg['ACS'].mean(), 2),
        'ROC-AUC': round(auc, 3),
        'p-value': f"{pval:.2e}" if pval < 0.001 else f"{pval:.3f}"
    })

path_df = pd.DataFrame(pathogen_summary)
path_df.to_csv(output_per_pathogen_csv, index=False)
print(f"[SUCCESS] Per-pathogen summary saved to: {output_per_pathogen_csv}")