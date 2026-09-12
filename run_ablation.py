# ========= Import Libraries ===========
import os
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

# ========= Define Paths =============
current_dir = Path(__file__).parent
results_csv = current_dir / "validation" / "output" / "validation_results_scored.csv"
output_summary_csv = current_dir / "validation" / "output" / "ablation_results_summary.csv"
output_plot = current_dir / "validation" / "output" / "ablation_analysis.png"

# ======== Default Weights for MemVax-Py Modules ==========
# Mod1: Antigenicity (30%), Mod2: Topology (25%), Mod3: HostSim (15%), Mod4: Conservation (20%), Mod5: UsableRegion (10%)
DEFAULT_WEIGHTS = {
    'Mod1_Antigenicity': 0.30,
    'Mod2_Topology': 0.25,
    'Mod3_HostSim': 0.15,
    'Mod4_Conservation': 0.20,
    'Mod5_UsableRegion': 0.10
}

def renormalize_weights(omitted_module, weights):
    """
    Removes one module and renormalizes the remaining weights so they sum to 1.0 (100%).
    """
    remaining_weights = {k: v for k, v in weights.items() if k != omitted_module}
    total_remaining = sum(remaining_weights.values())
    renormalized = {k: v / total_remaining for k, v in remaining_weights.items()}
    return renormalized

def calculate_ablated_acs(df, active_weights):
    """
    Calculates composite ACS given a set of active, normalized module weights.
    Each module column in the CSV is scored out of 25; multiplying by 4 scales to 0-100.
    """
    acs_series = pd.Series(0.0, index=df.index)
    for mod_col, w in active_weights.items():
        acs_series += (df[mod_col] * 4.0) * w
    return acs_series

def run_ablation_study():
    print("==================================================")
    print("     MEMVAX-PY ABLATION ANALYSIS EXECUTION       ")
    print("==================================================")
    
    if not results_csv.exists():
        raise FileNotFoundError(f"Could not find input CSV at {results_csv}. Run validation first.")

    df = pd.read_csv(results_csv)
    print(f"Loaded validation dataset: {len(df)} total proteins ({sum(df['Label']==1)} Positives, {sum(df['Label']==0)} Controls)\n")

    # Compute baseline full model AUC dynamically
    baseline_acs = calculate_ablated_acs(df, DEFAULT_WEIGHTS)
    baseline_auc = roc_auc_score(df['Label'], baseline_acs)

    # Define ablation scenarios
    scenarios = {
        'Full Model (A+T+H+C+U)': None,
        'Without Antigenicity (-A)': 'Mod1_Antigenicity',
        'Without Topology (-T)': 'Mod2_Topology',
        'Without Host Similarity (-H)': 'Mod3_HostSim',
        'Without Conservation (-C)': 'Mod4_Conservation',
        'Without Usable Region (-U)': 'Mod5_UsableRegion'
    }

    ablation_summary = []

    for name, omit_mod in scenarios.items():
        if omit_mod is None:
            active_weights = DEFAULT_WEIGHTS.copy()
            weight_str = "A:30%, T:25%, H:15%, C:20%, U:10%"
        else:
            active_weights = renormalize_weights(omit_mod, DEFAULT_WEIGHTS)
            weight_str = ", ".join([f"{k.split('_')[0]}:{v*100:.2f}%" for k, v in active_weights.items()])

        # Compute ACS for this scenario
        acs_scores = calculate_ablated_acs(df, active_weights)
        
        # Calculate ROC-AUC & Mann-Whitney U test
        auc = roc_auc_score(df['Label'], acs_scores)
        pos_scores = acs_scores[df['Label'] == 1]
        neg_scores = acs_scores[df['Label'] == 0]
        stat, pval = mannwhitneyu(pos_scores, neg_scores, alternative='greater')

        # Compute dynamic delta against baseline
        delta_auc = auc - baseline_auc

        ablation_summary.append({
            'Model Condition': name,
            'Omitted Module': omit_mod if omit_mod else 'None',
            'Renormalized Weights': weight_str,
            'Positives Mean ACS': round(pos_scores.mean(), 2),
            'Controls Mean ACS': round(neg_scores.mean(), 2),
            'ROC-AUC': round(auc, 4),
            'Delta AUC': round(delta_auc, 4),
            'p-value': f"{pval:.2e}" if pval < 0.001 else f"{pval:.4f}"
        })

    # Convert to DataFrame
    summary_df = pd.DataFrame(ablation_summary)
    
    # Save results
    summary_df.to_csv(output_summary_csv, index=False)
    print("Ablation Results Summary:")
    print("--------------------------------------------------")
    print(summary_df[['Model Condition', 'ROC-AUC', 'Delta AUC', 'p-value']].to_string(index=False))
    print("--------------------------------------------------")

    # Generate Visualization
    plt.figure(figsize=(10, 6))
    colors = ['#2b5c8f' if x == 'Full Model (A+T+H+C+U)' else '#d95f02' if 'Antigenicity' in x else '#7570b3' for x in summary_df['Model Condition']]
    
    bars = plt.barh(summary_df['Model Condition'], summary_df['ROC-AUC'], color=colors, edgecolor='black', alpha=0.85)
    plt.axvline(x=baseline_auc, color='black', linestyle='--', linewidth=1.5, label=f'Full Model Baseline ({baseline_auc:.3f})')
    plt.axvline(x=0.5000, color='gray', linestyle=':', linewidth=1.0, label='Random Chance (0.500)')
    
    plt.xlim(0.40, 0.90)
    plt.xlabel('ROC-AUC Score', fontsize=12, fontweight='bold')
    plt.ylabel('Ablation Condition', fontsize=12, fontweight='bold')
    plt.title('MemVax-Py Leave-One-Module-Out Ablation Study', fontsize=14, fontweight='bold', pad=15)
    
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.005, bar.get_y() + bar.get_height()/2, f'{width:.4f}', 
                 va='center', ha='left', fontsize=10, fontweight='bold')

    plt.legend(loc='lower right')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()

    print(f"\n[SUCCESS] Summary CSV exported to: {output_summary_csv}")
    print(f"[SUCCESS] Ablation plot saved to: {output_plot}\n")

if __name__ == "__main__":
    run_ablation_study()