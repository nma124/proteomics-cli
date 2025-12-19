#!/usr/bin/env python3
"""
Plotting Module for PRM Paired Ratio Analysis

Generates calibration curve plots showing:
- Area ratio vs concentration for each fragment
- Separate colors for each replicate group (1-4)
- Individual regression lines per group with R² annotations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path
from typing import Dict, List

# Use non-interactive backend for server environments
matplotlib.use('Agg')


def create_output_directories(output_dir: Path) -> Dict[str, Path]:
    """
    Create directory structure for plot outputs.
    
    Returns:
        Dictionary with paths for 'fragments' and 'summary' directories
    """
    fragments_dir = output_dir / "fragments"
    summary_dir = output_dir / "summary"
    
    fragments_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)
    
    return {
        'fragments': fragments_dir,
        'summary': summary_dir
    }


def plot_fragment_calibration(
    data: pd.DataFrame,
    fragment_name: str,
    peptide_name: str,
    output_path: Path,
    regression_data: pd.DataFrame
) -> None:
    """
    Create calibration curve plot for a single fragment showing all replicate groups.
    
    Args:
        data: DataFrame with area_ratio, heavy_conc, replicate_group columns
        fragment_name: Fragment identifier (e.g., "b5_2+")
        peptide_name: Peptide sequence
        output_path: Path to save the plot
        regression_data: DataFrame with r2, slope, intercept per replicate_group
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Color scheme for 4 replicate groups
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']  # blue, orange, green, red
    group_names = {1: 'Group 1 (7-14)', 2: 'Group 2 (17-24)', 
                   3: 'Group 3 (27-34)', 4: 'Group 4 (37-44)'}
    
    # Plot each replicate group
    for idx, group in enumerate(sorted(data['replicate_group'].unique())):
        if group == 0:  # Skip invalid groups
            continue
            
        group_data = data[data['replicate_group'] == group]
        
        # Get regression parameters for this group
        reg_info = regression_data[regression_data['replicate_group'] == group]
        if len(reg_info) == 0:
            continue
            
        r2 = reg_info['r2'].values[0]
        slope = reg_info['slope'].values[0]
        intercept = reg_info['intercept'].values[0]
        
        # Skip if invalid regression
        if not np.isfinite([r2, slope, intercept]).all():
            continue
        
        # Scatter plot
        color = colors[idx % len(colors)]
        ax.scatter(
            group_data['area_ratio'], 
            group_data['heavy_conc'],
            alpha=0.6,
            s=80,
            color=color,
            label=f"{group_names.get(group, f'Group {group}')} (R²={r2:.3f})",
            edgecolors='black',
            linewidth=0.5
        )
        
        # Regression line
        x_range = np.linspace(
            group_data['area_ratio'].min(),
            group_data['area_ratio'].max(),
            100
        )
        y_fit = slope * x_range + intercept
        ax.plot(x_range, y_fit, '--', color=color, linewidth=2, alpha=0.8)
    
    # Labels and title
    ax.set_xlabel('Area Ratio (Heavy/Light)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Heavy Concentration (ng/mL)', fontsize=12, fontweight='bold')
    ax.set_title(
        f'{peptide_name}\n{fragment_name}',
        fontsize=14,
        fontweight='bold',
        pad=20
    )
    
    # Legend
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=10)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Tight layout
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_summary_plot(
    summary_df: pd.DataFrame,
    output_path: Path
) -> None:
    """
    Create summary visualization showing R² distribution across all fragments.
    
    Args:
        summary_df: DataFrame with aggregated metrics per fragment
        output_path: Path to save the plot
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: R² distribution
    r2_values = summary_df['mean_r2'].dropna()
    ax1.hist(r2_values, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axvline(r2_values.mean(), color='red', linestyle='--', linewidth=2, 
                label=f'Mean R² = {r2_values.mean():.3f}')
    ax1.set_xlabel('Mean R²', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of Mean R² Across Fragments', fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: CV distribution
    cv_values = summary_df['cv_r2'].dropna()
    ax2.hist(cv_values, bins=20, color='darkorange', edgecolor='black', alpha=0.7)
    ax2.axvline(cv_values.median(), color='red', linestyle='--', linewidth=2,
                label=f'Median CV = {cv_values.median():.1f}%')
    ax2.set_xlabel('CV of R² (%)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title('Distribution of R² Coefficient of Variation', fontsize=13, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_plots(
    ms_file: str,
    concentration_file: str,
    results_csv: str,
    output_dir: str = "plots"
) -> None:
    """
    Generate all calibration plots from processed results.
    
    This function re-loads the raw data and results to create plots showing:
    1. Individual fragment calibration curves (one per fragment)
    2. Summary statistics visualization
    
    Args:
        ms_file: Path to original MS data CSV
        concentration_file: Path to concentration CSV
        results_csv: Path to processed results CSV
        output_dir: Directory to save plots (default: "plots")
    """
    print("\n" + "="*70)
    print("GENERATING CALIBRATION PLOTS")
    print("="*70)
    
    # Create output directories
    output_path = Path(output_dir)
    dirs = create_output_directories(output_path)
    
    # Load results for summary
    results_df = pd.read_csv(results_csv)
    print(f"Loaded results: {len(results_df)} fragments")
    
    # Re-process data to get individual points for plotting
    # (We need the raw area_ratio data, not just aggregated statistics)
    from proteomics_cli.process_prm_data import (
        get_replicate_number, get_replicate_group, get_dilution_from_replicate,
        get_fragment_ion_with_charge, calculate_area_ratio
    )
    
    # Load MS data
    ms_df = pd.read_csv(ms_file)
    ms_df.columns = ms_df.columns.str.lower().str.strip()
    
    # Load concentration data
    conc_df = pd.read_csv(concentration_file)
    conc_df.columns = conc_df.columns.str.strip()
    conc_long = pd.melt(
        conc_df,
        id_vars=['Peptides'],
        var_name='dilution',
        value_name='heavy_conc'
    )
    conc_long['dilution'] = conc_long['dilution'].str.extract(r'(D\d+)')[0]
    conc_long = conc_long.dropna()
    
    # Process MS data
    ms_df['fragment_charged'] = ms_df.apply(get_fragment_ion_with_charge, axis=1)
    ms_df['replicate_number'] = ms_df['replicate'].apply(get_replicate_number)
    ms_df['replicate_group'] = ms_df['replicate_number'].apply(get_replicate_group)
    ms_df['dilution'] = ms_df['replicate'].apply(get_dilution_from_replicate)
    
    # Calculate ratios
    ratio_df = ms_df.groupby(
        ['peptide', 'protein', 'replicate', 'replicate_number', 'replicate_group', 
         'dilution', 'fragment_charged']
    ).apply(calculate_area_ratio).reset_index()
    
    # Merge with concentrations
    ratio_df = ratio_df.merge(
        conc_long,
        left_on=['peptide', 'dilution'],
        right_on=['Peptides', 'dilution'],
        how='left'
    )
    
    # Create regression lookup
    # Re-calculate regressions for each fragment/group combination
    from proteomics_cli.process_prm_data import fit_linear_regression
    
    ratio_df['regression_cat'] = (
        ratio_df['peptide'] + '_' + 
        ratio_df['fragment_charged'] + '_' + 
        ratio_df['replicate_group'].astype(str)
    )
    
    regression_results = []
    for cat in ratio_df['regression_cat'].unique():
        cat_data = ratio_df[ratio_df['regression_cat'] == cat]
        peptide = cat_data['peptide'].iloc[0]
        fragment = cat_data['fragment_charged'].iloc[0]
        rep_group = cat_data['replicate_group'].iloc[0]
        
        x = cat_data['area_ratio'].values
        y = cat_data['heavy_conc'].values
        r2, intercept, slope = fit_linear_regression(x, y)
        
        regression_results.append({
            'peptide': peptide,
            'fragment': fragment,
            'replicate_group': rep_group,
            'r2': r2,
            'slope': slope,
            'intercept': intercept
        })
    
    regression_df = pd.DataFrame(regression_results)
    
    # Generate individual fragment plots
    print(f"\nGenerating fragment plots...")
    fragments = ratio_df.groupby(['peptide', 'fragment_charged'])
    
    plot_count = 0
    for (peptide, fragment), group_data in fragments:
        # Get regression data for this fragment
        frag_regression = regression_df[
            (regression_df['peptide'] == peptide) & 
            (regression_df['fragment'] == fragment)
        ]
        
        if len(frag_regression) == 0:
            continue
        
        # Create safe filename
        safe_fragment = fragment.replace('+', 'plus').replace('/', '_')
        filename = f"{peptide}_{safe_fragment}.png"
        output_file = dirs['fragments'] / filename
        
        # Generate plot
        plot_fragment_calibration(
            data=group_data,
            fragment_name=fragment,
            peptide_name=peptide,
            output_path=output_file,
            regression_data=frag_regression
        )
        
        plot_count += 1
    
    print(f"Created {plot_count} fragment plots in: {dirs['fragments']}")
    
    # Generate summary plot
    print(f"\nGenerating summary plot...")
    summary_file = dirs['summary'] / "r2_summary.png"
    create_summary_plot(results_df, summary_file)
    print(f"Created summary plot in: {dirs['summary']}")
    
    print("\n" + "="*70)
    print("PLOTTING COMPLETE")
    print("="*70)
    print(f"Total plots generated: {plot_count + 1}")
    print(f"Output directory: {output_path.absolute()}")
