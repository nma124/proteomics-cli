#!/usr/bin/env python3
"""
JPT Data Processing Module

Processes JPT mass spectrometry and peptide concentration data:
1. Loads Skyline-style PRM export (JPT_1_3_-_HHT_1.2_-_Cytosolic.csv)
2. Loads peptide concentration metadata (JPT1-3_peptide_conc.csv)
3. Extracts experimental conditions from replicate names
4. Aggregates MS intensities by peptide/condition/fragment
5. Merges with concentration data for calibration
6. Computes regression metrics (R², slope, intercept, CV)
7. Generates comprehensive QC report

Based on the heavy_1st processing pipeline.
"""

import pandas as pd
import numpy as np
from sklearn import linear_model, metrics
import pathlib
import warnings

warnings.filterwarnings("ignore", category=UserWarning)


def extract_conditions_from_replicate(replicate_col: pd.Series) -> pd.DataFrame:
    """
    Parse experimental metadata from replicate names.
    
    Example: D1_Soroush_cyto_W_column10_2uL_52 →
      - day: 1
      - condition: cyto
      - wash_type: W
    """
    result = pd.DataFrame()
    result['day'] = replicate_col.str.extract(r'D(\d+)')[0]
    result['condition'] = replicate_col.str.extract(r'(cyto|SDC_Sup)')
    result['wash_type'] = replicate_col.str.extract(r'_(W{1,2})_')
    return result


def get_fragment_ion_with_charge(row):
    """Annotate fragment ion with product charge (e.g., b7_2+)."""
    frag = row.get('fragment ion', '')
    pc = row.get('product charge', '')
    try:
        z = int(float(pc))
        return f"{frag}_{z}+"
    except (ValueError, TypeError):
        return frag


def aggregate_ms_by_condition(ms_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate MS measurements by peptide, condition, and fragment.
    
    Groups precursor ion areas by experimental replicate class,
    computing mean, std, count, and sum for each group.
    """
    # Extract conditions
    conditions = extract_conditions_from_replicate(ms_df['replicate'])
    ms_combined = pd.concat([ms_df, conditions], axis=1)
    
    # Filter for precursor ions (primary quantitative signal)
    precursor = ms_combined[ms_combined['fragment ion'] == 'precursor'].copy()
    
    # Add charged fragment name
    precursor['fragment_ion_charged'] = precursor.apply(get_fragment_ion_with_charge, axis=1)
    
    # Aggregate by peptide, day, condition, wash_type, fragment
    agg = precursor.groupby([
        'peptide', 'protein', 'day', 'condition', 'wash_type', 'fragment_ion_charged'
    ])['area'].agg(['mean', 'std', 'count', 'sum']).reset_index()
    
    agg.columns = [
        'peptide', 'protein', 'day', 'condition', 'wash_type', 'fragment_ion_charged',
        'area_mean', 'area_std', 'area_count', 'area_sum'
    ]
    
    return agg


def load_concentration_data(conc_file: str) -> pd.DataFrame:
    """
    Load peptide concentration metadata.
    
    Expected columns: Peptides, D0, D1, D2, D3, D4, D5, D6, D7 (or subset)
    with concentration values in ng/mL.
    """
    conc_df = pd.read_csv(conc_file)
    
    # Normalize column names
    conc_df.columns = conc_df.columns.str.strip()
    
    # Melt to long format
    conc_long = pd.melt(
        conc_df,
        id_vars=['Peptides'],
        var_name='timepoint',
        value_name='concentration_ng_mL'
    )
    
    # Extract day number
    conc_long['day'] = pd.to_numeric(
        conc_long['timepoint'].str.extract(r'D(\d+)')[0],
        errors='coerce'
    )
    
    conc_long = conc_long.dropna(subset=['day'])
    conc_long['day'] = conc_long['day'].astype(int)
    
    return conc_long


def create_regression_category(row):
    """
    Create a grouping key for regression analysis.
    
    Format: peptide_condition_fragment_ion_charged
    This ensures each unique fragment/charge/condition combo gets its own R².
    """
    return '_'.join([
        row['peptide'],
        row['condition'],
        row['fragment_ion_charged']
    ])


def fit_linear_regression(data: pd.DataFrame, cat: str, x_col: str, y_col: str):
    """
    Fit linear regression for a category.
    
    Returns: (r2, intercept, slope) or (NaN, NaN, NaN) if insufficient data.
    """
    cat_data = data[data['regression_category'] == cat].copy()
    
    if len(cat_data) < 2:
        return (np.nan, np.nan, np.nan)
    
    x = cat_data[x_col].values
    y = cat_data[y_col].values
    
    # Remove infinite and NaN values
    valid_mask = np.isfinite(x) & np.isfinite(y)
    x = x[valid_mask]
    y = y[valid_mask]
    
    if len(x) < 2:
        return (np.nan, np.nan, np.nan)
    
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    
    try:
        model = linear_model.LinearRegression()
        model.fit(x, y)
        y_fit = model.predict(x)
        
        r2 = metrics.r2_score(y, y_fit)
        intercept = model.intercept_.squeeze()
        slope = model.coef_[0].squeeze()
        
        return (r2, intercept, slope)
    except Exception:
        return (np.nan, np.nan, np.nan)


def aggregate_regression_metrics(grouped_df: pd.DataFrame):
    """
    Compute mean, std, CV, and Q-test for regression parameters within a group.
    """
    # Convert to numeric and filter
    r2_vals = pd.to_numeric(grouped_df['R2'], errors='coerce').dropna().values.astype(float)
    r2_vals = r2_vals[np.isfinite(r2_vals)]
    
    slope_vals = pd.to_numeric(grouped_df['slope'], errors='coerce').dropna().values.astype(float)
    slope_vals = slope_vals[np.isfinite(slope_vals)]
    
    intercept_vals = pd.to_numeric(grouped_df['intercept'], errors='coerce').dropna().values.astype(float)
    intercept_vals = intercept_vals[np.isfinite(intercept_vals)]
    
    def safe_stats(vals):
        if len(vals) == 0:
            return np.nan, np.nan, np.nan
        mean = float(np.mean(vals))
        std = float(np.std(vals))
        cv = std / mean if mean != 0 else np.nan
        return mean, std, cv
    
    r2_mean, r2_std, r2_cv = safe_stats(r2_vals)
    slope_mean, slope_std, slope_cv = safe_stats(slope_vals)
    intercept_mean, intercept_std, intercept_cv = safe_stats(intercept_vals)
    
    return pd.Series({
        'mean_r2': r2_mean, 'std_r2': r2_std, 'cv_r2': r2_cv,
        'mean_slope': slope_mean, 'std_slope': slope_std, 'cv_slope': slope_cv,
        'mean_intercept': intercept_mean, 'std_intercept': intercept_std, 'cv_intercept': intercept_cv,
    })


def process_jpt_data(ms_file: str, conc_file: str, output_file: str):
    """
    Main processing function for JPT data.
    
    Args:
        ms_file: Path to Skyline PRM export CSV
        conc_file: Path to peptide concentration CSV
        output_file: Path for output results CSV
    
    Returns:
        DataFrame with per-fragment regression results and QC metrics
    """
    print("📂 Loading JPT data files...")
    
    # Load MS data
    ms_df = pd.read_csv(ms_file)
    ms_df.columns = ms_df.columns.str.lower().str.strip()
    print(f"   ✓ MS data: {len(ms_df)} measurements")
    
    # Load concentration data
    conc_long = load_concentration_data(conc_file)
    print(f"   ✓ Concentration data: {len(conc_long)} timepoint measurements")
    
    # Aggregate MS data
    print("\n🔧 Processing MS data...")
    ms_agg = aggregate_ms_by_condition(ms_df)
    print(f"   ✓ Aggregated: {len(ms_agg)} peptide-condition-fragment combinations")
    
    # Merge with concentration data
    print("\n🔧 Merging with concentration data...")
    ms_agg['day'] = pd.to_numeric(ms_agg['day'], errors='coerce')
    
    conc_merge = conc_long.rename(columns={'Peptides': 'peptide'}).copy()
    conc_merge['day'] = pd.to_numeric(conc_merge['day'], errors='coerce')
    
    merged = ms_agg.merge(
        conc_merge[['peptide', 'day', 'concentration_ng_mL']],
        on=['peptide', 'day'],
        how='left'
    )
    print(f"   ✓ Merged: {len(merged)} rows with concentration data")
    
    # Remove rows with missing concentration
    merged = merged.dropna(subset=['concentration_ng_mL'])
    print(f"   ✓ After filtering: {len(merged)} rows with valid concentrations")
    
    # Create regression category
    merged['regression_category'] = merged.apply(create_regression_category, axis=1)
    
    # Fit regressions for each category
    print("\n🔧 Computing regression metrics...")
    regression_categories = merged['regression_category'].unique()
    regression_results = []
    
    for cat in regression_categories:
        r2, intercept, slope = fit_linear_regression(
            merged, cat, 'concentration_ng_mL', 'area_mean'
        )
        regression_results.append({
            'regression_category': cat,
            'R2': r2,
            'intercept': intercept,
            'slope': slope
        })
    
    regression_df = pd.DataFrame(regression_results)
    print(f"   ✓ Fitted {len(regression_df)} fragment/condition regressions")
    
    # Merge regression metrics back into main data
    final_df = merged.merge(regression_df, on='regression_category', how='left')
    
    # Aggregate regression metrics by fragment/condition
    print("\n🔧 Aggregating QC metrics...")
    agg_qc = final_df.groupby('regression_category').apply(
        aggregate_regression_metrics
    ).reset_index()
    
    final_df = final_df.merge(agg_qc, on='regression_category', how='left')
    
    # Save results
    print(f"\n💾 Writing output to {pathlib.Path(output_file).name}")
    final_df.to_csv(output_file, index=False)
    print(f"   ✓ Saved {len(final_df)} rows")
    
    return final_df
