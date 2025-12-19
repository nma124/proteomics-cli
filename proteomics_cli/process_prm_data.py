#!/usr/bin/env python3
"""
Proteomics PRM Paired Ratio Data Processing

Processes Skyline PRM export data with heavy/light peptide pairs.
Calculates area ratios and performs linear regression per replicate group.

Key concept: Each fragment has 4 replicate groups that are independently analyzed,
then aggregated for final statistics.
"""

import pandas as pd
import numpy as np
from sklearn import linear_model, metrics
import warnings

warnings.filterwarnings("ignore", category=UserWarning)


def get_replicate_number(replicate_name: str) -> int:
    """
    Extract the last number from replicate name.
    
    Example: 'D1_2ul_24' -> 24
    """
    try:
        return int(replicate_name.split('_')[-1])
    except (ValueError, IndexError):
        return 0


def get_replicate_group(replicate_number: int) -> int:
    """
    Assign replicate group based on replicate number.
    
    Groups:
    - Group 1: 7-14
    - Group 2: 17-24
    - Group 3: 27-34
    - Group 4: 37-44
    
    Returns:
        int: Group number (1-4), or 0 if outside known ranges
    """
    if 7 <= replicate_number <= 14:
        return 1
    elif 17 <= replicate_number <= 24:
        return 2
    elif 27 <= replicate_number <= 34:
        return 3
    elif 37 <= replicate_number <= 44:
        return 4
    else:
        # For other ranges, try to infer pattern
        return ((replicate_number - 7) // 10) + 1


def get_dilution_from_replicate(replicate_name: str) -> str:
    """
    Extract dilution identifier from replicate name.
    
    Example: 'D1_2ul_14' -> 'D1'
             'col_test_7' -> 'D0'
    """
    str_to_check = replicate_name.split('_')[0]
    if str_to_check.startswith('col'):
        return 'D0'
    return str_to_check


def get_fragment_ion_with_charge(row) -> str:
    """
    Create fragment ion identifier with charge.
    
    Example: fragment='b5', charge=2 -> 'b5_2+'
    """
    frag = row.get('fragment ion', '')
    pc = row.get('product charge', '')
    try:
        z = int(float(pc))
        return f"{frag}_{z}+"
    except (ValueError, TypeError):
        return str(frag)


def calculate_area_ratio(grouped_df: pd.DataFrame) -> pd.Series:
    """
    Calculate heavy/light area ratio for a paired measurement.
    
    Assumes two precursor m/z values per group:
    - Lower m/z = light
    - Higher m/z = heavy
    
    Ratio = heavy_area / light_area
    """
    df = grouped_df[['precursor mz', 'area']].copy()
    df = df.sort_values(by='precursor mz', ascending=True)
    
    if len(df) < 2:
        return pd.Series({'area_ratio': np.nan, 'light_area': np.nan, 'heavy_area': np.nan})
    
    light_area = df['area'].iloc[0]
    heavy_area = df['area'].iloc[1]
    
    area_ratio = heavy_area / light_area if light_area != 0 else np.inf
    
    return pd.Series({
        'area_ratio': area_ratio,
        'light_area': light_area,
        'heavy_area': heavy_area
    })


def fit_linear_regression(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Fit linear regression and return metrics.
    
    Returns:
        (r2, intercept, slope) or (NaN, NaN, NaN) if insufficient data
    """
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


def qtest(data: np.ndarray, right: bool = True) -> float:
    """
    Dixon's Q-test for outlier detection.
    
    Args:
        data: Array of values
        right: If True, test right tail; if False, test left tail
    
    Returns:
        Q-statistic value
    """
    sorted_data = sorted(data)
    
    if len(sorted_data) < 2:
        return 0.0
    
    if right:
        gap = sorted_data[-1] - sorted_data[-2]
    else:
        gap = sorted_data[1] - sorted_data[0]
    
    range_val = sorted_data[-1] - sorted_data[0]
    q_val = gap / range_val if range_val != 0 else gap
    
    return q_val


def aggregate_regression_metrics(grouped_df: pd.DataFrame) -> pd.Series:
    """
    Aggregate regression metrics across replicate groups.
    
    Computes mean, std, CV, and Q-test for R2, slope, and intercept.
    """
    # Convert to numeric and filter
    r2_vals = pd.to_numeric(grouped_df['r2'], errors='coerce').dropna().values
    r2_vals = r2_vals[np.isfinite(r2_vals)]
    
    slope_vals = pd.to_numeric(grouped_df['slope'], errors='coerce').dropna().values
    slope_vals = slope_vals[np.isfinite(slope_vals)]
    
    intercept_vals = pd.to_numeric(grouped_df['intercept'], errors='coerce').dropna().values
    intercept_vals = intercept_vals[np.isfinite(intercept_vals)]
    
    def safe_stats(vals):
        if len(vals) == 0:
            return np.nan, np.nan, np.nan, np.nan
        mean = float(np.mean(vals))
        std = float(np.std(vals))
        cv = (std / mean * 100) if mean != 0 else np.nan
        q = qtest(vals, right=True) if len(vals) >= 3 else np.nan
        return mean, std, cv, q
    
    r2_mean, r2_std, r2_cv, r2_q = safe_stats(r2_vals)
    slope_mean, slope_std, slope_cv, slope_q = safe_stats(slope_vals)
    intercept_mean, intercept_std, intercept_cv, intercept_q = safe_stats(intercept_vals)
    
    return pd.Series({
        'mean_r2': r2_mean,
        'std_r2': r2_std,
        'cv_r2': r2_cv,
        'q_r2': r2_q,
        'mean_slope': slope_mean,
        'std_slope': slope_std,
        'cv_slope': slope_cv,
        'q_slope': slope_q,
        'mean_intercept': intercept_mean,
        'std_intercept': intercept_std,
        'cv_intercept': intercept_cv,
        'q_intercept': intercept_q,
        'n_groups': len(grouped_df)
    })


def process_prm_data(ms_file: str, concentration_file: str, output_file: str) -> pd.DataFrame:
    """
    Main processing function for paired ratio PRM analysis.
    
    Pipeline:
    1. Load and normalize MS and concentration data
    2. Calculate area ratios (heavy/light) for each measurement
    3. Assign replicate groups based on replicate numbers
    4. Perform linear regression per fragment per replicate group
    5. Aggregate regression metrics across replicate groups
    6. Save results to CSV
    
    Args:
        ms_file: Path to Skyline export CSV
        concentration_file: Path to concentration CSV
        output_file: Path for output results
    
    Returns:
        Final results DataFrame
    """
    print("\n" + "="*70)
    print("STEP 1: Loading Data")
    print("="*70)
    
    # Load MS data
    ms_df = pd.read_csv(ms_file)
    ms_df.columns = ms_df.columns.str.lower().str.strip()
    print(f"Loaded MS data: {len(ms_df)} rows")
    print(f"Unique peptides: {ms_df['peptide'].nunique()}")
    print(f"Unique replicates: {ms_df['replicate'].nunique()}")
    
    # Load concentration data
    conc_df = pd.read_csv(concentration_file)
    conc_df.columns = conc_df.columns.str.strip()
    
    # Melt concentration data to long format
    conc_long = pd.melt(
        conc_df,
        id_vars=['Peptides'],
        var_name='dilution',
        value_name='heavy_conc'
    )
    conc_long['dilution'] = conc_long['dilution'].str.extract(r'(D\d+)')[0]
    conc_long = conc_long.dropna()
    print(f"Loaded concentrations: {len(conc_long)} dilution points")
    
    print("\n" + "="*70)
    print("STEP 2: Calculate Area Ratios")
    print("="*70)
    
    # Add fragment with charge
    ms_df['fragment_charged'] = ms_df.apply(get_fragment_ion_with_charge, axis=1)
    
    # Add replicate metadata
    ms_df['replicate_number'] = ms_df['replicate'].apply(get_replicate_number)
    ms_df['replicate_group'] = ms_df['replicate_number'].apply(get_replicate_group)
    ms_df['dilution'] = ms_df['replicate'].apply(get_dilution_from_replicate)
    
    # Calculate area ratios for each peptide/replicate/fragment
    ratio_df = ms_df.groupby(
        ['peptide', 'protein', 'replicate', 'replicate_number', 'replicate_group', 
         'dilution', 'fragment_charged']
    ).apply(calculate_area_ratio).reset_index()
    
    print(f"Calculated {len(ratio_df)} area ratios")
    print(f"Replicate groups found: {sorted(ratio_df['replicate_group'].unique())}")
    
    # Merge with concentrations
    ratio_df = ratio_df.merge(
        conc_long,
        left_on=['peptide', 'dilution'],
        right_on=['Peptides', 'dilution'],
        how='left'
    )
    
    print("\n" + "="*70)
    print("STEP 3: Linear Regression Per Fragment Per Replicate Group")
    print("="*70)
    
    # Create regression category: peptide_fragment_replicateGroup
    ratio_df['regression_cat'] = (
        ratio_df['peptide'] + '_' + 
        ratio_df['fragment_charged'] + '_' + 
        ratio_df['replicate_group'].astype(str)
    )
    
    # Perform regression for each category
    regression_results = []
    
    for cat in ratio_df['regression_cat'].unique():
        cat_data = ratio_df[ratio_df['regression_cat'] == cat]
        
        # Extract metadata
        peptide = cat_data['peptide'].iloc[0]
        fragment = cat_data['fragment_charged'].iloc[0]
        rep_group = cat_data['replicate_group'].iloc[0]
        protein = cat_data['protein'].iloc[0]
        
        # Get regression data
        x = cat_data['area_ratio'].values
        y = cat_data['heavy_conc'].values
        
        # Fit regression
        r2, intercept, slope = fit_linear_regression(x, y)
        
        regression_results.append({
            'peptide': peptide,
            'protein': protein,
            'fragment': fragment,
            'replicate_group': rep_group,
            'r2': r2,
            'slope': slope,
            'intercept': intercept,
            'n_points': len(cat_data)
        })
    
    regression_df = pd.DataFrame(regression_results)
    print(f"Completed {len(regression_df)} regressions")
    
    print("\n" + "="*70)
    print("STEP 4: Aggregate Metrics Per Fragment")
    print("="*70)
    
    # Aggregate across replicate groups for each fragment
    final_df = regression_df.groupby(['peptide', 'protein', 'fragment']).apply(
        aggregate_regression_metrics
    ).reset_index()
    
    print(f"Final results: {len(final_df)} fragments")
    
    # Save results
    final_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    return final_df
