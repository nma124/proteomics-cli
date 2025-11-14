#!/usr/bin/env python3
"""
Unified PRM Data Processor

General-purpose interface that handles multiple PRM input formats:
- Heavy/Light paired peptides (area ratio analysis)
- Single peptide intensity analysis (direct calibration)

The processor automatically detects the format and routes to appropriate logic.
"""

import pandas as pd
import numpy as np
from sklearn import linear_model, metrics
import warnings

warnings.filterwarnings("ignore", category=UserWarning)


class PRMDataProcessor:
    """
    General-purpose PRM data processor.
    
    Handles different experimental designs by detecting key features:
    - Paired measurements (light/heavy) vs single measurements
    - Different replicate naming conventions
    - Varying experimental conditions
    """
    
    def __init__(self, ms_file: str, concentration_file: str):
        """Load and inspect input files."""
        self.ms_file = ms_file
        self.concentration_file = concentration_file
        
        # Load data
        self.ms_df = pd.read_csv(ms_file)
        self.ms_df.columns = self.ms_df.columns.str.lower().str.strip()
        
        self.conc_df = pd.read_csv(concentration_file)
        self.conc_df.columns = self.conc_df.columns.str.strip()
        
        # Detect format
        self.format_type = self._detect_format()
        self.has_pairs = self._check_for_pairs()
        self.has_conditions = self._check_for_conditions()
        
    def _detect_format(self) -> dict:
        """
        Inspect data structure and return format characteristics.
        
        Returns dict with:
        - has_light_heavy_pairs: bool
        - has_multiple_conditions: bool
        - dilution_identifier: str (column name or replicate pattern)
        - concentration_columns: list
        """
        format_info = {
            'has_light_heavy_pairs': False,
            'has_multiple_conditions': False,
            'dilution_identifier': None,
            'concentration_columns': [],
            'replicate_pattern': None
        }
        
        # Check concentration file for day/dilution columns
        conc_cols = [col for col in self.conc_df.columns if col.startswith('D')]
        format_info['concentration_columns'] = conc_cols
        
        # Check MS file for pairs (multiple precursor m/z per peptide/replicate/fragment)
        if 'precursor mz' in self.ms_df.columns:
            grouped = self.ms_df.groupby(
                ['peptide', 'replicate', 'fragment ion']
            )['precursor mz'].nunique()
            format_info['has_light_heavy_pairs'] = (grouped > 1).any()
        
        # Check for experimental conditions in replicate names
        if 'replicate' in self.ms_df.columns:
            sample_rep = str(self.ms_df['replicate'].iloc[0])
            # Look for condition indicators
            has_conditions = any(cond in sample_rep for cond in ['cyto', 'SDC_Sup', 'W_', 'WW_'])
            format_info['has_multiple_conditions'] = has_conditions
            format_info['replicate_pattern'] = sample_rep
        
        return format_info
    
    def _check_for_pairs(self) -> bool:
        """Check if data contains light/heavy pairs."""
        return self.format_type['has_light_heavy_pairs']
    
    def _check_for_conditions(self) -> bool:
        """Check if data has multiple experimental conditions."""
        return self.format_type['has_multiple_conditions']
    
    def get_processing_mode(self) -> str:
        """
        Determine which processing pipeline to use.
        
        Returns:
            'paired_ratio': Heavy/light area ratio analysis
            'single_intensity': Direct intensity calibration with conditions
            'single_basic': Direct intensity calibration without conditions
        """
        if self.has_pairs:
            return 'paired_ratio'
        elif self.has_conditions:
            return 'single_intensity'
        else:
            return 'single_basic'
    
    def process(self, output_file: str) -> pd.DataFrame:
        """
        Process data using appropriate pipeline.
        
        Args:
            output_file: Path to save results
            
        Returns:
            Processed DataFrame with regression metrics
        """
        mode = self.get_processing_mode()
        
        print(f"🔧 Processing mode: {mode}")
        print(f"   • Light/Heavy pairs: {self.has_pairs}")
        print(f"   • Multiple conditions: {self.has_conditions}")
        print(f"   • Concentration columns: {len(self.format_type['concentration_columns'])}")
        
        if mode == 'paired_ratio':
            return self._process_paired_ratio(output_file)
        elif mode == 'single_intensity':
            return self._process_single_with_conditions(output_file)
        else:
            return self._process_single_basic(output_file)
    
    def _process_paired_ratio(self, output_file: str) -> pd.DataFrame:
        """Process heavy/light paired peptides (Heavy_1st format)."""
        from proteomics_cli.process_prm_data import process_prm_data
        return process_prm_data(self.ms_file, self.concentration_file, output_file)
    
    def _process_single_with_conditions(self, output_file: str) -> pd.DataFrame:
        """Process single peptides with experimental conditions (JPT format)."""
        from proteomics_cli.process_jpt_data import process_jpt_data
        return process_jpt_data(self.ms_file, self.concentration_file, output_file)
    
    def _process_single_basic(self, output_file: str) -> pd.DataFrame:
        """Process single peptides without conditions (basic format)."""
        # This is a simplified version for basic single-peptide analysis
        print("⚠️  Using basic single-peptide processing (no conditions)")
        print("    If you have experimental conditions, update replicate naming.")
        
        # Use JPT processor but it will treat all as one condition
        from proteomics_cli.process_jpt_data import process_jpt_data
        return process_jpt_data(self.ms_file, self.concentration_file, output_file)
    
    def print_format_summary(self):
        """Print detected format information."""
        print("\n" + "="*70)
        print("📋 DATA FORMAT SUMMARY")
        print("="*70)
        print(f"MS measurements: {len(self.ms_df)} rows")
        print(f"Peptides: {self.ms_df['peptide'].nunique()}")
        print(f"Replicates: {self.ms_df['replicate'].nunique()}")
        
        print(f"\n🔍 Format Detection:")
        print(f"   • Light/Heavy pairs: {'Yes ✓' if self.has_pairs else 'No ✗'}")
        print(f"   • Multiple conditions: {'Yes ✓' if self.has_conditions else 'No ✗'}")
        print(f"   • Dilution points: {len(self.format_type['concentration_columns'])}")
        print(f"   • Sample replicate: {self.format_type['replicate_pattern']}")
        
        mode = self.get_processing_mode()
        print(f"\n🎯 Selected processing mode: {mode.upper().replace('_', ' ')}")
        
        if mode == 'paired_ratio':
            print("   → Will calculate heavy/light area ratios")
            print("   → Regression: area_ratio ~ heavy_concentration")
        elif mode == 'single_intensity':
            print("   → Will aggregate intensities by condition")
            print("   → Regression: intensity ~ concentration (per condition)")
        else:
            print("   → Will process as single peptide intensities")
            print("   → Regression: intensity ~ concentration")


def process_prm_unified(ms_file: str, concentration_file: str, output_file: str) -> pd.DataFrame:
    """
    Unified entry point for PRM data processing.
    
    Automatically detects format and processes accordingly.
    
    Args:
        ms_file: Path to mass spectrometry CSV (Skyline export)
        concentration_file: Path to peptide concentration CSV
        output_file: Path for output results
        
    Returns:
        Processed DataFrame with regression analysis
    """
    processor = PRMDataProcessor(ms_file, concentration_file)
    processor.print_format_summary()
    
    print("\n" + "-"*70)
    result_df = processor.process(output_file)
    
    return result_df
