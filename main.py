#!/usr/bin/env python3
"""
Proteomics PRM Paired Ratio Analysis

Processes Skyline PRM export with heavy/light isotope pairs.
Calculates area ratios and performs linear regression per replicate group.
"""

import sys
import argparse
from pathlib import Path

from proteomics_cli.process_prm_data import process_prm_data


def main():
    """CLI for paired ratio PRM analysis."""
    parser = argparse.ArgumentParser(
        description="Proteomics PRM Paired Ratio Processor - Heavy/Light isotope analysis",
        epilog="Example: python main.py ms_data.csv concentrations.csv -o results.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("ms_file", 
                       help="Path to mass spectrometry data CSV (Skyline PRM export)")
    parser.add_argument("concentration_file", 
                       help="Path to peptide concentration/dilution CSV")
    parser.add_argument("-o", "--output", 
                       default="prm_analysis_output.csv", 
                       help="Output CSV file path (default: prm_analysis_output.csv)")
    parser.add_argument("-p", "--plot",
                       action="store_true",
                       help="Generate calibration curve plots after processing")
    parser.add_argument("--plot-dir",
                       default="plots",
                       help="Directory for plot outputs (default: plots)")
    parser.add_argument("-c", "--config",
                       default=None,
                       help="Path to YAML configuration file (optional)")
    parser.add_argument("--version", 
                       action="version", 
                       version="PRM Paired Ratio Processing v3.0.0")
    
    args = parser.parse_args()
    
    # Validate input files
    ms_path = Path(args.ms_file)
    conc_path = Path(args.concentration_file)
    output_path = Path(args.output)
    
    if not ms_path.exists():
        print(f"[Error] MS data file not found: {ms_path}")
        sys.exit(1)
    
    if not conc_path.exists():
        print(f"[Error] Concentration data file not found: {conc_path}")
        sys.exit(1)
    
    # Header
    print("="*70)
    print("PROTEOMICS PRM PAIRED RATIO ANALYSIS")
    print("="*70)
    print(f"MS data: {ms_path.name}")
    print(f"Concentration data: {conc_path.name}")
    print(f"Output: {output_path}")
    
    try:
        # Process using paired ratio analysis
        result_df = process_prm_data(
            str(ms_path), 
            str(conc_path), 
            str(output_path),
            config_file=args.config
        )
        
        # Success message
        print("\n" + "="*70)
        print("PROCESSING COMPLETE")
        print("="*70)
        print(f"Output: {result_df.shape[0]} rows x {result_df.shape[1]} columns")
        print(f"Saved to: {output_path}")
        
        # Generate plots if requested
        if args.plot:
            from proteomics_cli.plot_results import generate_plots
            generate_plots(
                ms_file=str(ms_path),
                concentration_file=str(conc_path),
                results_csv=str(output_path),
                output_dir=args.plot_dir
            )
        
        print("\nProcessing completed successfully!\n")
        
    except Exception as e:
        print(f"\n[Error] during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
