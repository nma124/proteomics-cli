#!/usr/bin/env python3
"""
Unified Proteomics PRM Data Processing CLI

General-purpose interface that automatically detects and processes:
- Heavy/Light paired peptide formats
- Single peptide intensity formats (with/without conditions)
- Any future format variations
"""

import sys
import argparse
from pathlib import Path

from proteomics_cli.processor import process_prm_unified


def main():
    """Unified CLI with automatic format detection."""
    parser = argparse.ArgumentParser(
        description="Unified Proteomics PRM Data Processor - Handles any Skyline export format",
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
    parser.add_argument("--version", 
                       action="version", 
                       version="Unified PRM Processing v2.0.0")
    
    args = parser.parse_args()
    
    # Validate input files
    ms_path = Path(args.ms_file)
    conc_path = Path(args.concentration_file)
    output_path = Path(args.output)
    
    if not ms_path.exists():
        print(f"❌ Error: MS data file not found: {ms_path}")
        sys.exit(1)
    
    if not conc_path.exists():
        print(f"❌ Error: Concentration data file not found: {conc_path}")
        sys.exit(1)
    
    # Header
    print("="*70)
    print("UNIFIED PROTEOMICS PRM DATA PROCESSING")
    print("="*70)
    print(f"📁 MS data: {ms_path.name}")
    print(f"📁 Concentration data: {conc_path.name}")
    print(f"💾 Output: {output_path}")
    
    try:
        # Process using unified processor (auto-detects format)
        result_df = process_prm_unified(str(ms_path), str(conc_path), str(output_path))
        
        # Success message
        print("\n" + "="*70)
        print("✅ PROCESSING COMPLETE")
        print("="*70)
        print(f"📊 Output: {result_df.shape[0]} rows × {result_df.shape[1]} columns")
        print(f"💾 Saved to: {output_path}")
        print("\n✨ Processing completed successfully!\n")
        
    except Exception as e:
        print(f"\n❌ Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
