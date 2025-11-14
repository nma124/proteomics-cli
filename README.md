# Proteomics PRM CLI

A unified command-line tool for automated proteomics PRM (Parallel Reaction Monitoring) data analysis with automatic format detection.

## Features

- 🔍 **Automatic Format Detection** - Intelligently detects data format and processing mode
- 📊 **Multiple Format Support** - Handles heavy/light paired peptides and single intensity measurements
- 📈 **Statistical Analysis** - Linear regression with R², CV, and Q-test metrics
- ⚡ **Fast Processing** - Optimized pipeline for large datasets
- 🎯 **No Configuration Required** - Just point to your data files

## Installation

```bash
# Clone the repository
git clone https://github.com/nma124/proteomics-cli.git
cd proteomics-cli

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

```bash
python main.py <ms_data.csv> <concentrations.csv> -o results.csv
```

### Example

```bash
python main.py \
  data/skyline_export.csv \
  data/peptide_concentrations.csv \
  -o analysis_results.csv
```

## Usage

```
python main.py [-h] [--version] [-o OUTPUT] ms_file concentration_file

Positional arguments:
  ms_file              Path to mass spectrometry data CSV (Skyline PRM export)
  concentration_file   Path to peptide concentration/dilution CSV

Optional arguments:
  -h, --help          Show help message
  --version           Show version
  -o, --output        Output CSV file path (default: prm_analysis_output.csv)
```

## Input File Requirements

### MS Data File (Skyline Export)
Required columns:
- `Peptide`
- `Protein`
- `Replicate`
- `Precursor Mz`
- `Precursor Charge`
- `Product Mz`
- `Product Charge`
- `Fragment Ion`
- `Area`

### Concentration File
Required columns:
- `Peptides`: Peptide sequences (must match MS file)
- `D0 (ng/mL)`, `D1 (ng/mL)`, etc.: Concentration values for each dilution point

## Processing Modes

The tool automatically selects the appropriate processing mode:

### 1. Paired Ratio Mode
**Triggered when:** Light/heavy isotope pairs detected

**Processing:**
- Calculates area ratios (heavy/light)
- Fits linear regression: `area_ratio ~ concentration`
- Computes per-fragment R², CV, and Q-test statistics

### 2. Single Intensity Mode
**Triggered when:** No pairs, but experimental conditions detected

**Processing:**
- Extracts conditions from replicate names
- Aggregates intensities by peptide/condition/fragment
- Fits regression: `intensity ~ concentration` per condition

### 3. Single Basic Mode
**Triggered when:** Simple single-peptide format

**Processing:**
- Direct intensity vs concentration regression
- Simplified QC metrics

## Output

Results are saved as a CSV file containing:
- Regression metrics (R², slope, intercept)
- Aggregated statistics (mean, std, CV)
- QC metrics (Q-test values)
- Per-fragment analysis results

## Examples

### Heavy/Light Paired Peptides
```bash
python main.py \
  data/heavy_light_export.csv \
  data/dilution_series.csv \
  -o paired_results.csv
```

### Single Intensity with Conditions
```bash
python main.py \
  data/jpt_cytosolic.csv \
  data/jpt_concentrations.csv \
  -o jpt_results.csv
```

## Requirements

- Python 3.8+
- pandas >= 1.3.0
- numpy >= 1.21.0
- scikit-learn >= 1.0.0

