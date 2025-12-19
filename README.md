# Proteomics PRM CLI

A command-line tool for proteomics PRM (Parallel Reaction Monitoring) data analysis with heavy/light isotope pairs.

## Features

- **Paired Ratio Analysis** - Heavy/light isotope quantification
- **Replicate Group Regression** - Independent R² calculations per replicate group (4 groups)
- **Statistical Analysis** - Linear regression with R², CV, and Q-test metrics
- **Calibration Plots** - Automated generation of scatter plots with regression lines
- **Fast Processing** - Optimized pipeline for large datasets

## Installation

```bash
# Clone the repository
git clone https://github.com/nma124/proteomics-cli.git
cd proteomics-cli

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### Basic Processing
```bash
python main.py <ms_data.csv> <concentrations.csv> -o results.csv
```

### With Calibration Plots
```bash
python main.py <ms_data.csv> <concentrations.csv> -o results.csv -p
```

### Example
```bash
python main.py \
  input/heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv \
  input/peptide_dilution_conc_peggy.csv \
  -o results.csv \
  -p
```

## Usage

```
python main.py [-h] [--version] [-o OUTPUT] [-p] [--plot-dir PLOT_DIR] [-c CONFIG] ms_file concentration_file

Positional arguments:
  ms_file              Path to mass spectrometry data CSV (Skyline PRM export)
  concentration_file   Path to peptide concentration/dilution CSV

Optional arguments:
  -h, --help           Show help message
  --version            Show version
  -o, --output OUTPUT  Output CSV file path (default: prm_analysis_output.csv)
  -p, --plot           Generate calibration curve plots
  --plot-dir PLOT_DIR  Directory for plot outputs (default: plots)
  -c, --config CONFIG  Path to YAML configuration file (optional)
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

### Replicate Naming Convention
Replicate names should follow the pattern: `D{dilution}_{description}_{number}`

Example: `D1_2ul_24`
- `D1` = Dilution level
- `2ul` = Description (optional)
- `24` = Replicate number (used for grouping)

**Replicate Groups:**
- Group 1: Numbers 7-14
- Group 2: Numbers 17-24
- Group 3: Numbers 27-34
- Group 4: Numbers 37-44

## Processing Pipeline

### 1. Data Loading
- Loads Skyline PRM export and concentration data
- Normalizes column names
- Validates required columns

### 2. Area Ratio Calculation
- Groups by peptide/replicate/fragment
- Identifies light (lower m/z) and heavy (higher m/z) pairs
- Calculates ratio: `heavy_area / light_area`

### 3. Replicate Grouping
- Extracts replicate number from replicate name
- Assigns to one of 4 groups based on number ranges
- Each group analyzed independently

### 4. Linear Regression
- **Per fragment per replicate group**: Individual R², slope, intercept
- Fits: `heavy_concentration ~ area_ratio`
- Filters out infinite/NaN values

### 5. Aggregation
- Combines results across 4 replicate groups
- Computes: mean, std, CV, Q-test for R², slope, intercept
- One row per fragment in output

## Output

### CSV Results File
One row per fragment with columns:
- **Identifiers**: `peptide`, `protein`, `fragment`
- **Regression Metrics** (aggregated across 4 groups):
  - `mean_r2`, `std_r2`, `cv_r2`, `q_r2`
  - `mean_slope`, `std_slope`, `cv_slope`, `q_slope`
  - `mean_intercept`, `std_intercept`, `cv_intercept`, `q_intercept`
- **Metadata**: `n_groups` (number of replicate groups with valid data)

### Plots (if `-p` flag used)
Generated in `plots/` directory:
```
plots/
├── fragments/
│   ├── PEPTIDE_fragment_1plus.png
│   ├── PEPTIDE_fragment_2plus.png
│   └── ...
└── summary/
    └── r2_summary.png
```

**Fragment plots** show:
- Scatter points colored by replicate group (4 colors)
- Individual regression line per group
- R² value for each group in legend

**Summary plot** shows:
- Distribution of mean R² values across fragments
- Distribution of CV values

## Examples

### Basic Analysis
```bash
python main.py \
  input/heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv \
  input/peptide_dilution_conc_peggy.csv \
  -o results.csv
```

### With Plots
```bash
python main.py \
  input/heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv \
  input/peptide_dilution_conc_peggy.csv \
  -o results.csv \
  -p
```

### Custom Output Locations
```bash
python main.py \
  input/heavy_1st_expanded_D0_op_AHS_3+_1+_2+_combined_w0.csv \
  input/peptide_dilution_conc_peggy.csv \
  -o output/my_results.csv \
  -p \
  --plot-dir output/my_plots
```

## Requirements

- Python 3.8+
- pandas >= 1.3.0
- numpy >= 1.21.0
- scikit-learn >= 1.0.0
- matplotlib >= 3.3.0 (for plotting)

## Key Concepts

### Heavy/Light Isotope Pairs
The tool expects data with isotope-labeled peptides:
- **Light**: Natural isotope distribution (lower m/z)
- **Heavy**: Isotope-labeled (higher m/z, typically +3-10 Da)
- Ratio = heavy_area / light_area

### Replicate Groups
Measurements are organized into 4 replicate groups based on the last number in the replicate name:
- **Group 1**: 7-14 → Independent regression
- **Group 2**: 17-24 → Independent regression  
- **Group 3**: 27-34 → Independent regression
- **Group 4**: 37-44 → Independent regression

Each group's R² is calculated separately, then aggregated.

### Fragment Charges
Fragment ions with different charges are treated as separate entities:
- `b7_1+` and `b7_2+` get different regressions
- Each charge state has its own R² values

## Configuration (Advanced)

The tool supports custom configurations for different replicate naming schemes:

```bash
# Use custom configuration
python main.py ms_data.csv concentrations.csv -o results.csv -c my_config.yaml
```

Configuration files allow you to customize:
- Replicate name delimiter (e.g., `_`, `-`, `.`)
- Position of replicate number in name
- Replicate group ranges
- Special baseline prefixes

**Example config for sequential numbering (1-4, 5-8, 9-12):**
```yaml
replicate_groups:
  - [1, 4]
  - [5, 8]
  - [9, 12]
```

See `CONFIG_GUIDE.md` for detailed documentation and more examples.

