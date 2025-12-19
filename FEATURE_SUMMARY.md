# Feature: Configurable Replicate Grouping

**Branch:** `feature/configurable-replicate-groups`  
**Status:** Complete  
**Commit:** acd5f91

## Overview

Added YAML-based configuration system that allows users to customize how replicate names are parsed and grouped, without modifying code.

## What Was Added

### 1. Configuration System (`proteomics_cli/config.py`)
- `PRMConfig` class for loading and managing configuration
- Methods for extracting replicate numbers and dilutions
- Configurable group assignment with fallback logic
- Pretty-print configuration summary

### 2. Default Configuration (`prm_config.yaml`)
- Matches original hard-coded behavior
- Delimiter: `_`
- Number position: Last element
- Groups: [7-14], [17-24], [27-34], [37-44]

### 3. Example Configurations (`config_examples/`)
- `sequential_numbering.yaml` - For simple 1-4, 5-8, 9-12 numbering
- `three_groups.yaml` - For experiments with 3 replicate groups

### 4. CLI Integration
- Added `-c/--config` argument to `main.py`
- Optional parameter to `process_prm_data()`
- Prints configuration summary when custom config used

### 5. Documentation
- `CONFIG_GUIDE.md` - Comprehensive guide with examples
- Troubleshooting section
- Best practices

## Usage

### Default Behavior (No Config)
```bash
python main.py ms_data.csv concentrations.csv -o results.csv
```
Works exactly as before - no changes to existing behavior.

### With Custom Config
```bash
python main.py ms_data.csv concentrations.csv -o results.csv \
  -c config_examples/sequential_numbering.yaml
```

### Creating Custom Config
```yaml
# my_experiment.yaml
replicate_parsing:
  delimiter: "_"
  number_position: -1  # or 0 for first element
  dilution_position: 0

replicate_groups:
  - [1, 5]
  - [6, 10]
  - [11, 15]
```

## Configuration Options

### Replicate Parsing
- `delimiter`: Character separating name parts (default: `_`)
- `number_position`: Where replicate number appears (default: `-1` for last)
- `dilution_position`: Where dilution ID appears (default: `0` for first)

### Replicate Groups
- List of `[start, end]` ranges
- Can have any number of groups (not limited to 4)
- Each group gets independent regression

### Special Cases
- `baseline_prefixes`: List of prefixes for baseline samples
- Example: `"col"` → treated as `D0`

## Backward Compatibility

**100% backward compatible**
- Existing commands work without changes
- Default config matches original behavior
- No breaking changes

## Testing Scenarios

### Scenario 1: Sequential Numbering
**Replicate names:** `D1_1`, `D1_2`, ..., `D1_12`  
**Config:** `config_examples/sequential_numbering.yaml`

### Scenario 2: Three Groups
**Replicate names:** `D1_test_5`, `D1_test_15`, `D1_test_25`  
**Config:** `config_examples/three_groups.yaml`

### Scenario 3: Different Delimiter
**Replicate names:** `D1-sample-24`  
**Custom config with:** `delimiter: "-"`

### Scenario 4: Number at Start
**Replicate names:** `001_D1_sample`  
**Custom config with:** `number_position: 0`, `dilution_position: 1`

## Benefits

1. **Flexibility**: Works with different experimental designs
2. **No Code Changes**: Configure via YAML files
3. **Version Control**: Track configs alongside data
4. **Reusability**: Share configs across experiments
5. **Validation**: Prints summary for verification

## Files Modified

- `main.py` - Added `-c/--config` argument
- `proteomics_cli/process_prm_data.py` - Uses config for parsing
- `requirements.txt` - Added PyYAML

## Files Added

- `proteomics_cli/config.py` - Configuration management
- `prm_config.yaml` - Default configuration
- `config_examples/sequential_numbering.yaml` - Example 1
- `config_examples/three_groups.yaml` - Example 2
- `CONFIG_GUIDE.md` - User documentation
- `FEATURE_SUMMARY.md` - This file

## Next Steps

To merge this feature into master:

```bash
# Review changes
git log master..feature/configurable-replicate-groups

# Switch to master
git checkout master

# Merge feature branch
git merge feature/configurable-replicate-groups

# Push to remote
git push origin master
```

## Example Output

When using a custom config:

```
======================================================================
CONFIGURATION SUMMARY
======================================================================
Replicate delimiter: '_'
Number position: -1
Dilution position: 0

Replicate groups (3):
  Group 1: 1-4
  Group 2: 5-8
  Group 3: 9-12

Baseline prefixes: col, control
======================================================================

======================================================================
STEP 1: Loading Data
======================================================================
...
```

## Maintenance

To add support for new configuration options:

1. Update `DEFAULT_CONFIG` in `config.py`
2. Add property method to `PRMConfig` class
3. Update `prm_config.yaml` with new option
4. Document in `CONFIG_GUIDE.md`
