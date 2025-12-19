# Configuration Guide

This guide explains how to use configuration files to customize replicate grouping and parsing behavior.

## Quick Start

### Using Default Configuration

By default, the tool uses these settings:
- Replicate delimiter: `_`
- Number position: Last element
- Dilution position: First element
- Groups: [7-14], [17-24], [27-34], [37-44]

```bash
# No config file needed - uses defaults
python main.py ms_data.csv concentrations.csv -o results.csv
```

### Using Custom Configuration

```bash
# Specify custom config file
python main.py ms_data.csv concentrations.csv -o results.csv -c my_config.yaml
```

## Configuration File Format

Configuration files use YAML format. Here's the structure:

```yaml
replicate_parsing:
  delimiter: "_"           # Character that separates replicate name parts
  number_position: -1      # Position of replicate number (0=first, -1=last)
  dilution_position: 0     # Position of dilution identifier

replicate_groups:
  - [7, 14]    # Group 1 range
  - [17, 24]   # Group 2 range
  - [27, 34]   # Group 3 range
  - [37, 44]   # Group 4 range

special_cases:
  baseline_prefixes:
    - "col"      # Names starting with "col" -> D0
    - "control"  # Names starting with "control" -> D0
```

## How Parsing Works

### Replicate Name Structure

Given a replicate name like: `D1_2ul_24`

1. **Split by delimiter**: `["D1", "2ul", "24"]`
2. **Extract number**: Position -1 → `24`
3. **Extract dilution**: Position 0 → `D1`
4. **Assign group**: Number 24 falls in [17-24] → Group 2

### Position Indexing

- `0`: First element
- `1`: Second element
- `-1`: Last element
- `-2`: Second-to-last element

### Examples

**Standard format:** `D1_2ul_24`
```yaml
number_position: -1    # Extracts "24"
dilution_position: 0   # Extracts "D1"
```

**Alternative format:** `001_D1_sample`
```yaml
number_position: 0     # Extracts "001"
dilution_position: 1   # Extracts "D1"
```

**Simple format:** `D1-24` (with hyphen delimiter)
```yaml
delimiter: "-"
number_position: -1    # Extracts "24"
dilution_position: 0   # Extracts "D1"
```

## Example Configurations

### Example 1: Sequential Numbering (1-4, 5-8, 9-12)

For experiments with simple sequential replicate numbers:

```yaml
# config_examples/sequential_numbering.yaml
replicate_parsing:
  delimiter: "_"
  number_position: -1
  dilution_position: 0

replicate_groups:
  - [1, 4]
  - [5, 8]
  - [9, 12]
```

**Usage:**
```bash
python main.py ms_data.csv concentrations.csv -o results.csv \
  -c config_examples/sequential_numbering.yaml
```

### Example 2: Three Replicate Groups

For experiments with only 3 technical replicates:

```yaml
# config_examples/three_groups.yaml
replicate_groups:
  - [1, 10]
  - [11, 20]
  - [21, 30]
```

### Example 3: Different Delimiter

For replicate names like `D1-sample-24`:

```yaml
replicate_parsing:
  delimiter: "-"
  number_position: -1
  dilution_position: 0

replicate_groups:
  - [7, 14]
  - [17, 24]
  - [27, 34]
  - [37, 44]
```

### Example 4: Number at Beginning

For replicate names like `001_D1_test`:

```yaml
replicate_parsing:
  delimiter: "_"
  number_position: 0     # First position
  dilution_position: 1   # Second position

replicate_groups:
  - [1, 10]
  - [11, 20]
  - [21, 30]
  - [31, 40]
```

## Validation

When you run the tool with a custom config, it will print a summary:

```
======================================================================
CONFIGURATION SUMMARY
======================================================================
Replicate delimiter: '_'
Number position: -1
Dilution position: 0

Replicate groups (4):
  Group 1: 7-14
  Group 2: 17-24
  Group 3: 27-34
  Group 4: 37-44

Baseline prefixes: col, control
======================================================================
```

Check this summary to ensure your configuration is correct before processing.

## Troubleshooting

### Problem: "Replicate number is 0"

**Cause:** The tool couldn't extract a number from your replicate names.

**Solution:** Check your `number_position` setting:
- If numbers are at the end: use `-1`
- If numbers are at the beginning: use `0`
- Make sure the `delimiter` is correct

### Problem: "No replicate groups found"

**Cause:** Replicate numbers don't fall within any configured range.

**Solution:** 
1. Check actual replicate numbers in your data
2. Adjust `replicate_groups` ranges to include those numbers

### Problem: "Dilution not found in concentration file"

**Cause:** Extracted dilution doesn't match concentration file columns.

**Solution:**
- Check `dilution_position` setting
- Verify your concentration file has columns like `D0 (ng/mL)`, `D1 (ng/mL)`, etc.

## Best Practices

1. **Start with defaults**: Try running without a config file first
2. **Inspect your data**: Look at actual replicate names before writing config
3. **Test with small subset**: Validate config on a few samples first
4. **Keep configs in version control**: Track different experimental designs
5. **Name configs descriptively**: `seq_1to12.yaml`, `three_groups.yaml`, etc.

## Default Configuration Location

The default configuration is at: `prm_config.yaml`

You can edit this file directly, or create project-specific configs in `config_examples/`.

## Advanced: Fallback Behavior

If a replicate number doesn't match any group range, the tool attempts to infer the group based on the pattern of your configured ranges. This fallback may work for minor numbering variations but is not guaranteed.

For best results, explicitly define all replicate number ranges.
