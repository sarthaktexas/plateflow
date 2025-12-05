# plateflow
A robust fluorescence plate-reader analysis pipeline for kinetic assays. Includes per-well baseline fitting, normalization, replicate aggregation, mean ± SEM statistics, and kinetic modeling with automated parameter extraction. Designed for multi-plate datasets with noisy baselines and variable well performance. Supports different plate sizes but defaults to 384-well plates.
<img width="1792" height="1032" alt="image" src="https://github.com/user-attachments/assets/884ce123-503a-4c06-bfdc-d758f3550b1f" />
## Overview

Plate Viewer processes Excel files containing fluorescence time-course data from multi-well plates and generates:
- **CSV files**: Individual well data and condition-specific aggregated data (mean ± SEM)
- **Interactive web viewer**: A modern browser-based interface for exploring and analyzing plate data

The tool supports advanced data processing workflows including baseline normalization, triplicate grouping, regression curve fitting, and statistical analysis.

## Features

### Data Processing Workflow

1. **Baseline Normalization** (Step 1)
   - Fit baseline functions using multiple methods
   - Normalize data using fitted baselines
   - Filter data points before first timepoint and after cutoff time

2. **Triplicate Grouping** (Step 2)
   - Automatically group wells into triplicate groups
   - Calculate mean and SEM across replicates
   - Requires Step 1 (normalization) to be enabled

3. **Data Cutoff** (Step 3)
   - Exclude data points before the first real timepoint
   - Exclude data points after a specified cutoff time
   - Helps focus analysis on relevant time windows

4. **Regression Curve Fitting** (Step 4)
   - Fit various regression models to data points
   - Supports multiple curve types including mono-exponential (rise to plateau)
   - Display fitted curves overlaid on data

### Baseline Fitting Methods

- **LOWESS** (Locally Weighted Scatterplot Smoothing)
  - Best for: Non-linear baselines with smooth trends
  - Parameters: Smoothing fraction (frac) - higher = smoother, more global fit
  - Formula: g(t) = LOWESS(F_baseline, t_baseline, frac)

- **CONSTANT** (0th-order polynomial)
  - Best for: Stable/flat baselines with minimal drift
  - Parameters: None
  - Formula: g(t) = mean(F_baseline) for all t

- **POLYNOMIAL** (1st-order or higher)
  - Best for: Baselines with linear or polynomial drift
  - Parameters: Polynomial order (1 = linear, 2 = quadratic, etc.)
  - Formula: g(t) = c₀ + c₁t + c₂t² + ...

### Normalization Modes

- **ΔF/F** (Delta F over F)
  - Formula: F'(t) = (F(t) - g(t)) / g(t)
  - Best for: Detecting relative changes from baseline
  - Values can be negative (below baseline) or positive (above baseline)
  - Interpretation: 0 = no change, 0.5 = 50% increase, -0.3 = 30% decrease

- **MULTIPLICATIVE**
  - Formula: F_norm(t) = F(t) / g(t)
  - Best for: Fold-change analysis
  - Values are always positive (ratio to baseline)
  - Interpretation: 1.0 = no change, 2.0 = 2-fold increase, 0.5 = 2-fold decrease

### Regression Curve Types

- **Linear**: y = a + b·x
- **Polynomial**: y = c₀ + c₁x + c₂x² + ... (order 2-5)
- **Exponential**: y = a · e^(b·x)
- **Logarithmic**: y = a + b·ln(x)
- **Mono-exponential (rise to plateau)**: y(t) = 1 + A(1 - e^(-k(t-t₀)))
  - Parameters: A (amplitude), k (rate constant), t₀ (start time), τ (time constant), t_half (half-time), plateau

## Installation

### Prerequisites

- Python 3.x
- Required Python packages:
  - `pandas`
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `openpyxl`
  - `statsmodels` (for LOWESS fitting)

Install dependencies:
```bash
pip install pandas numpy scipy matplotlib openpyxl statsmodels
```

## Usage

### Basic Usage

1. Place your `.xlsx` files in the same directory as `plate_viewer.py`
2. Run the script:
   ```bash
   python3 plate_viewer.py
   ```
3. The script will:
   - Validate filename formats
   - Process all Excel files
   - Generate CSV files in the `csv/` directory
   - Create web viewer files in the `web/` directory
4. Launch the web viewer:
   - Double-click `web.command` in Finder, or
   - Run `./web.command` in Terminal

### Filename Format

Excel files must follow this naming convention:
```
protein_genotype_buffer[_scientist][_number].xlsx
```

Examples:
- `PLCb_WT_Ca_MC.xlsx` (4 parts)
- `PLCb_WT_Ca_MC_13.xlsx` (5 parts)
- `PLCb_H332A_Ca_SM_12.xlsx` (5 parts)

The script extracts:
- **Column group**: First 3 parts (protein_genotype_buffer) - used for grouping datasets
- **Protein**: First part
- **Genotype**: Second part
- **Buffer**: Third part
- **Scientist**: Fourth part (optional)
- **Number**: Fifth part (optional)

### Excel File Format

Excel files should contain:
- A sheet named "Table All Data points" (or the first sheet will be used)
- A row with "Well" in column A and "Time [s]" in column B
- Data rows with:
  - Column A: Well ID (e.g., "B13", "C14")
  - Column B: Content/label
  - Column C+: Time values (in seconds)
  - Subsequent columns: Fluorescence values for each timepoint

## Configuration

### Plate Configuration File

Create or edit `plate_config.json` to customize:

```json
{
  "well_labels": {
    "B": "250-10-2",
    "C": "250-10-2",
    "D": "250-10-2",
    "E": "500-10-2",
    ...
  },
  "triplicate_groups": [
    {
      "name": "250-10-2",
      "wells": ["B", "C", "D"]
    },
    {
      "name": "500-10-2",
      "wells": ["E", "F", "G"]
    },
    ...
  ],
  "control_rows": ["N", "O"]
}
```

**Configuration Options:**

- **well_labels**: Map row letters to condition labels (applies to all columns)
- **triplicate_groups**: Define groups of wells that are replicates
  - `name`: Condition name for the group
  - `wells`: Array of row letters (e.g., ["B", "C", "D"])
- **control_rows**: Array of row letters that are control wells (independently selectable, not grouped)

If no configuration file exists, the script will auto-generate one with defaults based on the data.

## Web Viewer Options

### Dataset Selection

- **Select Dataset Group**: Choose a column group (protein_genotype_buffer combination)
- **Column Legend**: Enable/disable specific columns to filter data
- **Well Selector**: Click wells on the plate grid to include/exclude them from analysis

### Data Processing Options

#### Step 1: Normalize using fitted baseline

**Baseline Fitting Method:**
- **LOWESS**: For non-linear baselines with smooth trends
  - Smoothing fraction (frac): 0.1-1.0 (default: 0.5)
    - Higher (0.7-0.9) = smoother, more global fit
    - Lower (0.2-0.4) = more local, follows data closely
- **CONSTANT**: For stable/flat baselines
- **POLYNOMIAL**: For linear or polynomial drift
  - Polynomial order: 1-5 (default: 1 = linear)

**Baseline Window End Time**: Maximum time (seconds) for baseline window (default: 24.0)

**Normalization Mode:**
- **ΔF/F**: Relative change from baseline (can be negative)
- **MULTIPLICATIVE**: Fold-change (always positive)

#### Step 2: Group wells into triplicate groups

- Groups wells that are part of triplicate groups (from configuration)
- Calculates mean ± SEM across replicates
- **Requires Step 1 to be enabled**

#### Step 3: Cutoff data points

- **Cutoff Time**: Maximum time (seconds) to include (default: 360)
- Excludes data before first real timepoint and after cutoff time

#### Step 4: Fit regression curve to data points

**Regression Type:**
- **Linear**: y = a + b·x
- **Polynomial**: y = c₀ + c₁x + c₂x² + ... (order 2-5)
- **Exponential**: y = a · e^(b·x)
- **Logarithmic**: y = a + b·ln(x)
- **Mono-exponential (rise to plateau)**: y(t) = 1 + A(1 - e^(-k(t-t₀)))
  - Displays parameters: A (amplitude), k (rate constant), τ (time constant), t₀ (start time), t_half (half-time), plateau

**Polynomial Order**: For polynomial regression (2-5, default: 2)

### Chart Features

- **Interactive legend**: Click to show/hide datasets
- **Error bars**: Display SEM for triplicate groups
- **Tooltips**: Hover over data points to see values
- **Zoom/Pan**: Use mouse wheel and drag to navigate
- **Baseline labels**: Display baseline values on chart (when applicable)

### Well Information

- Click wells on the plate grid to see:
  - Well ID
  - Content/label
  - Row-based label (from configuration)
  - Dataset count (if multiple datasets have data for this well)

## Compiling to AppleScript

The `run_plate_viewer.applescript` file embeds the Python code as base64-encoded text. You can compile it into a quick action which appears in Finder when you right-click a folder.
Please use Automator to create this.

### What the AppleScript Does

The AppleScript:
1. Prompts user to select a folder containing `.xlsx` files (or uses folder from Automator/Finder)
2. Automatically installs required Python packages if needed:
   - Tries `pip install --user` first
   - Falls back to `--break-system-packages` for Python 3.11+
   - Tries without flags as last resort
3. Decodes and executes the embedded Python code
4. Processes all Excel files in the selected folder
5. Generates CSV files and web viewer
6. Shows a completion notification

**Note**: The Python code is embedded as base64 in the AppleScript. To update it, you would need to:
1. Encode `plate_viewer.py` to base64
2. Replace the `pythonCodeBase64` variable in the AppleScript

How do I do that?
`base64 < plate_viewer.py > /tmp/python_base64.txt && awk -v base64_file="/tmp/python_base64.txt" 'NR==84 {getline base64_content < base64_file; close(base64_file); gsub(/"/, "\\\"", base64_content); print "\t\tset pythonCodeBase64 to \"" base64_content "\""; next} {print}' run_plate_viewer.applescript > /tmp/updated_applescript.applescript && mv /tmp/updated_applescript.applescript run_`

## Output Files

### CSV Files

Generated in `csv/<plate_id>/`:

1. **Per-well CSVs**: `{plate_id}_{well_id}.csv`
   - Columns: `time_s`, `value`
   - One file per well

2. **Condition-specific CSVs**: `{plate_id}_{condition_name}_col{column}.csv`
   - Columns: `time_s`, `mean`, `sem`
   - One file per condition-column combination
   - Contains mean ± SEM across triplicate wells

### Web Viewer Files

Generated in `web/`:

- `index.html`: Main web viewer interface
- `viewer_data.json`: All plate data in JSON format
- `web.command`: Double-clickable launcher for local web server

## Statistics Computation

For each condition and timepoint:

- **Mean**: μ(t) = (1/n) · Σ F'(t)
  - Average normalized value across all wells in the condition

- **SD**: σ(t) = sqrt(Σ(F'(t) - μ(t))² / (n-1))
  - Sample standard deviation (dof=1) across wells

- **SEM**: SEM(t) = σ(t) / sqrt(n)
  - Standard error of the mean
  - Note: SEM is computed across wells, not propagated from baseline error

## Troubleshooting

### "No .xlsx files found"
- Ensure Excel files are in the same directory as `plate_viewer.py`
- Check that files have `.xlsx` extension (not `.xls`)

### "Invalid filename format"
- Rename files to match: `protein_genotype_buffer[_scientist][_number].xlsx`
- Ensure at least 3 underscore-separated parts before `.xlsx`

### "Duplicate datasets found"
- Multiple files are producing the same plate_id
- Delete duplicate files and rerun

### "Could not find 'Well' header row"
- Ensure Excel file has "Well" in column A of the header row
- Check that "Time [s]" appears in column B of the same or next row

### Web server won't start
- Ensure Python 3 is installed: `python3 --version`
- Try running manually: `cd web && python3 -m http.server 8000`
- Check if port 8000 is already in use

### Baseline fitting fails
- Ensure `statsmodels` is installed for LOWESS: `pip install statsmodels`
- Check that baseline window contains sufficient data points
- Try a different baseline method (CONSTANT or POLYNOMIAL)

## Advanced Usage

### Command Line Options

Currently, the script processes all `.xlsx` files in the current directory. No command-line arguments are supported yet.

### Customizing Defaults

Edit constants in `plate_viewer.py`:

- `PLATE_ROWS`: Row letters (default: "ABCDEFGHIJKLMNOP")
- `PLATE_COLS`: Number of columns (default: 24)
- `SHEET_NAME_PREFERRED`: Preferred Excel sheet name (default: "Table All Data points")
- `DEFAULT_WELL_LABELS`: Default well label mappings
- `DEFAULT_TRIPLICATE_GROUPS`: Default triplicate group definitions
- `TRIPLICATE_COLORS`: Color palette for triplicate groups

## License

This is not licensed. If you want to use, you may but please credit me unless you belong to the [Falzone Lab](https://falzonelab.com).

Made with thought 
