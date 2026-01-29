# Plate Viewer
A robust fluorescence plate-reader analysis pipeline for kinetic assays. Includes per-well baseline fitting, baseline subtraction, normalization, replicate aggregation, mean ± SEM statistics, and kinetic modeling with automated parameter extraction. Designed for multi-plate datasets with noisy baselines and variable well performance. Supports different plate sizes but defaults to 384-well plates.
<img width="1792" height="1032" alt="image" src="https://github.com/user-attachments/assets/884ce123-503a-4c06-bfdc-d758f3550b1f" />
## Overview

Plate Viewer processes Excel files containing fluorescence time-course data from multi-well plates and generates:
- **CSV files**: Individual well data and condition-specific aggregated data (mean ± SEM)
- **Interactive web viewer**: A modern browser-based interface for exploring and analyzing plate data

The tool supports baseline subtraction (no division), baseline normalization, triplicate grouping, regression curve fitting, and statistical analysis. You can copy or download the chart as an image, combine datasets by scientist, and filter by column per scientist.

## Features

### Data Processing Workflow

0. **Baseline Subtraction** (Step 0, optional)
   - Subtract all points from a reference: **lowest point in baseline range** or **first point after baseline**
   - Same scale as raw fluorescence (no division). Step 2 (grouping) can be enabled after Step 0 only, or after Step 1, or both.

1. **Normalization** (Step 1, optional)
   - Fit baseline functions (LOWESS, constant, polynomial) and normalize (e.g. ΔF/F, multiplicative, lowest point, first point, or **none**)
   - Step 2 (grouping) requires **Step 0 or Step 1** (or both) to be enabled.

2. **Triplicate Grouping** (Step 2)
   - Automatically group wells into triplicate groups
   - Calculate mean and SEM across replicates
   - **Requires Step 0 (baseline subtraction) or Step 1 (normalization)** to be enabled

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

### Normalization Modes (Step 1)

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

- **LOWEST POINT IN RANGE** — Normalize by the minimum value in the baseline window: F_norm(t) = F(t) / min(F_baseline).

- **FIRST POINT AFTER BASELINE** — Normalize by the first data point after the baseline window ends (per well).

- **NONE** — No normalization; raw fluorescence values are used. Step 2 (triplicate grouping) still requires Step 0 or Step 1 to be enabled (e.g. enable Step 0 only for baseline-subtracted raw data with grouping).

When **Step 0** (baseline subtraction) and **Step 1** (normalization) are both enabled, the combined pipeline applies (F − ref) / ref (ΔF/F), where ref is from Step 0 (lowest point in baseline range or first point after baseline).

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

- **Select Dataset Groups**: Choose one or more column groups (protein_genotype_buffer combination).
  - **Click**: Select a single group (clears previous selection).
  - **Shift+Click**: Select a range of groups (from last clicked to current).
  - **Ctrl/Cmd+Click**: Toggle a group on or off without clearing others.
- **Column Legend**: Per-scientist column filter — enable/disable specific columns for each scientist’s datasets. Checkboxes apply per scientist so you can show only certain columns per plate group.
- **Well Selector**: Expand “Select Wells” and click wells on the mini grid to **exclude** them from analysis (excluded wells show ✕). Click again to include. Applies across the selected dataset groups.

### Combining Scientist Well Plates

When multiple plates share the same column group but have different scientist initials (from the filename’s 4th part), the viewer:

- **Groups plates by scientist**: Plates are grouped by scientist initials; each scientist gets a separate plate grid section.
- **Horizontal scrolling**: Grids are laid out side-by-side with scroll-snap; use the horizontal scroll or the **plate indicators** (dots below the grids) to jump to a scientist’s section.
- **Plate indicators**: Dots below the plate grids show which scientist section is in view; click a dot to scroll to that scientist’s grid. Tooltips show scientist name (or “No Scientist”).
- **Per-scientist column filter**: In the Column Legend, columns are listed with checkboxes per dataset; enabling/disabling is tracked per scientist so you can filter columns independently for each scientist’s plates.

### Data Processing Options (Workflow)

#### Step 0: Baseline subtraction (optional)

- **Baseline window end (s)**: End time for the baseline range (default: 24).
- **Reference**: **Lowest point in baseline range** or **First point after baseline**.
- Effect: Subtract every point from the chosen reference (same units as raw fluorescence; no division). Step 2 (grouping) can be enabled with Step 0 only, with Step 1 only, or with both. When Step 1 is also enabled, Step 1’s baseline window is used and the pipeline becomes (F − ref) / ref (ΔF/F).

#### Step 1: Normalize using fitted baseline (optional)

- **Baseline window end (s)**: Maximum time for baseline fitting (default: 24).
- **Baseline Fitting Method** and **Normalization mode** are set in the expandable **Baseline Fitting & Normalization Options** section (click to expand). Methods: LOWESS, CONSTANT, or POLYNOMIAL. Modes: ΔF/F, MULTIPLICATIVE, Lowest Point in Range, First Point After Baseline, or **None** (no normalization; raw values).
- Step 2 requires **Step 0 or Step 1** (or both) to be enabled.

#### Step 2: Group wells into triplicate groups (optional)

- Groups wells that belong to triplicate groups (from configuration).
- Calculates mean ± SEM across replicates.
- **Requires Step 0 (baseline subtraction) or Step 1 (normalization)** to be enabled.

#### Step 3: Cutoff data points (optional)

- **Cutoff time (s)**: Maximum time to include (default: 360).
- Excludes points before the first real timepoint and after the cutoff time.

#### Step 4: Fit regression curve (optional)

- **Regression type**: Linear, Polynomial, Exponential, Logarithmic, or Mono-exponential (rise to plateau).
- **Polynomial order**: For polynomial regression (2–5, default: 2).
- Mono-exponential fits show an **info panel** with A, k, τ, t₀, t₁/₂, plateau. Control wells use linear regression by default.

### Chart Features

- **Copy as image**: Copies the current chart as a PNG to the clipboard (falls back to download if copy is not supported).
- **Download as image**: Saves the current chart as a PNG file (`plate-viewer-chart-YYYY-MM-DD-HHMMSS.png`).
- **Interactive legend**: Click to show/hide datasets.
- **Error bars**: SEM for triplicate (or n-plicate) groups.
- **Tooltips**: Hover over points to see values (and ± error when available).
- **Baseline labels**: Baseline values shown on the chart when applicable.

### Well Information and Interactions

- **Plate grid**: Click wells to add/remove them from the chart. Triplicate groups: one click selects/deselects the whole group for that column. Control wells are selected independently.
- **Duplicate wells**: If the same well (same scientist + well ID) appears in multiple files, the well is marked and clicking it selects all matching datasets for comparison.
- **Well info area**: Shows number of selected wells and datasets; prompts to click wells or “Update Chart” as needed.
- **Well selector**: Expand “Select Wells” to exclude specific wells from analysis via the mini grid (excluded = ✕).

## Compiling to AppleScript

The `run_plate_viewer.applescript` file embeds the Python code as base64-encoded text. You can compile it into a quick action which appears in Finder when you right-click a folder.
Please use Automator to create this.

### What the AppleScript Does

The AppleScript:
1. Prompts user to select a folder containing `.xlsx` files (or uses folder from Automator/Finder)
2. Optionally shows an introductory dialog explaining baseline and normalization settings (Step 0 / Step 1 / Step 2) for CSV generation and the web viewer
3. Prompts for **normalization / processing mode**: None (raw values), ΔF/F, Multiplicative, Lowest point, First point after baseline, **Baseline subtract: lowest point (Step 0)**, or **Baseline subtract: first point after baseline (Step 0)** — matching the web viewer’s Step 0 and Step 1
4. Prompts for **baseline fitting method** (when applicable): LOWESS, constant, or polynomial, with clear labels
5. Automatically installs required Python packages if needed:
   - Tries `pip install --user` first
   - Falls back to `--break-system-packages` for Python 3.11+
   - Tries without flags as last resort
6. Decodes and executes the embedded Python code
7. Processes all Excel files in the selected folder
8. Generates CSV files and web viewer
9. Shows a completion notification

**Note**: The Python code is embedded as base64 in the AppleScript. To update it, you would need to:
1. Encode `plate_viewer.py` to base64
2. Replace the `pythonCodeBase64` variable in the AppleScript

How do I do that?
`base64 < plate_viewer.py > /tmp/python_base64.txt && awk -v base64_file="/tmp/python_base64.txt" 'NR==84 {getline base64_content < base64_file; close(base64_file); gsub(/"/, "\\\"", base64_content); print "\t\tset pythonCodeBase64 to \"" base64_content "\""; next} {print}' run_plate_viewer.applescript > /tmp/updated_applescript.applescript && mv /tmp/updated_applescript.applescript run_plate_viewer.applescript`

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

When run interactively, `plate_viewer.py` may prompt for:

- **Normalization / processing mode** (for CSV generation): `none` (raw values), `delta_f_over_f`, `multiplicative`, `lowest_point`, `first_point`, `baseline_subtract_lowest` (Step 0: subtract from lowest point in baseline range), or `baseline_subtract_first` (Step 0: subtract from first point after baseline). These align with the web viewer’s Step 0 and Step 1; baseline-subtract modes produce baseline-subtracted CSVs without division; `none` skips normalization so CSVs contain raw fluorescence.
- **Baseline fitting method** (when a normalization mode that uses a fitted baseline is chosen): LOWESS, constant, or polynomial, with optional parameters (e.g. smoothing fraction, polynomial order).
- **Baseline window end time**: End time (seconds) for the baseline window used in normalization or Step 0.

The script processes all `.xlsx` files in the current directory. Optional command-line arguments (e.g. `--normalization`, `--baseline-method`) may be available; run `python3 plate_viewer.py --help` to check.

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
