# Plate Viewer
A fluorescence plate-reader analysis tool for kinetic assays. The Python script parses Excel plate-reader exports and builds an interactive web viewer. Baseline fitting, baseline subtraction, normalization, replicate aggregation, mean ± SEM statistics, regression curve fitting, and kinetic modeling all run in the browser. Designed for multi-plate datasets with noisy baselines and variable well performance. Supports different plate sizes but defaults to 384-well plates.
<img width="1792" height="1032" alt="image" src="https://github.com/user-attachments/assets/884ce123-503a-4c06-bfdc-d758f3550b1f" />
## Overview

Plate Viewer processes Excel files containing fluorescence time-course data from multi-well plates and generates an **interactive web viewer** for exploring and analyzing plate data.

**What the script does:** validates filenames, parses Excel files, resolves multi-plate layout conflicts, loads or auto-generates `plate_config.json`, and writes `web/viewer_data.json` plus `web/index.html`.

**What the web viewer does:** baseline subtraction, normalization, triplicate grouping, statistics, regression fitting, well exclusion, and chart export (PNG or CSV). There is no separate CSV export step or `csv/` output folder.

## Features

### Configuration and Run Logs

- **Auto plate config**: On first run (no `plate_config.json`), well labels are detected from the Excel **Content** column and merged across all files.
- **Existing config preserved**: `plate_config.json` is never overwritten with defaults or re-detected from Excel once it exists.
- **Forgiving JSON**: Hand-edited config files tolerate common mistakes (e.g. trailing commas after the last entry). Minor syntax issues are auto-repaired and saved as valid JSON on load.
- **Run logs**: Warnings and errors are written to `WARNINGS.txt` and `ERRORS.txt` in the project folder.

### Data Processing Workflow

All steps below are configured in the **web viewer** (not at script run time):

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
2. Run the script (non-interactive — no baseline or normalization prompts):
   ```bash
   python3 plate_viewer.py
   ```
3. The script will:
   - Validate filename formats
   - Process all Excel files
   - Create web viewer files in the `web/` directory
   - Write `WARNINGS.txt` and/or `ERRORS.txt` when issues are detected (see [Run Logs](#run-logs))
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
- A header row with **Well** in column A and **Content** in column B
- An optional **Group** column in column C (BMG Voyager exports; ignored by the parser)
- A **Time [s]** row directly below the header row (typically in column B)
- Numeric time values (seconds) starting in the first column after any non-numeric cells on the Time row (column D when a Group column is present)
- Fluorescence channel headers (e.g. `Raw Data (485-20/535-20)`) may appear on the header row above the time values

**Layout example (with optional Group column):**

| | A | B | C | D | E | … |
|---|---|---|---|---|---|---|
| Header | Well | Content | Group | Raw Data … | Raw Data … | … |
| Time | | Time [s] | | 0 | 20.5 | … |
| Data | B02 | SEC only | A | 1191 | 989421 | … |

**Data columns:**
- Column A: Well ID (e.g., "B13", "C14")
- Column B: Content/label (used for auto-detecting `plate_config.json` on first run)
- Column C: Group (optional; not used by the viewer)
- Column D+: Fluorescence values for each timepoint (aligned with the time row)

## Configuration

### Plate Configuration File

Create or edit `plate_config.json` to customize well labels and control rows:

```json
{
  "well_labels": {
    "B": "SOS1 - 2uM",
    "C": "SOS1 - 2uM",
    "D": "SOS1 - 2uM",
    "E": "SOS1 - 4uM",
    "F": "SOS1 - 4uM",
    "G": "SOS1 - 4uM"
  },
  "control_rows": ["N", "O"]
}
```

**Configuration Options:**

- **well_labels**: Map row letters to condition labels (applies to all columns). Triplicate groups in the web viewer are **auto-generated** from rows that share the same label.
- **control_rows**: Row letters for control wells (independently selectable, not grouped). When using a plate config, no control rows are assumed unless you specify them.

To exclude individual wells from analysis, use the **Well Selector** in the web viewer (not `plate_config.json`). Per-well exclusion works per dataset in the viewer and is the supported approach when the same well position appears on multiple plates.

#### Editing the config safely

You can edit `plate_config.json` in any text editor. The loader is lenient about small syntax mistakes that strict JSON does not allow:

- **Trailing commas** after the last item in an object or array (e.g. `"M": "label",` before `}`) are accepted and automatically fixed.
- On repair, the script rewrites the file as valid JSON and prints: `Repaired minor JSON syntax in plate_config.json (e.g. trailing commas)`.
- Your labels and settings are preserved — repairs only fix formatting, not content.

If the file has a **serious JSON error** that cannot be auto-repaired, the script:

1. Leaves your file **unchanged** on disk
2. Records the error in `WARNINGS.txt`
3. Uses empty defaults for **that run only** (no well labels or control rows)

Fix the JSON and rerun; nothing is reset to Excel-detected or built-in defaults.

#### Auto-detection from Excel (first run only)

If **`plate_config.json` does not exist**, the script creates it automatically by reading the **Content** column (column B) in your Excel files:

1. For each well, the row letter (e.g. `B` from `B13`) is mapped to that well’s content label.
2. **Multiple Excel files are merged**: if one file has content in rows B, C, D and another has B, C, E, the generated config includes all of those rows.
3. **Conflicting labels** (same row, different content in different files, or multiple labels within one file for the same row) are recorded as warnings in `WARNINGS.txt`.

**Important:** An existing `plate_config.json` is **never replaced or regenerated from Excel**. Auto-detection runs only when the file is missing. To regenerate from Excel, delete `plate_config.json` and rerun the script.

## Run Logs

During each run, the script collects warnings and errors and writes them to the project folder:

| File | Contents |
|------|----------|
| `WARNINGS.txt` | Non-fatal issues (label conflicts, duplicate column overlaps, missing preferred Excel sheet, config load issues, etc.) |
| `ERRORS.txt` | Fatal or serious issues (invalid filenames, parse failures, duplicate plate IDs, no data found, etc.) |

- Messages are also printed to the terminal.
- If a run completes with no warnings or errors, the corresponding file is removed (so stale logs from a previous run are not left behind).
- On fatal exit, log files are written before the script stops.

**Examples of what appears in `WARNINGS.txt`:**
- Content label conflicts when auto-generating `plate_config.json`
- Unrecoverable `plate_config.json` parse errors (file left unchanged; fix JSON and rerun)
- Duplicate wells on the same layout across files (only the first file alphabetically is used in the viewer)
- Excel sheet `"Table All Data points"` not found (first sheet used instead)

**Examples of what appears in `ERRORS.txt`:**
- Invalid filename format
- Failed Excel parse for a specific file
- Duplicate datasets (same plate ID from multiple files)
- No `.xlsx` files or no plates successfully parsed

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
- **Download as CSV**: Exports the currently displayed chart series (and SEM columns when shown) as a CSV file.
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

The AppleScript:
1. Prompts user to select a folder containing `.xlsx` files (or uses folder from Automator/Finder)
2. Automatically installs required Python packages if needed:
   - Tries `pip install --user` first
   - Falls back to `--break-system-packages` for Python 3.11+
   - Tries without flags as last resort
3. Decodes and executes the embedded Python code
4. Processes all Excel files in the selected folder and refreshes the web viewer data
5. Shows a completion notification

**Note**: The Python code is embedded as base64 in the AppleScript. To update it, you would need to:
1. Encode `plate_viewer.py` to base64
2. Replace the `pythonCodeBase64` variable in the AppleScript

How do I do that?
`base64 < plate_viewer.py > /tmp/python_base64.txt && awk -v base64_file="/tmp/python_base64.txt" 'NR==84 {getline base64_content < base64_file; close(base64_file); gsub(/"/, "\\\"", base64_content); print "\t\tset pythonCodeBase64 to \"" base64_content "\""; next} {print}' run_plate_viewer.applescript > /tmp/updated_applescript.applescript && mv /tmp/updated_applescript.applescript run_plate_viewer.applescript`

## Output Files

### Run Logs

Written in the same directory as `plate_viewer.py` (when applicable):

- `WARNINGS.txt` — non-fatal issues from the current run
- `ERRORS.txt` — errors from the current run (including runs that exit early)

See [Run Logs](#run-logs) for details.

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

Check **`WARNINGS.txt`** and **`ERRORS.txt`** in the project folder first — they list issues from the most recent run with full detail.

### "No .xlsx files found"
- Ensure Excel files are in the same directory as `plate_viewer.py`
- Check that files have `.xlsx` extension (not `.xls`)

### "Invalid filename format"
- Rename files to match: `protein_genotype_buffer[_scientist][_number].xlsx`
- Ensure at least 3 underscore-separated parts before `.xlsx`

### "Duplicate datasets found"
- Multiple files are producing the same plate_id
- Delete duplicate files and rerun
- See `ERRORS.txt` for which plate ID and files are involved

### Content label conflicts
- Occurs when auto-generating `plate_config.json` and different Excel files assign different labels to the same row
- Details are in `WARNINGS.txt`; the first-seen label is used in the generated config
- Fix by editing `plate_config.json` manually, or delete it and resolve conflicts in the Excel content column before rerunning

### `plate_config.json` could not be loaded
- Usually caused by invalid JSON that auto-repair cannot fix (unclosed quotes, missing commas between keys, etc.)
- Trailing commas after the last entry are fixed automatically — you should not see a warning for those
- Check `WARNINGS.txt` for the exact parse error and line number
- Your file on disk is **not** modified; fix the JSON and rerun
- If you want to start over, delete `plate_config.json` and rerun to regenerate from Excel

### "Could not find 'Well' header row"
- Ensure Excel file has "Well" in column A of the header row
- Check that "Time [s]" appears on the row directly below the header (usually column B)

### "No timepoints parsed"
- Usually means time values were not found on the Time row — often because an optional **Group** column in column C shifted the data without being accounted for
- Confirm the Time row has numeric seconds (e.g. `0`, `20.5`) starting in column D when column C is **Group**, or column C when there is no Group column
- Ensure the header row uses **Group** (not a different label) in column C if that column is present

### Web server won't start
- Ensure Python 3 is installed: `python3 --version`
- Try running manually: `cd web && python3 -m http.server 8000`
- Check if port 8000 is already in use

### Baseline fitting fails
- Ensure `statsmodels` is installed for LOWESS: `pip install statsmodels`
- Check that baseline window contains sufficient data points
- Try a different baseline method (CONSTANT or POLYNOMIAL)

## Advanced Usage

`plate_viewer.py` is a one-shot, non-interactive pipeline: it processes all `.xlsx` files in the current directory and regenerates the web viewer. Baseline subtraction, normalization, well exclusion, and chart export are configured entirely in the web interface.

There is no `csv/` output directory and no command-line flags for baseline or normalization settings.

### Customizing Defaults

Edit constants in `plate_viewer.py`:

- `PLATE_ROWS`: Row letters (default: "ABCDEFGHIJKLMNOP")
- `PLATE_COLS`: Number of columns (default: 24)
- `SHEET_NAME_PREFERRED`: Preferred Excel sheet name (default: "Table All Data points")
- `TRIPLICATE_COLORS`: Color palette for triplicate groups

Well labels are normally defined in `plate_config.json` (auto-detected from the Excel content column on first run). The script does not use hardcoded default well labels.

## License

This is not licensed. If you want to use, you may but please credit me unless you belong to the [Falzone Lab](https://falzonelab.com).

Made with thought 
