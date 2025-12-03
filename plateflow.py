#!/usr/bin/env python3
import os
import sys
import json
from typing import Dict, Any, List

import pandas as pd
import numpy as np

# ----- CONFIG -----
# 384-well plate: 16 rows (A-P) x 24 columns; change if needed.
PLATE_ROWS = "ABCDEFGHIJKLMNOP"  # 16 rows
PLATE_COLS = 24

SHEET_NAME_PREFERRED = "Table All Data points"  # from your example

# Default well labels and triplicate groups
DEFAULT_WELL_LABELS = {
    "B13": "250-10-2",
    "C13": "250-10-2",
    "D13": "250-10-2",
    "E13": "500-10-2",
    "F13": "500-10-2",
    "G13": "500-10-2",
    "H13": "250-5-2",
    "I13": "250-5-2",
    "J13": "250-5-2",
    "K13": "500-5-2",
    "L13": "500-5-2",
    "M13": "500-5-2",
    "N13": "500-10-2",
    "O13": "500-5-2"
}

DEFAULT_TRIPLICATE_GROUPS = [
    {
        "name": "250-10-2",
        "wells": ["B", "C", "D"]
    },
    {
        "name": "500-10-2",
        "wells": ["E", "F", "G"]
    },
    {
        "name": "250-5-2",
        "wells": ["H", "I", "J"]
    },
    {
        "name": "500-5-2",
        "wells": ["K", "L", "M"]
    }
]

# Color palette for triplicate groups - matches web interface
# Each group gets a distinct color (using border color for plots)
TRIPLICATE_COLORS = [
    "#ff6b6b",  # Red
    "#4ecdc4",  # Teal
    "#45b7d1",  # Blue
    "#f9ca24",  # Yellow
    "#6c5ce7",  # Purple
    "#a29bfe",  # Light Purple
    "#fd79a8",  # Pink
    "#00b894",  # Green
    "#e17055",  # Orange
    "#74b9ff",  # Light Blue
    "#55efc4",  # Mint
    "#fdcb6e"   # Light Orange
]

def is_control_well(well_id: str, control_rows: List[str] = None) -> bool:
    """
    Check if a well is a control well.
    Default control rows are N and O, but this can be configured via control_rows parameter.
    Control wells should be independently selectable and not grouped as triplicates.
    """
    if not well_id:
        return False
    if control_rows is None:
        control_rows = ["N", "O"]  # Default control rows
    row_letter = well_id[0] if well_id[0].isalpha() else ""
    return row_letter in control_rows


def validate_filename_format(filename: str):
    """
    Validate that filename matches the expected format:
    protein_genotype_buffer[_scientist][_number].xlsx
    
    Expected format: At least 3 underscore-separated parts before the .xlsx extension.
    Examples:
    - PLCb_WT_Ca_MC.xlsx (4 parts)
    - PLCb_WT_Ca_MC_13.xlsx (5 parts)
    - PLCb_H332A_Ca_SM_12.xlsx (5 parts)
    
    Raises ValueError if format is invalid.
    """
    if not filename.lower().endswith('.xlsx'):
        raise ValueError(
            f"Invalid filename format: '{filename}'\n"
            f"Expected format: protein_genotype_buffer[_scientist][_number].xlsx\n"
            f"Filename must end with .xlsx extension."
        )
    
    # Remove .xlsx extension
    base_name = os.path.splitext(filename)[0]
    
    # Split by underscore
    parts = base_name.split("_")
    
    if len(parts) < 3:
        raise ValueError(
            f"Invalid filename format: '{filename}'\n"
            f"Expected format: protein_genotype_buffer[_scientist][_number].xlsx\n"
            f"The filename must have at least 3 underscore-separated parts before the .xlsx extension.\n"
            f"Found {len(parts)} part(s): {parts}\n"
            f"Example valid formats:\n"
            f"  - PLCb_WT_Ca_MC.xlsx\n"
            f"  - PLCb_WT_Ca_MC_13.xlsx\n"
            f"  - PLCb_H332A_Ca_SM_12.xlsx"
        )


def parse_column_group_from_filename(filename: str) -> Dict[str, str]:
    """
    Parse column group information from filename format:
    protein_genotype_buffer_initials-of-scientist_arbitrary-number.xlsx
    
    Returns a dictionary with:
    - "column_group": "protein_genotype_buffer" (used for grouping)
    - "protein": protein name
    - "genotype": genotype
    - "buffer": buffer
    - "scientist": initials-of-scientist
    - "number": arbitrary-number
    """
    # Validate filename format first
    validate_filename_format(filename)
    
    # Remove .xlsx extension
    base_name = os.path.splitext(filename)[0]
    
    # Split by underscore
    parts = base_name.split("_")
    
    # At this point we know we have at least 3 parts due to validation
    # If we have at least 3 parts, use first 3 for column_group (protein_genotype_buffer)
    # This handles cases like "PLCb_WT_EGTA_MC" (4 parts) or "PLCb_WT_EGTA_MC_13" (5 parts)
    protein = parts[0]
    genotype = parts[1]
    buffer = parts[2]
    scientist = parts[3] if len(parts) > 3 else ""
    number = parts[4] if len(parts) > 4 else ""
    
    # Column group is always protein_genotype_buffer (first 3 parts)
    column_group = f"{protein}_{genotype}_{buffer}"
    
    return {
        "column_group": column_group,
        "protein": protein,
        "genotype": genotype,
        "buffer": buffer,
        "scientist": scientist,
        "number": number
    }
    
    # If we have at least 3 parts, use first 3 for column_group (protein_genotype_buffer)
    # This handles cases like "PLCb_WT_EGTA_MC" (4 parts) or "PLCb_WT_EGTA_MC_13" (5 parts)
    protein = parts[0]
    genotype = parts[1]
    buffer = parts[2]
    scientist = parts[3] if len(parts) > 3 else ""
    number = parts[4] if len(parts) > 4 else ""
    
    # Column group is always protein_genotype_buffer (first 3 parts)
    column_group = f"{protein}_{genotype}_{buffer}"
    
    return {
        "column_group": column_group,
        "protein": protein,
        "genotype": genotype,
        "buffer": buffer,
        "scientist": scientist,
        "number": number
    }
    


def find_header_and_time_rows(df: pd.DataFrame) -> (int, int):
    """
    Find the 'Well' header row and the 'Time [s]' row.
    Returns (header_row_index, time_row_index).
    """
    header_row_idx = None
    time_row_idx = None

    for i in range(len(df)):
        cell0 = str(df.iat[i, 0]) if not pd.isna(df.iat[i, 0]) else ""
        cell1 = str(df.iat[i, 1]) if not pd.isna(df.iat[i, 1]) else ""
        if cell0.strip() == "Well":
            header_row_idx = i
        if "Time" in cell1 and "[s]" in cell1:
            time_row_idx = i

    if header_row_idx is None or time_row_idx is None:
        raise ValueError("Could not find 'Well' header row or 'Time [s]' row in Excel sheet.")

    if time_row_idx != header_row_idx + 1:
        # Not fatal, but warn in case format drifts
        print(f"Warning: Time row (index {time_row_idx}) is not directly after header row (index {header_row_idx}).",
              file=sys.stderr)

    return header_row_idx, time_row_idx


def parse_plate_file(xlsx_path: str) -> pd.DataFrame:
    """
    Parse a single Excel file into a long-format DataFrame with columns:
    ['plate_id', 'well', 'content', 'time_s', 'value']
    """
    print(f"Parsing {os.path.basename(xlsx_path)}...", flush=True)
    try:
        print(f"  Opening Excel file...", flush=True)
        # Explicitly use openpyxl engine for .xlsx files to avoid hanging
        xls = pd.ExcelFile(xlsx_path, engine='openpyxl')
        if SHEET_NAME_PREFERRED in xls.sheet_names:
            sheet_name = SHEET_NAME_PREFERRED
        else:
            sheet_name = xls.sheet_names[0]
            print(f"  Sheet '{SHEET_NAME_PREFERRED}' not found; using first sheet '{sheet_name}' instead.", flush=True)
        print(f"  Reading sheet '{sheet_name}'...", flush=True)
        # Use the same engine for read_excel
        df = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None, engine='openpyxl')
        print(f"  Sheet loaded: {len(df)} rows, {len(df.columns)} columns", flush=True)
        # Close the ExcelFile to free resources
        xls.close()
    except FileNotFoundError:
        raise RuntimeError(f"File not found: {xlsx_path}")
    except PermissionError:
        raise RuntimeError(f"Permission denied: {xlsx_path} (file may be open in another program)")
    except Exception as e:
        error_type = type(e).__name__
        raise RuntimeError(f"Failed to read {xlsx_path}: {error_type}: {e}")

    print(f"  Finding header and time rows...", flush=True)
    header_row_idx, time_row_idx = find_header_and_time_rows(df)

    # Time row: column 1 is "Time [s]", measurement times start from column 2 onward.
    time_row = df.iloc[time_row_idx]
    time_values = time_row.iloc[2:].tolist()
    # Clean up times
    time_values_clean = []
    for t in time_values:
        if pd.isna(t):
            break
        try:
            time_values_clean.append(float(t))
        except Exception:
            # Stop at the first non-numeric / weird thing
            break

    num_time_points = len(time_values_clean)
    if num_time_points == 0:
        raise ValueError(f"No timepoints parsed in {xlsx_path}.")

    plate_id = os.path.splitext(os.path.basename(xlsx_path))[0]

    records: List[Dict[str, Any]] = []

    # Data rows: from time_row_idx + 1 downwards until we hit empty "Well" cell
    print(f"  Processing data rows (starting from row {time_row_idx + 1})...", flush=True)
    rows_processed = 0
    for i in range(time_row_idx + 1, len(df)):
        well = df.iat[i, 0]
        content = df.iat[i, 1]

        if pd.isna(well) or str(well).strip() == "":
            # Assume we've reached the end of data
            continue

        well_id = str(well).strip()
        content_str = "" if pd.isna(content) else str(content).strip()

        row_vals = df.iloc[i, 2:2 + num_time_points].tolist()
        # Ensure same length as time list
        row_vals = row_vals[:num_time_points]

        for t, v in zip(time_values_clean, row_vals):
            if pd.isna(v):
                # Skip missing values, but keep timepoint in overall scheme
                continue
            try:
                v_float = float(v)
            except Exception:
                continue

            records.append(
                {
                    "plate_id": plate_id,
                    "well": well_id,
                    "content": content_str,
                    "time_s": float(t),
                    "value": v_float,
                }
            )
        
        rows_processed += 1
        if rows_processed % 100 == 0:
            print(f"    Processed {rows_processed} rows...", flush=True)

    print(f"  Processed {rows_processed} rows, created {len(records)} data points", flush=True)

    if not records:
        raise ValueError(f"No data rows found in {xlsx_path} after parsing.")

    print(f"  Creating DataFrame...", flush=True)
    long_df = pd.DataFrame.from_records(records)
    print(f"  ✓ Completed parsing {os.path.basename(xlsx_path)}", flush=True)
    return long_df


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def write_csvs_for_plate(long_df: pd.DataFrame, csv_root: str, config: Dict[str, Any] = None, filename: str = None):
    """
    Write CSV files:
      - csv_root/<plate_id>/<plate_id>_<well>.csv (per well, without SEM)
      - csv_root/<condition>_<triplicate_name>.csv (per condition-triplicate combination, with mean and SEM)
    """
    plate_id = long_df["plate_id"].iloc[0]
    print(f"  Writing CSVs for {plate_id}...", flush=True)

    # Load config if not provided
    if config is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "plate_config.json")
        config = load_config(config_path)
    
    # Parse column info from filename if provided
    column_info = None
    if filename:
        try:
            column_info = parse_column_group_from_filename(filename)
        except Exception:
            pass
    
    # Get triplicate groups from config
    triplicate_groups = config.get("triplicate_groups", DEFAULT_TRIPLICATE_GROUPS)
    control_rows = config.get("control_rows", ["N", "O"])
    
    # Build mapping: well_id -> triplicate group name
    well_to_group = {}
    for group in triplicate_groups:
        group_name = group.get("name", "")
        group_wells = group.get("wells", [])
        for well_pattern in group_wells:
            # Extract row letter from well pattern (e.g., "B" from "B" or "B13")
            row_letter = well_pattern[0] if well_pattern and well_pattern[0].isalpha() else ""
            if row_letter:
                # Find all wells in this row across all columns
                for well_id in long_df["well"].unique():
                    if well_id and well_id[0] == row_letter:
                        well_to_group[well_id] = group_name
    
    # Per-well CSVs (without SEM - SEM only belongs in triplicate CSVs)
    plate_folder = os.path.join(csv_root, plate_id)
    ensure_dir(plate_folder)
    
    wells = long_df.groupby("well")
    print(f"    Writing {len(wells)} per-well CSVs...", flush=True)
    
    # Write per-well CSVs
    for idx, (well_id, g) in enumerate(wells, 1):
        safe_well = well_id.replace(" ", "").replace(":", "-")
        out_path = os.path.join(plate_folder, f"{plate_id}_{safe_well}.csv")
        g_sorted = g.sort_values("time_s")
        
        # Create output DataFrame with time and value only (no SEM)
        output_df = pd.DataFrame({
            "time_s": g_sorted["time_s"].values,
            "value": g_sorted["value"].values
        })
        output_df.to_csv(out_path, index=False)
        
        if idx % 5 == 0:
            print(f"      Written {idx}/{len(wells)} wells...", flush=True)
    
    # Write condition-specific CSVs (genotype-buffer-triplicate)
    # One CSV per condition-triplicate-column combination
    if column_info:
        print(f"    Writing condition-specific CSVs...", flush=True)
        genotype = column_info.get("genotype", "")
        buffer = column_info.get("buffer", "")
        
        # Group by triplicate group name and column
        for group in triplicate_groups:
            group_name = group.get("name", "")
            if not group_name:
                continue
            
            # Get all wells in this triplicate group (across all columns)
            group_wells = []
            for well_id in long_df["well"].unique():
                if well_id in well_to_group and well_to_group[well_id] == group_name:
                    group_wells.append(well_id)
            
            if not group_wells:
                continue
            
            # Group wells by column number
            wells_by_column = {}
            for well_id in group_wells:
                # Extract column number from well ID (e.g., "B13" -> 13)
                try:
                    column_num = int(''.join(filter(str.isdigit, str(well_id))))
                    if column_num not in wells_by_column:
                        wells_by_column[column_num] = []
                    wells_by_column[column_num].append(well_id)
                except (ValueError, AttributeError):
                    continue
            
            # Create one CSV per column for this triplicate group
            for column_num in sorted(wells_by_column.keys()):
                column_wells = wells_by_column[column_num]
                
                # Get all time points for this column's wells
                all_times = sorted(long_df[long_df["well"].isin(column_wells)]["time_s"].unique())
                
                # Calculate mean and SEM for each time point across wells in this triplicate group for this column
                mean_values = []
                sem_values = []
                
                for time in all_times:
                    values_at_time = []
                    for well_id in column_wells:
                        well_data = long_df[(long_df["well"] == well_id) & (long_df["time_s"] == time)]
                        if not well_data.empty:
                            values_at_time.append(well_data["value"].iloc[0])
                    
                    if len(values_at_time) >= 2:
                        mean = sum(values_at_time) / len(values_at_time)
                        variance = sum((v - mean) ** 2 for v in values_at_time) / (len(values_at_time) - 1) # Sample variance using Bessel's correction/sample standard deviation
                        std_dev = variance ** 0.5
                        sem = std_dev / (len(values_at_time) ** 0.5)
                        mean_values.append(mean)
                        sem_values.append(sem)
                    elif len(values_at_time) == 1:
                        mean_values.append(values_at_time[0])
                        sem_values.append(None)
                    else:
                        mean_values.append(None)
                        sem_values.append(None)
                
                # Create condition CSV filename: plate_id_triplicate_name_col{column}
                # Plate ID already contains genotype and buffer, so just use triplicate name
                # Keep hyphens in group name (e.g., "250-10-2"), just sanitize spaces and colons
                safe_group_name = group_name.replace(" ", "_").replace(":", "-")
                condition_path = os.path.join(csv_root, f"{plate_id}_{safe_group_name}_col{column_num}.csv")
                
                # Write condition CSV
                condition_df = pd.DataFrame({
                    "time_s": all_times,
                    "mean": mean_values,
                    "sem": sem_values
                })
                condition_df.to_csv(condition_path, index=False)
        
        print(f"      Written condition-specific CSVs", flush=True)
    
    print(f"  ✓ Completed writing CSVs for {plate_id}", flush=True)


def write_default_config(config_path: str):
    """Write default configuration file from defaults."""
    # Convert DEFAULT_WELL_LABELS to row-based format
    well_labels = {}
    for well_id, label in DEFAULT_WELL_LABELS.items():
        if well_id and len(well_id) > 0 and well_id[0].isalpha():
            row_letter = well_id[0]
            if row_letter not in well_labels:
                well_labels[row_letter] = label
    
    default_config = {
        "well_labels": well_labels,
        "triplicate_groups": DEFAULT_TRIPLICATE_GROUPS,
        "control_rows": ["N", "O"]
    }
    
    try:
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(default_config, f, indent=2)
        print(f"  Generated default configuration file: {os.path.basename(config_path)}", flush=True)
    except Exception as e:
        print(f"  Warning: Could not write default config file: {e}", flush=True)


def load_config(config_path: str) -> Dict[str, Any]:
    """Load plate configuration file if it exists, otherwise generate from defaults."""
    if os.path.exists(config_path):
        try:
            with open(config_path, "r", encoding="utf-8") as f:
                config = json.load(f)
                print(f"  Loaded configuration from {os.path.basename(config_path)}", flush=True)
                return config
        except Exception as e:
            print(f"  Warning: Could not load config file: {e}", flush=True)
            # Generate default config if loading fails
            write_default_config(config_path)
            return {}
    else:
        # Generate default config file if it doesn't exist
        write_default_config(config_path)
        return {}


def build_viewer_json(all_plates: Dict[str, pd.DataFrame], json_path: str, config: Dict[str, Any] = None, plate_id_to_filename: Dict[str, str] = None):
    """
    Build viewer_data.json with structure:
    {
      "config": {
        "well_labels": {...},
        "triplicate_groups": [...],  # row-letter groups (applies to all columns)
        "column_groups": {...}  # mapping of column_group to list of plate_ids
      },
      "plates": [
        {
          "id": "<plate_id>",
          "file": "<original_file_name>",
          "column_group": "<protein_genotype_buffer>",
          "column_info": {
            "protein": "...",
            "genotype": "...",
            "buffer": "...",
            "scientist": "...",
            "number": "..."
          },
          "wells": {
            "B15": {
              "content": "Sample X1",
              "label": "Sample X1",  # from config if available
              "time_s": [...],
              "values": [...]
            },
            ...
          }
        },
        ...
      ]
    }
    """
    if config is None:
        config = {}
    if plate_id_to_filename is None:
        plate_id_to_filename = {}
    
    # Get control rows from config (default: ["N", "O"])
    control_rows = config.get("control_rows", ["N", "O"])
    if not isinstance(control_rows, list):
        control_rows = ["N", "O"]  # Fallback to default if invalid
    
    # Merge default well labels with config labels (config takes precedence)
    # Well labels are now row-based (e.g., "B": "label") instead of well-specific (e.g., "B13": "label")
    # Convert old format to new format for backward compatibility
    config_well_labels = config.get("well_labels", {})
    well_labels = {}
    
    # Process config labels - convert well-specific to row-based if needed
    for key, value in config_well_labels.items():
        if len(key) > 1 and key[0].isalpha() and key[1:].isdigit():
            # Old format: "B13" -> extract row "B"
            row_letter = key[0]
            # Use the last non-empty value for each row (in case of duplicates)
            if row_letter not in well_labels:
                well_labels[row_letter] = value
            elif not well_labels[row_letter] or not well_labels[row_letter].strip():
                well_labels[row_letter] = value
        else:
            # New format: "B" -> use as is
            well_labels[key] = value
    
    # Also convert DEFAULT_WELL_LABELS to row-based format
    default_row_labels = {}
    for key, value in DEFAULT_WELL_LABELS.items():
        if len(key) > 1 and key[0].isalpha():
            row_letter = key[0]
            if row_letter not in default_row_labels:
                default_row_labels[row_letter] = value
    
    # Merge defaults with config (config takes precedence)
    well_labels = {**default_row_labels, **well_labels}
    
    # Parse column groups from filenames
    column_groups: Dict[str, List[str]] = {}  # column_group -> list of plate_ids
    plate_column_info: Dict[str, Dict[str, str]] = {}  # plate_id -> column_info
    
    for plate_id, df in all_plates.items():
        filename = plate_id_to_filename.get(plate_id, plate_id + ".xlsx")
        column_info = parse_column_group_from_filename(filename)
        column_group = column_info["column_group"]
        
        plate_column_info[plate_id] = column_info
        
        if column_group not in column_groups:
            column_groups[column_group] = []
        if plate_id not in column_groups[column_group]:
            column_groups[column_group].append(plate_id)
    
    # Process triplicate groups - row letters only, applies to all columns
    # Config format: {"name": "...", "wells": ["B", "C", "D"]}
    # Column field is ignored if present (for backward compatibility)
    config_triplicate_groups = config.get("triplicate_groups", [])
    
    if config_triplicate_groups and len(config_triplicate_groups) > 0:
        # Process triplicate groups - accept row letters only (e.g., ["B", "C", "D"])
        # Also support backward compatibility with full well IDs (e.g., ["B13", "C13", "D13"])
        triplicate_groups = []
        for group in config_triplicate_groups:
            wells = group.get("wells", [])
            # Extract row letters from wells (handle both formats)
            row_letters = []
            for w in wells:
                # If it's a full well ID (e.g., "B13"), extract row letter
                if len(w) > 1 and w[0].isalpha() and w[1:].isdigit():
                    row_letter = w[0]
                # If it's just a row letter (e.g., "B"), use as is
                elif len(w) == 1 and w.isalpha():
                    row_letter = w
                else:
                    continue  # Skip invalid entries
                
                # Filter out control wells
                if not is_control_well(row_letter, control_rows):
                    if row_letter not in row_letters:
                        row_letters.append(row_letter)
            
            if row_letters:
                # Store row letters only - applies to all columns
                triplicate_groups.append({
                    "name": group.get("name", ""),
                    "wells": sorted(row_letters)  # Store just row letters
                })
        
        if triplicate_groups:
            print(f"  Using {len(triplicate_groups)} triplicate group(s) from config (row letters only, applies to all columns)", flush=True)
    else:
        # Auto-generate or use defaults
        # Collect all wells from all plates to auto-generate triplicate groups
        content_to_wells: Dict[str, List[str]] = {}
        
        for plate_id, df in all_plates.items():
            for well_id, g in df.groupby("well"):
                g_sorted = g.sort_values("time_s")
                content = g_sorted["content"].iloc[0]
                if content and str(content).strip():
                    if content not in content_to_wells:
                        content_to_wells[content] = []
                    if well_id not in content_to_wells[content]:
                        content_to_wells[content].append(well_id)
        
        # Auto-generate triplicate groups: group wells with the same content
        # Extract row pattern (e.g., BCD, EFG) and apply to all columns
        # Exclude control wells from triplicate groups
        auto_triplicate_groups = []
        for content, wells in content_to_wells.items():
            # Filter out control wells
            non_control_wells = [w for w in wells if not is_control_well(w, control_rows)]
            if len(non_control_wells) >= 2:  # At least 2 non-control wells with same content
                # Extract row letters (e.g., ["B13", "C13", "D13"] -> ["B", "C", "D"])
                row_letters = sorted(set([w[0] for w in non_control_wells if w and w[0].isalpha()]))
                
                # If we have a pattern of 3 consecutive rows, create a triplicate group
                # Check for patterns like BCD, EFG, HIJ, etc.
                if len(row_letters) >= 3:
                    # Sort row letters and check if they form a consecutive pattern
                    row_letters_sorted = sorted(row_letters)
                    # Create group with row pattern (applies to all columns)
                    auto_triplicate_groups.append({
                        "name": content,
                        "wells": row_letters_sorted[:3]  # Store row letters only
                    })
                elif len(row_letters) >= 2:
                    # For 2 wells, still create a group
                    auto_triplicate_groups.append({
                        "name": content,
                        "wells": sorted(row_letters)[:2]  # Store row letters only
                    })
        
        # Sort triplicate groups by name
        auto_triplicate_groups.sort(key=lambda x: x.get("name", ""))
        
        if auto_triplicate_groups:
            triplicate_groups = auto_triplicate_groups
            print(f"  Auto-generated {len(triplicate_groups)} triplicate group(s) from data (applies to all columns)", flush=True)
        else:
            # Fall back to defaults - already in row letter format
            # Exclude control wells from default groups
            triplicate_groups = []
            for group in DEFAULT_TRIPLICATE_GROUPS:
                wells = group.get("wells", [])
                # Filter out control wells
                row_letters = [w for w in wells if not is_control_well(w, control_rows)]
                if row_letters:
                    triplicate_groups.append({
                        "name": group.get("name", ""),
                        "wells": sorted(row_letters)  # Already row letters
                    })
            if triplicate_groups:
                print(f"  Using {len(triplicate_groups)} default triplicate group(s) (applies to all columns, control wells excluded)", flush=True)
    
    # Detect duplicate well IDs (same well ID appears in multiple files)
    # This means duplicate columns in the plate grid
    well_id_to_files: Dict[str, List[str]] = {}  # well_id -> list of filenames
    
    for plate_id, df in all_plates.items():
        filename = plate_id_to_filename.get(plate_id, plate_id + ".xlsx")
        # Get all unique well IDs from this plate
        unique_wells = df["well"].unique()
        for well_id in unique_wells:
            well_id_str = str(well_id).strip()
            if well_id_str not in well_id_to_files:
                well_id_to_files[well_id_str] = []
            if filename not in well_id_to_files[well_id_str]:
                well_id_to_files[well_id_str].append(filename)
    
    # Group duplicates by file pairs for cleaner messaging
    duplicate_column_warnings = []
    file_pairs_to_wells: Dict[str, List[str]] = {}  # "file1|file2" -> [well_ids]
    
    for well_id, files in well_id_to_files.items():
        if len(files) > 1:
            # Sort files alphabetically for consistent ordering
            sorted_files = sorted(files)
            file_key = "|".join(sorted_files)
            if file_key not in file_pairs_to_wells:
                file_pairs_to_wells[file_key] = []
            file_pairs_to_wells[file_key].append(well_id)
    
    # Build warnings grouped by file pairs
    for file_key, well_ids in file_pairs_to_wells.items():
        files = file_key.split("|")
        sorted_wells = sorted(well_ids)
        duplicate_column_warnings.append({
            "files": files,
            "well_ids": sorted_wells
        })
        
        # Console message - shorter version
        files_str = " and ".join(files)
        wells_str = ", ".join(sorted_wells[:5])  # Show first 5
        if len(sorted_wells) > 5:
            wells_str += f", ... ({len(sorted_wells)} total)"
        print(f"  ⚠️  WARNING: {files_str} have duplicate wells: {wells_str}", flush=True)
        print(f"     Only data from '{files[0]}' (first alphabetically) is used.", flush=True)
    
    data = {
        "config": {
            "well_labels": well_labels,
            "triplicate_groups": triplicate_groups,
            "column_groups": column_groups,
            "control_rows": control_rows,
            "warnings": {
                "duplicate_columns": duplicate_column_warnings
            }
        },
        "plates": []
    }
    
    for plate_id, df in all_plates.items():
        wells_dict: Dict[str, Any] = {}
        for well_id, g in df.groupby("well"):
            g_sorted = g.sort_values("time_s")
            content = g_sorted["content"].iloc[0]
            # Use row-based label from config if available, otherwise use content
            # Extract row letter from well_id (e.g., "B13" -> "B")
            row_letter = well_id[0] if well_id and well_id[0].isalpha() else ""
            label = well_labels.get(row_letter, content) if row_letter else content
            wells_dict[well_id] = {
                "content": content,
                "label": label,
                "time_s": g_sorted["time_s"].tolist(),
                "values": g_sorted["value"].tolist(),
            }
        
        filename = plate_id_to_filename.get(plate_id, plate_id + ".xlsx")
        column_info = plate_column_info.get(plate_id, parse_column_group_from_filename(filename))
        
        data["plates"].append(
            {
                "id": plate_id,
                "file": filename,
                "column_group": column_info["column_group"],
                "column_info": {
                    "protein": column_info["protein"],
                    "genotype": column_info["genotype"],
                    "buffer": column_info["buffer"],
                    "scientist": column_info["scientist"],
                    "number": column_info["number"]
                },
                "wells": wells_dict,
            }
        )

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Plate Viewer</title>
  <style>
    body {
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      margin: 0;
      padding: 0;
      display: flex;
      flex-direction: column;
      height: 100vh;
      box-sizing: border-box;
    }
    #content-wrapper {
      display: grid;
      grid-template-columns: 320px 1fr;
      flex: 1;
      overflow: hidden;
    }
    aside {
      border-right: 1px solid #ccc;
      padding: 16px;
      box-sizing: border-box;
      background: #f7f7f7;
      overflow-y: auto;
    }
    main {
      padding: 16px;
      box-sizing: border-box;
      overflow: auto;
    }
    h1 {
      font-size: 1.2rem;
      margin-top: 0;
    }
    label {
      display: block;
      margin-bottom: 8px;
      font-weight: 600;
    }
    select {
      width: 100%%;
      padding: 8px;
      border: 1px solid #ddd;
      border-radius: 4px;
      background: white;
      font-size: 0.9rem;
      margin-bottom: 16px;
    }
    #well-info {
      margin-top: 8px;
      font-size: 0.9rem;
      padding: 8px;
      background: white;
      border-radius: 4px;
      border: 1px solid #ddd;
    }
    #plate-grid {
      display: grid;
      grid-template-columns: auto repeat(%(cols)d, 1fr);
      grid-auto-rows: 32px;
      gap: 2px;
      font-size: 11px;
      margin-bottom: 16px;
    }
    .grid-header {
      font-weight: 700;
      text-align: center;
    }
    .row-label, .col-label {
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .well {
      border: 1px solid #ddd;
      border-radius: 3px;
      display: flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      background: #fff;
      position: relative;
      font-size: 10px;
      pointer-events: auto;
      z-index: 1;
    }
    .well[data-has-data="true"]:not(.triplicate-group):not(.control) {
      background: #e4f0ff;
      border-color: #6c9cff;
    }
    .well.selected {
      outline: 2px solid #2a6cff;
      outline-offset: -2px;
      background: #e4f0ff;
    }
    .well.replicate {
      border-style: dashed;
    }
    .well.triplicate-group {
      /* Borders are added via inline styles only when selected */
    }
    /* Individual colors will be applied via inline styles */
    .well.control {
      background: #e8f5e9;
      border: 1px solid #ddd;
    }
    .well.control[data-has-data="true"] {
      background: #c8e6c9;
      border-color: #6c9cff;
    }
    .well.control.selected {
      background: #c8e6c9;
      outline: 2px solid #2a6cff;
      outline-offset: -2px;
    }
    .well-label {
      font-size: 8px;
      position: absolute;
      bottom: 1px;
      left: 1px;
      right: 1px;
      text-align: center;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      pointer-events: none;
      z-index: 2;
    }
    #chart-container {
      margin-top: 12px;
      height: 800px;
      min-height: 800px;
      position: relative;
    }
    canvas {
      max-width: 100%%;
      height: 800px !important;
    }
    .btn {
      padding: 6px 12px;
      margin: 4px 0;
      border: 1px solid #ccc;
      border-radius: 4px;
      background: white;
      cursor: pointer;
      font-size: 0.9rem;
    }
    .btn:hover {
      background: #f0f0f0;
    }
    .btn-primary {
      background: #2a6cff;
      color: white;
      border-color: #2a6cff;
    }
    .btn-primary:hover {
      background: #1a5cff;
    }
    #warning-banner {
      background: #fff3cd;
      border: 2px solid #ffc107;
      border-radius: 4px;
      padding: 12px 16px;
      margin: 16px;
      margin-bottom: 0;
      font-size: 0.9rem;
      color: #856404;
      display: none;
    }
    #warning-banner.show {
      display: block;
    }
    #warning-banner h3 {
      margin: 0 0 8px 0;
      font-size: 1rem;
      color: #856404;
    }
    #warning-banner ul {
      margin: 8px 0 0 0;
      padding-left: 20px;
    }
    #warning-banner li {
      margin: 4px 0;
    }
    input[type="checkbox"]:disabled {
      opacity: 0.5;
      cursor: not-allowed !important;
    }
    input[type="checkbox"]:disabled + span {
      opacity: 0.6;
      color: #999;
    }
  </style>
  <!-- Chart.js via CDN -->
  <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
</head>
<body>
  <div id="warning-banner"></div>
  <div id="content-wrapper">
  <aside>
    <h1>Plate Viewer</h1>
    <label>Select Dataset Group:</label>
    <select id="dataset-select"></select>
    <div id="dataset-info" style="margin-top: 16px; padding: 12px; background: white; border-radius: 4px; border: 1px solid #ddd; font-size: 0.85rem;">
      <div style="font-weight: 600; margin-bottom: 8px;">Dataset Information</div>
      <div id="dataset-details">Select a dataset to view details.</div>
    </div>
    <div id="column-legend" style="margin-top: 16px; padding: 12px; background: white; border-radius: 4px; border: 1px solid #ddd; font-size: 0.85rem;">
      <div style="font-weight: 600; margin-bottom: 8px;">Select Datasets by Column</div>
      <div id="column-legend-content" style="font-size: 0.8rem; color: #666;">No column information available.</div>
    </div>
    <div id="well-selector" style="margin-top: 16px; padding: 12px; background: white; border-radius: 4px; border: 1px solid #ddd;">
      <div style="display: flex; align-items: center; cursor: pointer; font-weight: 600; margin-bottom: 8px;" onclick="toggleWellSelector()">
        <span id="well-selector-arrow" style="margin-right: 8px;">▼</span>
        <span>Select Wells (click to exclude)</span>
      </div>
      <div id="well-selector-grid-container" style="display: none;">
        <div id="well-selector-grid" style="display: grid; grid-template-columns: auto repeat(24, 1fr); grid-auto-rows: 18px; gap: 1px; font-size: 8px; margin-top: 8px;"></div>
        <div style="font-size: 0.75rem; color: #666; margin-top: 16px; text-align: center;">
          Click wells to exclude them (✕ = excluded)
        </div>
      </div>
    </div>
    <div style="margin-top: 16px; padding: 12px; background: white; border-radius: 4px; border: 1px solid #ddd;">
      <div style="font-weight: 600; margin-bottom: 8px; font-size: 0.95rem; color: #333;">Data Processing Workflow</div>
      <div style="font-size: 0.75rem; color: #666; margin-bottom: 12px; padding: 8px; background: #f9f9f9; border-radius: 3px;">
        <strong>Step 1:</strong> Normalize using fitted baseline (if enabled)<br>
        <strong>Step 2:</strong> Group wells into triplicate groups (if enabled)<br>
        <strong>Step 3:</strong> Cutoff data points before first timepoint and after time point (if enabled)<br>
        <strong>Step 4:</strong> Fit regression curve to data points (if enabled)
      </div>
      <label style="display: flex; align-items: center; cursor: pointer; font-size: 0.9rem;">
        <input type="checkbox" id="normalize-baseline-toggle" style="margin-right: 8px; cursor: pointer;" onchange="updateChart();">
        <span><strong>Step 1:</strong> Normalize using fitted baseline</span>
      </label>
      <div style="font-size: 0.75rem; color: #666; margin-top: 4px; margin-bottom: 8px; margin-left: 24px;">
        Fit a baseline function (LOWESS, constant, or polynomial) on the baseline window and normalize each well using the fitted baseline. Select method and parameters in the "Baseline Fitting & Normalization Options" section below.
      </div>
      <div id="normalize-baseline-parameters" style="font-size: 0.8rem; padding: 8px; background: #f9f9f9; border-radius: 3px; display: none; margin-left: 24px; margin-bottom: 8px;">
        <div style="margin-bottom: 8px;">
          <label style="display: flex; align-items: center; margin-bottom: 4px;">
            <span style="min-width: 120px;">Baseline window end (s):</span>
            <input type="number" id="normalize-baseline-end-time" value="24" min="0" step="0.1" style="width: 80px; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-family: monospace;" onchange="updateChart();">
          </label>
        </div>
        <div style="font-size: 0.75rem; color: #666; margin-top: 4px;">
          Baseline fitting method and normalization mode are selected in the "Baseline Fitting & Normalization Options" section below.
        </div>
        <div id="current-settings-display" style="font-size: 0.75rem; color: #333; margin-top: 8px; padding: 6px; background: #f0f0f0; border-radius: 3px; border-left: 3px solid #2a6cff;">
          <div style="font-weight: 600; margin-bottom: 4px; color: #2a6cff;">Current Settings:</div>
          <div id="current-settings-content" style="font-size: 0.7rem; color: #666;">
            Baseline: <span id="current-baseline-method">Constant</span> | Normalization: <span id="current-normalization-mode">Multiplicative</span> | Curve Fit: <span id="current-regression-type">Mono-exponential plateau</span>
          </div>
        </div>
      </div>
      <label style="display: flex; align-items: center; cursor: pointer; font-size: 0.9rem; margin-top: 12px;">
        <input type="checkbox" id="group-triplicates-toggle" style="margin-right: 8px; cursor: pointer;" onchange="updateChart();">
        <span><strong>Step 2:</strong> Group wells into triplicate groups</span>
      </label>
      <div style="font-size: 0.75rem; color: #666; margin-top: 4px; margin-left: 24px;">
        When enabled, wells that are part of triplicate groups (e.g., rows B, C, D) are grouped together. When disabled, each well is shown individually. <strong>Note:</strong> Requires Step 1 (normalization) to be enabled.
      </div>
      <label style="display: flex; align-items: center; cursor: pointer; font-size: 0.9rem; margin-top: 12px;">
        <input type="checkbox" id="cutoff-baseline-toggle" style="margin-right: 8px; cursor: pointer;" onchange="updateChart();">
        <span><strong>Step 3:</strong> Cutoff data points before first timepoint and after time point</span>
      </label>
      <div style="font-size: 0.75rem; color: #666; margin-top: 4px; margin-left: 24px;">
        When enabled, data points before the first real timepoint and after the cutoff time will be excluded from the graph display.
      </div>
      <div id="cutoff-baseline-parameters" style="font-size: 0.8rem; padding: 8px; background: #f9f9f9; border-radius: 3px; display: none; margin-left: 24px; margin-bottom: 8px;">
        <div style="margin-bottom: 8px;">
          <label style="display: flex; align-items: center; margin-bottom: 4px;">
            <span style="min-width: 120px;">Cutoff time (s):</span>
            <input type="number" id="cutoff-baseline-time" value="360" min="0" step="1" style="width: 80px; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-family: monospace;" onchange="updateChart();">
          </label>
        </div>
        <div style="font-size: 0.75rem; color: #666; margin-top: 4px;">
          Data points before the first real timepoint and after this time will be excluded (default: 6 minutes = 360 seconds).
        </div>
      </div>
      <label style="display: flex; align-items: center; cursor: pointer; font-size: 0.9rem; margin-top: 12px;">
        <input type="checkbox" id="fit-regression-toggle" style="margin-right: 8px; cursor: pointer;" onchange="updateRegressionTypeParams(true); updateChart();">
        <span><strong>Step 4:</strong> Fit regression curve to data points</span>
      </label>
      <div style="font-size: 0.75rem; color: #666; margin-top: 4px; margin-left: 24px;">
        When enabled, a regression curve will be fitted to the data points and displayed on the graph.
      </div>
      <div id="fit-regression-parameters" style="font-size: 0.8rem; padding: 8px; background: #f9f9f9; border-radius: 3px; display: none; margin-left: 24px; margin-bottom: 8px;">
        <div style="margin-bottom: 8px;">
          <label style="display: flex; align-items: center; margin-bottom: 4px;">
            <span style="min-width: 140px;">Regression type:</span>
            <select id="regression-type-select" style="flex: 1; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-size: 0.85rem;" onchange="updateRegressionTypeParams();">
              <option value="linear">Linear</option>
              <option value="polynomial">Polynomial</option>
              <option value="exponential">Exponential</option>
              <option value="logarithmic">Logarithmic</option>
              <option value="mono-exponential" selected>Mono-exponential (rise to plateau)</option>
            </select>
          </label>
        </div>
        <div id="regression-polynomial-order" style="margin-bottom: 8px; display: none;">
          <label style="display: flex; align-items: center; margin-bottom: 4px;">
            <span style="min-width: 140px;">Polynomial order:</span>
            <input type="number" id="polynomial-order-input" value="2" min="2" max="5" step="1" style="width: 80px; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-family: monospace;" onchange="updateChart();">
          </label>
        </div>
        <div style="font-size: 0.75rem; color: #666; margin-top: 4px;">
          The regression curve will be overlaid on the data points.
        </div>
        <div id="regression-info-panel" style="margin-top: 12px; padding: 8px; background: #e8f4fd; border-radius: 3px; border: 1px solid #2a6cff; display: none;">
          <div style="font-weight: 600; margin-bottom: 8px; color: #2a6cff;">Mono-exponential Fit Parameters</div>
          <div id="regression-info-content" style="font-size: 0.75rem; font-family: monospace;">
            <!-- Parameters will be populated here -->
          </div>
        </div>
      </div>
    </div>
    <div style="margin-top: 16px; padding: 12px; background: white; border-radius: 4px; border: 1px solid #ddd;">
      <div style="display: flex; align-items: center; cursor: pointer; font-weight: 600; margin-bottom: 8px; font-size: 0.95rem; color: #333;" onclick="toggleBaselineOptions()">
        <span id="baseline-options-arrow" style="margin-right: 8px;">▼</span>
        <span>Baseline Fitting & Normalization Options</span>
      </div>
      <div id="baseline-options-content" style="display: none; font-size: 0.8rem;">
        <div style="margin-bottom: 16px; padding: 8px; background: #e8f4fd; border-radius: 3px; border: 1px solid #2a6cff;">
          <div style="font-weight: 600; margin-bottom: 8px; color: #2a6cff;">Select Baseline Fitting Method</div>
          <div style="margin-bottom: 8px;">
            <label style="display: flex; align-items: center; margin-bottom: 4px; font-size: 0.85rem;">
              <span style="min-width: 140px;">Method:</span>
              <select id="baseline-method-select" style="flex: 1; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-size: 0.85rem;" onchange="updateBaselineMethodParams();">
                <option value="lowess">LOWESS</option>
                <option value="constant" selected>Constant</option>
                <option value="polynomial">Polynomial</option>
              </select>
            </label>
          </div>
          <div id="baseline-lowess-params" style="margin-left: 12px; margin-bottom: 4px; display: none;">
            <label style="display: flex; align-items: center; margin-bottom: 4px; font-size: 0.85rem;">
              <span style="min-width: 140px;">Smoothing fraction (frac):</span>
              <input type="number" id="baseline-frac-input" value="0.5" min="0.1" max="1.0" step="0.1" style="flex: 1; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-family: monospace; font-size: 0.85rem;">
            </label>
            <div style="font-size: 0.75rem; color: #666; margin-left: 140px; margin-top: 2px;">
              Higher = smoother (0.7-0.9), Lower = more local (0.2-0.4)
            </div>
          </div>
          <div id="baseline-polynomial-params" style="margin-left: 12px; margin-bottom: 4px; display: none;">
            <label style="display: flex; align-items: center; margin-bottom: 4px; font-size: 0.85rem;">
              <span style="min-width: 140px;">Polynomial order:</span>
              <input type="number" id="baseline-poly-order-input" value="1" min="1" max="5" step="1" style="flex: 1; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-family: monospace; font-size: 0.85rem;">
            </label>
            <div style="font-size: 0.75rem; color: #666; margin-left: 140px; margin-top: 2px;">
              1 = linear, 2 = quadratic, 3+ = higher order
            </div>
          </div>
        </div>
        <div style="margin-bottom: 16px; padding: 8px; background: #e8f4fd; border-radius: 3px; border: 1px solid #2a6cff;">
          <div style="font-weight: 600; margin-bottom: 8px; color: #2a6cff;">Select Normalization Method</div>
          <div>
            <label style="display: flex; align-items: center; margin-bottom: 4px; font-size: 0.85rem;">
              <span style="min-width: 140px;">Normalization mode:</span>
              <select id="normalization-mode-select" style="flex: 1; padding: 4px; border: 1px solid #ddd; border-radius: 3px; font-size: 0.85rem;" onchange="updateCurrentSettingsDisplay();">
                <option value="delta_f_over_f">ΔF/F (Delta F over F)</option>
                <option value="multiplicative" selected>Multiplicative</option>
              </select>
            </label>
          </div>
        </div>
        <div style="margin-bottom: 16px;">
          <div style="font-weight: 600; margin-bottom: 8px; color: #2a6cff;">Baseline Fitting Methods</div>
          <div style="margin-bottom: 12px; padding: 8px; background: #f9f9f9; border-radius: 3px;">
            <div style="font-weight: 600; margin-bottom: 4px;">1. LOWESS (Locally Weighted Scatterplot Smoothing)</div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Formula:</strong> g(t) = LOWESS(F_baseline, t_baseline, frac)
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>When to use:</strong> Baseline has gradual drift or non-linear trends
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Best for:</strong> Non-linear baselines with smooth trends
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Parameters:</strong>
              <ul style="margin: 4px 0 0 20px; padding-left: 0;">
                <li><strong>frac:</strong> Fraction of data used for smoothing (default: 0.5)
                  <ul style="margin: 2px 0 0 20px; padding-left: 0;">
                    <li>Higher frac (0.7-0.9) = smoother, more global fit</li>
                    <li>Lower frac (0.2-0.4) = more local, follows data closely</li>
                  </ul>
                </li>
              </ul>
            </div>
          </div>
          <div style="margin-bottom: 12px; padding: 8px; background: #f9f9f9; border-radius: 3px;">
            <div style="font-weight: 600; margin-bottom: 4px;">2. CONSTANT (0th-order polynomial)</div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Formula:</strong> g(t) = mean(F_baseline) for all t
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>When to use:</strong> Baseline is stable/flat with minimal drift
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Best for:</strong> Stable baselines with no time-dependent trends
            </div>
            <div style="margin-left: 12px;">
              <strong>Parameters:</strong> None
            </div>
          </div>
          <div style="margin-bottom: 12px; padding: 8px; background: #f9f9f9; border-radius: 3px;">
            <div style="font-weight: 600; margin-bottom: 4px;">3. POLYNOMIAL (1st-order or higher)</div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Formula:</strong> g(t) = c₀ + c₁t + c₂t² + ...
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>When to use:</strong> Baseline has linear or polynomial drift
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Best for:</strong> Baselines with known polynomial trends
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Parameters:</strong>
              <ul style="margin: 4px 0 0 20px; padding-left: 0;">
                <li><strong>poly_order:</strong> Polynomial order (default: 1 = linear)
                  <ul style="margin: 2px 0 0 20px; padding-left: 0;">
                    <li>Order 1: Linear drift (g(t) = c₀ + c₁t)</li>
                    <li>Order 2: Quadratic drift (g(t) = c₀ + c₁t + c₂t²)</li>
                    <li>Higher orders: More complex trends</li>
                  </ul>
                </li>
              </ul>
            </div>
          </div>
        </div>
        <div style="margin-bottom: 16px;">
          <div style="font-weight: 600; margin-bottom: 8px; color: #2a6cff;">Normalization Methods</div>
          <div style="margin-bottom: 12px; padding: 8px; background: #f9f9f9; border-radius: 3px;">
            <div style="font-weight: 600; margin-bottom: 4px;">1. DELTA_F_OVER_F (ΔF/F)</div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Formula:</strong> F'(t) = (F(t) - g(t)) / g(t)
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>When to use:</strong> Measure relative change from baseline
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Best for:</strong> Detecting increases/decreases relative to baseline
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Interpretation:</strong>
              <ul style="margin: 4px 0 0 20px; padding-left: 0;">
                <li>F'(t) = 0: No change from baseline</li>
                <li>F'(t) > 0: Increase above baseline (e.g., 0.5 = 50%% increase)</li>
                <li>F'(t) < 0: Decrease below baseline (e.g., -0.3 = 30%% decrease)</li>
              </ul>
            </div>
            <div style="margin-left: 12px;">
              Values can be negative (below baseline) or positive (above baseline)
            </div>
          </div>
          <div style="margin-bottom: 12px; padding: 8px; background: #f9f9f9; border-radius: 3px;">
            <div style="font-weight: 600; margin-bottom: 4px;">2. MULTIPLICATIVE</div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Formula:</strong> F_norm(t) = F(t) / g(t)
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>When to use:</strong> Normalize by baseline without subtracting
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Best for:</strong> Fold-change analysis
            </div>
            <div style="margin-left: 12px; margin-bottom: 4px;">
              <strong>Interpretation:</strong>
              <ul style="margin: 4px 0 0 20px; padding-left: 0;">
                <li>F_norm(t) = 1.0: No change from baseline</li>
                <li>F_norm(t) = 2.0: 2-fold increase (100%% increase)</li>
                <li>F_norm(t) = 0.5: 2-fold decrease (50%% of baseline)</li>
              </ul>
            </div>
            <div style="margin-left: 12px;">
              Values are always positive (ratio to baseline)
            </div>
          </div>
        </div>
        <div>
          <div style="font-weight: 600; margin-bottom: 8px; color: #2a6cff;">Statistics Computation</div>
          <div style="padding: 8px; background: #f9f9f9; border-radius: 3px;">
            <div style="margin-bottom: 8px;">
              For each condition and timepoint, the following statistics are computed:
            </div>
            <div style="margin-bottom: 6px;">
              <strong>Mean:</strong> μ(t) = (1/n) × Σ F'(t)
              <div style="margin-left: 12px; font-size: 0.75rem; color: #666;">
                Average normalized value across all wells in the condition
              </div>
            </div>
            <div style="margin-bottom: 6px;">
              <strong>SD:</strong> σ(t) = sqrt(Σ(F'(t) - μ(t))² / (n-1))
              <div style="margin-left: 12px; font-size: 0.75rem; color: #666;">
                Sample standard deviation (ddof=1) across wells
              </div>
            </div>
            <div>
              <strong>SEM:</strong> SEM(t) = σ(t) / sqrt(n)
              <div style="margin-left: 12px; font-size: 0.75rem; color: #666;">
                Standard error of the mean = SD / sqrt(n_wells)
              </div>
              <div style="margin-left: 12px; font-size: 0.75rem; color: #666; margin-top: 4px;">
                <strong>Note:</strong> SEM is computed across wells, not propagated from baseline error
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
    <button class="btn btn-primary" onclick="updateChart()" style="margin-top: 16px;">Update Chart</button>
    <button class="btn" onclick="clearSelection()">Clear Selection</button>
    <div id="well-info" style="margin-top: 16px;">Click wells on the plate to add them to the chart.</div>
  </aside>
  <main>
    <div id="plate-grid"></div>
    <div id="chart-container">
      <canvas id="well-chart"></canvas>
    </div>
  </main>
  </div>

<script>
const PLATE_ROWS = "%(rows)s".split("");
const PLATE_COLS = %(cols)d;

let viewerData = null;
let chart = null;
let selectedWells = new Set(); // Set of {plateId:wellId} strings
let allWellsData = {}; // Combined wells from all plates: {plateId: {wellId: {...}}}
let selectedDatasetId = null; // Currently selected dataset ID
let enabledDatasets = new Set(); // Set of enabled dataset IDs (empty = all enabled)
let enabledWells = new Set(); // Set of enabled well IDs (empty = all enabled)
let enabledColumns = new Set(); // Set of enabled column numbers (empty = all enabled)

// Register error bar plugin globally
Chart.register({
  id: "errorBars",
  afterDatasetsDraw: (chart) => {
    const ctx = chart.ctx;
    chart.data.datasets.forEach((dataset, datasetIndex) => {
      if (!dataset.data || dataset.data.length === 0) return;
      
      // Check if any data point has error information
      const hasErrors = dataset.data.some(d => d && d.error !== undefined && d.error !== null);
      if (!hasErrors) return;
      
      const meta = chart.getDatasetMeta(datasetIndex);
      ctx.save();
      ctx.strokeStyle = dataset.borderColor || "rgba(0,0,0,0.5)";
      ctx.lineWidth = 1;
      
      meta.data.forEach((point, index) => {
        if (!point || !dataset.data[index] || dataset.data[index].error === undefined || dataset.data[index].error === null) return;
        
        const x = point.x;
        const y = point.y;
        const error = dataset.data[index].error;
        const yScale = chart.scales.y;
        const errorPixels = Math.abs(yScale.getPixelForValue(y + error) - yScale.getPixelForValue(y));
        
        // Draw vertical error bar
        ctx.beginPath();
        ctx.moveTo(x, y - errorPixels);
        ctx.lineTo(x, y + errorPixels);
        ctx.stroke();
        
        // Draw horizontal caps
        const capWidth = 4;
        ctx.beginPath();
        ctx.moveTo(x - capWidth, y - errorPixels);
        ctx.lineTo(x + capWidth, y - errorPixels);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(x - capWidth, y + errorPixels);
        ctx.lineTo(x + capWidth, y + errorPixels);
        ctx.stroke();
      });
      
      ctx.restore();
    });
  }
});

// Register baseline label plugin to display baseline numbers on the chart
Chart.register({
  id: "baselineLabels",
  afterDatasetsDraw: (chart) => {
    const ctx = chart.ctx;
    const xScale = chart.scales.x;
    const yScale = chart.scales.y;
    
    chart.data.datasets.forEach((dataset, datasetIndex) => {
      // Only process baseline lines (identified by having exactly 2 points and order === 1000)
      if (!dataset.data || dataset.data.length !== 2 || dataset.order !== 1000) return;
      
      // Extract baseline value from label (format: "label (baseline: value)")
      const baselineMatch = dataset.label.match(/\(baseline:\s*([\d.]+)\)/);
      if (!baselineMatch) return;
      
      const baselineValue = parseFloat(baselineMatch[1]);
      if (isNaN(baselineValue)) return;
      
      // Get the y position of the baseline
      const baselineY = dataset.data[0].y;
      const yPixel = yScale.getPixelForValue(baselineY);
      
      // Position text at the right edge of the chart, slightly offset from the line
      const chartArea = chart.chartArea;
      const textX = chartArea.right - 10; // 10px from right edge
      const textY = yPixel - 5; // 5px above the line
      
      ctx.save();
      ctx.font = "12px Arial";
      ctx.fillStyle = dataset.borderColor || "rgba(128, 128, 128, 0.7)";
      ctx.textAlign = "right";
      ctx.textBaseline = "bottom";
      
      // Draw background rectangle for better readability
      const text = baselineValue.toFixed(2);
      const textMetrics = ctx.measureText(text);
      const padding = 4;
      ctx.fillStyle = "rgba(255, 255, 255, 0.8)";
      ctx.fillRect(
        textX - textMetrics.width - padding,
        textY - 12 - padding,
        textMetrics.width + padding * 2,
        12 + padding * 2
      );
      
      // Draw the text
      ctx.fillStyle = dataset.borderColor || "rgba(128, 128, 128, 0.9)";
      ctx.fillText(text, textX, textY);
      
      ctx.restore();
    });
  }
});

function displayWarnings() {
  const warningBanner = document.getElementById("warning-banner");
  if (!warningBanner || !viewerData || !viewerData.config) {
    return;
  }
  
  const warnings = viewerData.config.warnings || {};
  const duplicateColumns = warnings.duplicate_columns || [];
  
  if (duplicateColumns.length === 0) {
    warningBanner.classList.remove("show");
    return;
  }
  
  // Build warning message - shorter version
  let html = '<h3>⚠️ Warning: Duplicate Wells Detected</h3>';
  
  duplicateColumns.forEach(warning => {
    const files = warning.files || [];
    const wellIds = warning.well_ids || [];
    const filesStr = files.join(" and ");
    const wellsStr = wellIds.length <= 10 
      ? wellIds.join(", ")
      : wellIds.slice(0, 10).join(", ") + `, ... (${wellIds.length} total)`;
    
    html += `<div style="margin-bottom: 8px;"><strong>${filesStr}</strong> have duplicate wells: ${wellsStr}</div>`;
    html += `<div style="margin-left: 16px; font-size: 0.85rem; color: #666; margin-bottom: 12px;">Only data from '<strong>${files[0]}</strong>' (first alphabetically) is used.</div>`;
  });
  
  warningBanner.innerHTML = html;
  warningBanner.classList.add("show");
}

async function loadData() {
  const res = await fetch("viewer_data.json");
  viewerData = await res.json();
  combineAllWells();
  // Initialize triplicate group colors before rendering
  initializeTriplicateGroupColors();
  displayWarnings();
  initDatasetSelect();
  renderPlateGrid();
}

function combineAllWells() {
  // Combine all wells from all plates into a single structure
  allWellsData = {};
  viewerData.plates.forEach(plate => {
    allWellsData[plate.id] = plate.wells;
  });
}

function initDatasetSelect() {
  const select = document.getElementById("dataset-select");
  select.innerHTML = "";
  
  // Group plates by column_group
  const platesByGroup = {};
  viewerData.plates.forEach(plate => {
    const group = plate.column_group || "Other";
    if (!platesByGroup[group]) {
      platesByGroup[group] = [];
    }
    platesByGroup[group].push(plate);
  });
  
  // Sort groups alphabetically
  const sortedGroups = Object.keys(platesByGroup).sort();
  
  // Create options for each group (the group itself is selectable)
  sortedGroups.forEach(groupName => {
    const plates = platesByGroup[groupName];
    
    const option = document.createElement("option");
    // Use special prefix to identify group selections
    option.value = `group:${groupName}`;
    
    // Create friendly name from first plate's column info
    const firstPlate = plates[0];
    const columnInfo = firstPlate.column_info || {};
    const parts = [];
    if (columnInfo.protein) parts.push(columnInfo.protein);
    if (columnInfo.genotype) parts.push(columnInfo.genotype);
    if (columnInfo.buffer) parts.push(columnInfo.buffer);
    
    // Add count if multiple datasets
    if (plates.length > 1) {
      option.textContent = parts.length > 0 
        ? `${parts.join(" - ")} (${plates.length} datasets)`
        : `${groupName} (${plates.length} datasets)`;
    } else {
      option.textContent = parts.length > 0 
        ? parts.join(" - ")
        : groupName;
    }
    
    select.appendChild(option);
  });
  
  // Select first group by default
  if (sortedGroups.length > 0) {
    selectedDatasetId = `group:${sortedGroups[0]}`;
    select.value = selectedDatasetId;
  }
  
  select.addEventListener("change", () => {
    selectedDatasetId = select.value;
    // Clean up selected wells - only keep wells from selected datasets
    const selectedDatasets = getSelectedDatasets();
    const newSelectedWells = new Set();
    selectedWells.forEach(key => {
      const [plateId] = key.split(":");
      if (selectedDatasets.includes(plateId)) {
        newSelectedWells.add(key);
      }
    });
    selectedWells = newSelectedWells;
    // Reset enabled datasets/wells/columns when changing dataset group
    enabledDatasets.clear();
    enabledWells.clear();
    enabledColumns.clear();
    
    renderPlateGrid();
    updateWellInfo();
    updateDatasetInfo();
    renderWellSelectorGrid();
  });
  
  // Update info on initial load
  updateDatasetInfo();
}

function updateDatasetInfo() {
  const infoDiv = document.getElementById("dataset-details");
  if (!infoDiv || !selectedDatasetId) {
    if (infoDiv) {
      infoDiv.textContent = "Select a dataset to view details.";
    }
    return;
  }
  
  const selectedDatasets = getSelectedDatasets();
  if (selectedDatasets.length === 0) {
    infoDiv.textContent = "No dataset selected.";
    return;
  }
  
  // Get metadata from first dataset
  const firstPlate = viewerData.plates.find(p => p.id === selectedDatasets[0]);
  if (!firstPlate) {
    infoDiv.textContent = "Dataset information not available.";
    return;
  }
  
  const columnInfo = firstPlate.column_info || {};
  const allPlates = viewerData.plates.filter(p => selectedDatasets.includes(p.id));
  
  // Build info HTML
  let html = "";
  
  if (columnInfo.protein) {
    html += `<div style="margin-bottom: 4px;"><strong>Protein:</strong> ${columnInfo.protein}</div>`;
  }
  if (columnInfo.genotype) {
    html += `<div style="margin-bottom: 4px;"><strong>Genotype:</strong> ${columnInfo.genotype}</div>`;
  }
  if (columnInfo.buffer) {
    html += `<div style="margin-bottom: 4px;"><strong>Buffer:</strong> ${columnInfo.buffer}</div>`;
  }
  if (columnInfo.scientist) {
    html += `<div style="margin-bottom: 4px;"><strong>Scientist:</strong> ${columnInfo.scientist}</div>`;
  }
  
  html += `<div style="margin-top: 8px; padding-top: 8px; border-top: 1px solid #ddd;">`;
  html += `<div style="margin-bottom: 4px;"><strong>Datasets:</strong> ${allPlates.length}</div>`;
  html += `<div style="font-size: 0.8rem; color: #666;">`;
  allPlates.forEach((plate, idx) => {
    html += `${idx + 1}. ${plate.file || plate.id}`;
    if (plate.column_info && plate.column_info.number) {
      html += ` (${plate.column_info.number})`;
    }
    html += "<br>";
  });
  html += `</div></div>`;
  
  infoDiv.innerHTML = html;
}

function getSelectedDatasets() {
  // Return all datasets in the selected group, filtered by enabledDatasets
  if (!selectedDatasetId) {
    return [];
  }
  
  let datasets = [];
  
  // Check if it's a group selection (format: "group:GROUP_NAME")
  if (selectedDatasetId.startsWith("group:")) {
    const groupName = selectedDatasetId.substring(6); // Remove "group:" prefix
    const selectedPlates = viewerData.plates.filter(p => p.column_group === groupName);
    datasets = selectedPlates.map(p => p.id);
  } else {
    // Otherwise, it's a single dataset selection (legacy support)
    const plate = viewerData.plates.find(p => p.id === selectedDatasetId);
    if (plate) {
      datasets = [selectedDatasetId];
    }
  }
  
  // Filter by enabledDatasets if any are set
  if (enabledDatasets.size > 0) {
    return datasets.filter(id => enabledDatasets.has(id));
  }
  
  return datasets;
}

function toggleWellSelector() {
  const container = document.getElementById("well-selector-grid-container");
  const arrow = document.getElementById("well-selector-arrow");
  if (container.style.display === "none") {
    container.style.display = "block";
    arrow.textContent = "▲";
    renderWellSelectorGrid();
  } else {
    container.style.display = "none";
    arrow.textContent = "▼";
  }
}

function toggleBaselineOptions() {
  const container = document.getElementById("baseline-options-content");
  const arrow = document.getElementById("baseline-options-arrow");
  if (container.style.display === "none") {
    container.style.display = "block";
    arrow.textContent = "▲";
  } else {
    container.style.display = "none";
    arrow.textContent = "▼";
  }
}

function updateBaselineMethodParams() {
  const method = document.getElementById("baseline-method-select").value;
  const lowessParams = document.getElementById("baseline-lowess-params");
  const polynomialParams = document.getElementById("baseline-polynomial-params");
  
  if (method === "lowess") {
    lowessParams.style.display = "block";
    polynomialParams.style.display = "none";
  } else if (method === "polynomial") {
    lowessParams.style.display = "none";
    polynomialParams.style.display = "block";
  } else {
    // constant
    lowessParams.style.display = "none";
    polynomialParams.style.display = "none";
  }
  
  updateCurrentSettingsDisplay();
}

function updateRegressionTypeParams(skipChartUpdate) {
  const regressionType = document.getElementById("regression-type-select").value;
  const polynomialOrderDiv = document.getElementById("regression-polynomial-order");
  const infoPanel = document.getElementById("regression-info-panel");
  const fitRegressionToggle = document.getElementById("fit-regression-toggle");
  
  if (regressionType === "polynomial" && polynomialOrderDiv) {
    polynomialOrderDiv.style.display = "block";
  } else if (polynomialOrderDiv) {
    polynomialOrderDiv.style.display = "none";
  }
  
  // Show/hide info panel for mono-exponential (only if regression toggle is checked)
  if (infoPanel) {
    const isRegressionEnabled = fitRegressionToggle && fitRegressionToggle.checked;
    if (regressionType === "mono-exponential" && isRegressionEnabled) {
      infoPanel.style.display = "block";
    } else {
      infoPanel.style.display = "none";
    }
  }
  
  updateCurrentSettingsDisplay();
  
  // Update chart when regression type changes (unless skipped for initialization)
  if (!skipChartUpdate) {
    updateChart();
  }
}

function updateCurrentSettingsDisplay() {
  const baselineMethodSelect = document.getElementById("baseline-method-select");
  const normalizationModeSelect = document.getElementById("normalization-mode-select");
  const regressionTypeSelect = document.getElementById("regression-type-select");
  
  const baselineMethodEl = document.getElementById("current-baseline-method");
  const normalizationModeEl = document.getElementById("current-normalization-mode");
  const regressionTypeEl = document.getElementById("current-regression-type");
  
  if (baselineMethodEl && baselineMethodSelect) {
    const method = baselineMethodSelect.value;
    baselineMethodEl.textContent = method.charAt(0).toUpperCase() + method.slice(1);
  }
  
  if (normalizationModeEl && normalizationModeSelect) {
    const mode = normalizationModeSelect.value;
    normalizationModeEl.textContent = mode === "delta_f_over_f" ? "ΔF/F" : "Multiplicative";
  }
  
  if (regressionTypeEl && regressionTypeSelect) {
    const type = regressionTypeSelect.value;
    const typeLabels = {
      "linear": "Linear",
      "polynomial": "Polynomial",
      "exponential": "Exponential",
      "logarithmic": "Logarithmic",
      "mono-exponential": "Mono-exponential plateau"
    };
    regressionTypeEl.textContent = typeLabels[type] || type;
  }
}

function parseLabelForLegend(label) {
  // Parse label format: "E17, F17, G17, E14, G14 - 500-10-2 (PLCe_WT_Ca) Col 17 (5-plicate, 2 datasets, mean ± SE)"
  // Extract main heading (e.g., "500-10-2") and return rest as details
  
  if (!label) return { mainHeading: label, details: "" };
  
  // Pattern to match: number-number-number (the main heading like "500-10-2")
  const mainHeadingMatch = label.match(/\b(\d+-\d+-\d+)\b/);
  
  if (mainHeadingMatch) {
    const mainHeading = mainHeadingMatch[1];
    const matchIndex = mainHeadingMatch.index;
    
    // Get parts before and after the main heading
    const before = label.substring(0, matchIndex).trim();
    const after = label.substring(matchIndex + mainHeading.length).trim();
    
    // Combine before and after
    let details = (before + " " + after).trim();
    // Clean up extra spaces and dashes
    details = details.replace(/\s*-\s*/g, " - ").replace(/\s+/g, " ").trim();
    // Remove leading/trailing dashes and clean up
    details = details.replace(/^-\s*/, "").replace(/\s*-$/, "").trim();
    
    return { mainHeading: mainHeading, details: details };
  }
  
  // If no pattern match, return original as main heading
  return { mainHeading: label, details: "" };
}

function updateMonoExponentialInfoPanel(monoExpParams) {
  const infoPanel = document.getElementById("regression-info-panel");
  const infoContent = document.getElementById("regression-info-content");
  
  if (!infoPanel || !infoContent) {
    return;
  }
  
  if (monoExpParams.length === 0) {
    infoContent.innerHTML = "<div style='color: #666;'>No fits available. Select wells and enable regression.</div>";
    return;
  }
  
  // Build HTML for parameters
  let html = "";
  
  // Model equation
  html += "<div style='margin-bottom: 8px; padding: 4px; background: white; border-radius: 2px;'>";
  html += "<strong>Model:</strong> y(t) = 1 + A(1 - e<sup>-k(t-t₀)</sup>)";
  html += "</div>";
  
  // Parameters for each fit
  monoExpParams.forEach((fit, idx) => {
    const p = fit.params;
    const borderColor = fit.color || '#2a6cff'; // Use dataset color or default blue
    
    // Parse label to extract main heading and details
    const parsed = parseLabelForLegend(fit.label);
    
    html += `<div style='margin-bottom: 12px; padding: 8px; background: white; border-radius: 2px; border-left: 3px solid ${borderColor};'>`;
    
    // Main heading in card color
    html += `<div style='color: ${borderColor}; font-weight: 600; font-size: 0.9rem; margin-bottom: 4px;'>${parsed.mainHeading}</div>`;
    
    // Details in gray below
    if (parsed.details) {
      html += `<div style='color: #666; font-size: 0.7rem; margin-bottom: 8px; line-height: 1.3;'>${parsed.details}</div>`;
    }
    
    // Parameters grid
    html += `<div style='display: grid; grid-template-columns: 1fr 1fr; gap: 4px; font-size: 0.7rem; margin-top: 8px;'>`;
    
    // Left column
    html += `<div><strong>A (amplitude):</strong> ${p.A.toFixed(4)}</div>`;
    html += `<div><strong>k (rate constant):</strong> ${p.k.toFixed(6)} s⁻¹</div>`;
    html += `<div><strong>τ (time constant):</strong> ${p.tau.toFixed(2)} s</div>`;
    
    // Right column
    html += `<div><strong>t₀ (start time):</strong> ${p.t0.toFixed(2)} s</div>`;
    html += `<div><strong>t₁/₂ (half-time):</strong> ${p.tHalf.toFixed(2)} s</div>`;
    html += `<div><strong>Plateau:</strong> ${p.plateau.toFixed(4)}</div>`;
    
    html += `</div></div>`;
  });
  
  infoContent.innerHTML = html;
}

function renderWellSelectorGrid() {
  const grid = document.getElementById("well-selector-grid");
  if (!grid || !selectedDatasetId) {
    return;
  }
  
  const selectedDatasets = getSelectedDatasets();
  if (selectedDatasets.length === 0) {
    grid.innerHTML = "<div style='color: #666; font-size: 0.8rem; padding: 8px;'>No datasets selected.</div>";
    return;
  }
  
  grid.innerHTML = "";
  
  // Collect all unique well IDs from selected datasets
  const allWellsSet = new Set();
  selectedDatasets.forEach(plateId => {
    const wells = allWellsData[plateId] || {};
    Object.keys(wells).forEach(wellId => {
      allWellsSet.add(wellId);
    });
  });
  
  // Top-left corner
  const corner = document.createElement("div");
  corner.className = "well-selector-header";
  grid.appendChild(corner);
  
  // Column headers
  for (let c = 1; c <= PLATE_COLS; c++) {
    const header = document.createElement("div");
    header.className = "well-selector-header";
    header.textContent = c;
    grid.appendChild(header);
  }
  
  // Rows
  PLATE_ROWS.forEach(rowLetter => {
    // Row label
    const rowLabel = document.createElement("div");
    rowLabel.className = "well-selector-header";
    rowLabel.textContent = rowLetter;
    grid.appendChild(rowLabel);
    
    // Wells
    for (let col = 1; col <= PLATE_COLS; col++) {
      const wellId = rowLetter + String(col);
      const cell = document.createElement("div");
      cell.className = "well-selector-cell";
      
      const hasData = allWellsSet.has(wellId);
      const isExcluded = enabledWells.size > 0 && !enabledWells.has(wellId);
      
      if (hasData) {
        cell.classList.add("has-data");
        if (isExcluded) {
          // Use inline styles with rgb to avoid hex escaping issues
          cell.style.background = "rgb(255, 204, 204)";
          cell.style.borderColor = "rgb(255, 107, 107)";
          cell.style.color = "rgb(204, 0, 0)";
          cell.style.fontWeight = "bold";
          cell.textContent = "✕";
          cell.title = `Well ${wellId} - EXCLUDED (click to include)`;
        } else {
          cell.textContent = col;
          cell.title = `Well ${wellId} - Click to exclude`;
        }
        
        cell.addEventListener("click", () => {
          toggleWellExclusion(wellId);
        });
      } else {
        cell.style.background = "#f2f2f2";
        cell.style.cursor = "default";
        cell.style.opacity = "0.5";
      }
      
      grid.appendChild(cell);
    }
  });
}

function toggleWellExclusion(wellId) {
  if (enabledWells.size === 0) {
    // If no exclusions yet, first collect all available wells
    const selectedDatasets = getSelectedDatasets();
    const allWellsSet = new Set();
    selectedDatasets.forEach(plateId => {
      const wells = allWellsData[plateId] || {};
      Object.keys(wells).forEach(wId => {
        allWellsSet.add(wId);
      });
    });
    // Enable all except the one being clicked
    allWellsSet.forEach(wId => {
      if (wId !== wellId) {
        enabledWells.add(wId);
      }
    });
  } else {
    // Toggle this specific well
    if (enabledWells.has(wellId)) {
      enabledWells.delete(wellId);
    } else {
      enabledWells.add(wellId);
    }
    
    // If all wells are enabled, clear the set (meaning all are enabled)
    const selectedDatasets = getSelectedDatasets();
    const allWellsSet = new Set();
    selectedDatasets.forEach(plateId => {
      const wells = allWellsData[plateId] || {};
      Object.keys(wells).forEach(wId => {
        allWellsSet.add(wId);
      });
    });
    
    // Check if all wells are enabled
    let allEnabled = true;
    allWellsSet.forEach(wId => {
      if (!enabledWells.has(wId)) {
        allEnabled = false;
      }
    });
    
    if (allEnabled) {
      enabledWells.clear(); // Clear means all enabled
    }
  }
  
  renderWellSelectorGrid();
  renderPlateGrid();
  updateWellInfo();
}

function getAllDatasetsInGroup() {
  // Get all datasets in the selected group (before filtering by enabledDatasets)
  if (!selectedDatasetId) {
    return [];
  }
  
  if (selectedDatasetId.startsWith("group:")) {
    const groupName = selectedDatasetId.substring(6);
    const selectedPlates = viewerData.plates.filter(p => p.column_group === groupName);
    return selectedPlates.map(p => p.id);
  } else {
    const plate = viewerData.plates.find(p => p.id === selectedDatasetId);
    return plate ? [selectedDatasetId] : [];
  }
}

function initializeTriplicateGroupColors() {
  // Initialize color map if needed (don't reset - keep consistent colors)
  // This ensures all wells with the same triplicate group name get the same color
  // For example, all wells in rows B, C, D with group name "250-10-2" get the same color
  if (Object.keys(triplicateGroupColorMap).length === 0 && viewerData && viewerData.config) {
    const triplicateGroups = viewerData.config.triplicate_groups || [];
    let colorIndex = 0;
    triplicateGroups.forEach(group => {
      const name = String(group.name || "").trim();
      if (name && !triplicateGroupColorMap[name]) {
        triplicateGroupColorMap[name] = TRIPLICATE_COLORS[colorIndex %% TRIPLICATE_COLORS.length];
        colorIndex++;
      }
    });
  }
}

function renderPlateGrid() {
  const grid = document.getElementById("plate-grid");
  grid.innerHTML = "";
  
  // Initialize color map if needed (don't reset - keep consistent colors)
  // This ensures all wells with the same triplicate group name get the same color
  if (Object.keys(triplicateGroupColorMap).length === 0) {
    initializeTriplicateGroupColors();
  }

  const selectedDatasets = getSelectedDatasets();
  const allGroupDatasets = getAllDatasetsInGroup();
  
  // Build column info map: column -> { genotype, buffer, scientist }
  // Use all datasets in group to show all available column groups
  const columnInfoMap = {}; // column -> { genotype, buffer, scientist }
  
  allGroupDatasets.forEach(plateId => {
    const plate = viewerData.plates.find(p => p.id === plateId);
    if (!plate || !plate.column_info) return;
    
    const columnInfo = plate.column_info;
    // Find which columns have data for this plate
    const wells = allWellsData[plateId] || {};
    Object.keys(wells).forEach(wellId => {
      const col = getColumnFromWellId(wellId);
      if (!columnInfoMap[col]) {
        columnInfoMap[col] = {
          genotype: columnInfo.genotype || "",
          buffer: columnInfo.buffer || "",
          scientist: columnInfo.scientist || ""
        };
      }
    });
  });

  // Top-left blank
  const corner = document.createElement("div");
  corner.className = "grid-header";
  grid.appendChild(corner);

  // Column number labels
  for (let c = 1; c <= PLATE_COLS; c++) {
    const div = document.createElement("div");
    div.className = "grid-header col-label";
    div.textContent = c;
    grid.appendChild(div);
  }

  const wellCounts = {}; // Track how many datasets have data for each well

  // Count datasets per well, filtering by enabled wells
  // Normalize well IDs to ensure consistent matching
  selectedDatasets.forEach(plateId => {
    const wells = allWellsData[plateId] || {};
    Object.keys(wells).forEach(wellId => {
      // Normalize well ID (remove spaces, ensure consistent format)
      const normalizedWellId = String(wellId).trim().toUpperCase();
      
      // Filter by enabledWells if any are set (check both original and normalized)
      if (enabledWells.size > 0 && !enabledWells.has(wellId) && !enabledWells.has(normalizedWellId)) {
        return;
      }
      
      // Filter by enabledColumns if any are set
      const column = getColumnFromWellId(wellId);
      if (enabledColumns.size > 0 && !enabledColumns.has(column)) {
        return;
      }
      
      if (!wellCounts[normalizedWellId]) {
        wellCounts[normalizedWellId] = new Set(); // Use Set to avoid duplicates
      }
      wellCounts[normalizedWellId].add(plateId);
    });
  });
  
  // Convert Sets to arrays for easier use
  Object.keys(wellCounts).forEach(wellId => {
    wellCounts[wellId] = Array.from(wellCounts[wellId]);
  });

  // Rows
  PLATE_ROWS.forEach(rowLetter => {
    // Row label
    const label = document.createElement("div");
    label.className = "grid-header row-label";
    label.textContent = rowLetter;
    grid.appendChild(label);

    // Wells
    for (let col = 1; col <= PLATE_COLS; col++) {
      const wellId = rowLetter + String(col);
      const normalizedWellId = wellId.toUpperCase(); // Normalize for matching
      const div = document.createElement("div");
      div.className = "well";
      const hasData = wellCounts.hasOwnProperty(normalizedWellId);
      div.dataset.wellId = wellId;
      div.dataset.hasData = hasData ? "true" : "false";
      
      if (hasData) {
        const plates = wellCounts[normalizedWellId];
        
        // Check if this well is a control well (configurable via control_rows)
        const isControl = isControlWell(wellId);
        
        // Get label from row-based config or content
        // Try to find the actual well ID in the data (might be different case/format)
        const firstPlateId = plates[0];
        const wells = allWellsData[firstPlateId] || {};
        // Try normalized ID first, then original ID, then search for matching well
        let actualWellId = normalizedWellId;
        if (!wells[actualWellId]) {
          actualWellId = wellId;
          if (!wells[actualWellId]) {
            // Search for matching well ID (case-insensitive)
            const matchingKey = Object.keys(wells).find(key => 
              key.toUpperCase() === normalizedWellId
            );
            if (matchingKey) {
              actualWellId = matchingKey;
            }
          }
        }
        const wellData = wells[actualWellId];
        if (!wellData) {
          // Well data not found, skip this well
          div.style.cursor = "default";
          div.style.background = "#f2f2f2";
          grid.appendChild(div);
          continue;
        }
        
        const config = viewerData.config || {};
        const wellLabels = config.well_labels || {};
        // Extract row letter from wellId (e.g., "B13" -> "B")
        const rowLetter = getRowLetter(wellId);
        
        // Count unique datasets that actually have data for this well
        const uniquePlates = plates.filter(plateId => {
          const plateWells = allWellsData[plateId] || {};
          // Check if this plate has data for this well (try different ID formats)
          return plateWells[actualWellId] || plateWells[wellId] || 
                 Object.keys(plateWells).some(key => key.toUpperCase() === normalizedWellId);
        });
        const datasetCount = uniquePlates.length;
        
        // For control wells: main text is "Control", label below is the control label
        // For non-control wells: main text is plate count or well ID, label below is the configured label
        let mainText, label;
        if (isControl) {
          mainText = "Control";
          // Get the control label from config or well data and format it
          label = wellLabels[rowLetter] || wellData.label || wellData.content || "";
          if (label) {
            label = formatLabel(label);
          }
        } else {
          // Show dataset count if more than 1, otherwise show well ID
          mainText = datasetCount > 1 ? `${datasetCount}` : wellId;
          // Use row-based label if available, otherwise use well's label or content
          label = wellLabels[rowLetter] || wellData.label || wellData.content || "";
          // Format label (remove uM, format as number-number-number)
          if (label) {
            label = formatLabel(label);
          }
        }
        
        div.textContent = mainText;
        div.style.cursor = "pointer";
        div.style.pointerEvents = "auto";
        div.style.position = "relative";
        div.style.zIndex = "1";
        
        // Add visual indicator for duplicate wells
        if (isDuplicateWell(wellId)) {
          div.style.borderWidth = "2px";
          div.style.borderStyle = "double";
          div.title = (div.title || "") + " ⚠️ Duplicate well - click to compare datasets";
        }
        
        if (label) {
          const labelDiv = document.createElement("div");
          labelDiv.className = "well-label";
          labelDiv.textContent = label;
          labelDiv.style.pointerEvents = "none";
          div.appendChild(labelDiv);
        }
        
        div.addEventListener("click", (e) => {
          e.stopPropagation();
          toggleWellSelection(wellId, plates);
        });
        if (isControl) {
          div.classList.add("control");
          div.title = "Control well (independently selectable)";
        }
        
        // Check if this well is part of a triplicate group (controls are never in triplicate groups)
        const column = getColumnFromWellId(wellId);
        const triplicateGroup = getTriplicateGroup(wellId, column);
        
        // Add visual indicator for triplicate groups with distinct colors
        // All wells in the same triplicate group (same name) get the same color
        // For example, all wells in rows B, C, D with group name "250-10-2" get the same color
        // Borders/strokes only appear when selected
        if (triplicateGroup && triplicateGroup.name) {
          div.classList.add("triplicate-group");
          div.title = `Triplicate group: ${triplicateGroup.name} (click to select all)`;
          
          // Apply distinct color for this triplicate group (background only, no border unless selected)
          // All wells with the same group name will get the same color
          // This ensures B, C, D rows all get the same color for the same triplicate group
          const groupName = String(triplicateGroup.name || "").trim(); // Normalize group name
          const color = getTriplicateGroupColor(groupName);
          if (color && color.bgHasData) {
            // Set background color only - borders will be added when selected
            const bgColor = hasData ? color.bgHasData : color.bg;
            div.style.setProperty("background-color", bgColor, "important");
          }
        }
        
        // Highlight if selected (check if any selected dataset has this well selected)
        // Also check if this well is part of a triplicate group that's selected
        const selectedDatasets = getSelectedDatasets();
        let isSelected = false;
        
        if (triplicateGroup) {
          // Check if any well in this triplicate group (same row pattern, same column) is selected
          const groupRowLetters = getRowLettersFromWells(triplicateGroup.wells);
          selectedDatasets.forEach(plateId => {
            groupRowLetters.forEach(rowLetter => {
              const wId = rowLetter + column;
              if (plates.includes(plateId)) {
                const key = `${plateId}:${wId}`;
                if (selectedWells.has(key)) {
                  isSelected = true;
                }
              }
            });
          });
        } else {
          // Check if this specific well is selected
          isSelected = selectedDatasets.some(plateId => {
            if (plates.includes(plateId)) {
              const key = `${plateId}:${wellId}`;
              return selectedWells.has(key);
            }
            return false;
          });
        }
        
        if (isSelected) {
          div.classList.add("selected");
          // Update background color and add border for selected triplicate groups
          if (triplicateGroup && triplicateGroup.name) {
            const groupName = String(triplicateGroup.name || "").trim(); // Normalize group name
            const color = getTriplicateGroupColor(groupName);
            if (color) {
              div.style.setProperty("background-color", color.bgHasData, "important");
              // Add border when selected
              div.style.setProperty("border-color", color.border, "important");
              div.style.setProperty("border-width", "2px", "important");
              div.style.setProperty("border-style", "solid", "important");
              div.style.setProperty("outline", `2px solid ${color.border}`, "important");
              div.style.setProperty("outline-offset", "-2px", "important");
            }
          }
        }
      } else {
        div.style.cursor = "default";
        div.style.background = "#f2f2f2";
      }
      grid.appendChild(div);
    }
  });
  
  // Update column legend with dataset mapping (use all datasets in group)
  updateColumnLegend(columnInfoMap, allGroupDatasets);
}

function updateColumnLegend(columnInfoMap, allDatasets) {
  const legendDiv = document.getElementById("column-legend-content");
  if (!legendDiv) return;
  
  // Build map: column -> { datasets: Set, info: { genotype, buffer, scientist } }
  const columnData = {}; // column -> { datasets: Set, info: {} }
  
  // Find which datasets have data in each column
  for (let col = 1; col <= PLATE_COLS; col++) {
    const info = columnInfoMap[col];
    if (!info) continue;
    
    columnData[col] = {
      datasets: new Set(),
      info: {
        genotype: info.genotype || "",
        buffer: info.buffer || "",
        scientist: info.scientist || ""
      }
    };
    
    // Find all datasets that have data in this column
    allDatasets.forEach(plateId => {
      const wells = allWellsData[plateId] || {};
      const hasDataInColumn = Object.keys(wells).some(wellId => {
        const wellCol = getColumnFromWellId(wellId);
        return wellCol === col;
      });
      
      if (hasDataInColumn) {
        columnData[col].datasets.add(plateId);
      }
    });
  }
  
  if (Object.keys(columnData).length === 0) {
    legendDiv.innerHTML = "<div style='color: #666; font-size: 0.8rem;'>No column information available.</div>";
    return;
  }
  
  // Build HTML with individual column checkboxes
  legendDiv.innerHTML = "";
  
  // Sort columns numerically
  const sortedColumns = Object.keys(columnData).map(Number).sort((a, b) => a - b);
  
  sortedColumns.forEach(col => {
    const colData = columnData[col];
    const colDatasets = Array.from(colData.datasets);
    
    if (colDatasets.length === 0) return; // Skip columns with no data
    
    // Check if this column is enabled (empty set = all enabled, or column is in the set)
    const columnEnabled = enabledColumns.size === 0 || enabledColumns.has(col);
    
    const div = document.createElement("div");
    div.style.marginBottom = "6px";
    div.style.padding = "6px";
    div.style.background = "#f9f9f9";
    div.style.borderRadius = "3px";
    
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.id = `column-${col}`;
    checkbox.checked = columnEnabled;
    checkbox.style.marginRight = "8px";
    checkbox.style.cursor = "pointer";
    checkbox.dataset.column = col;
    checkbox.dataset.datasets = JSON.stringify(colDatasets);
    
    checkbox.addEventListener("change", (e) => {
      const column = parseInt(e.target.dataset.column);
      const allGroupDatasets = getAllDatasetsInGroup();
      
      // Get all columns that have data
      const allColumnsWithData = new Set();
      allGroupDatasets.forEach(plateId => {
        const wells = allWellsData[plateId] || {};
        Object.keys(wells).forEach(wellId => {
          const col = getColumnFromWellId(wellId);
          allColumnsWithData.add(col);
        });
      });
      
      if (e.target.checked) {
        // Enabling: add column to enabled set
        enabledColumns.add(column);
        
        // If all columns are now enabled, clear the set (meaning all enabled)
        if (Array.from(allColumnsWithData).every(col => enabledColumns.has(col))) {
          enabledColumns.clear();
        }
      } else {
        // Disabling: if set is empty (all enabled), populate it with all columns except the one being unchecked
        if (enabledColumns.size === 0) {
          // Currently all are enabled, so enable all columns except the one being unchecked
          allColumnsWithData.forEach(col => {
            if (col !== column) {
              enabledColumns.add(col);
            }
          });
        } else {
          // Some columns are already disabled, remove this one from enabled set
          enabledColumns.delete(column);
        }
      }
      
      renderPlateGrid();
      updateWellInfo();
      updateChart();
    });
    
    const label = document.createElement("label");
    label.htmlFor = `column-${col}`;
    label.style.cursor = "pointer";
    label.style.display = "flex";
    label.style.alignItems = "center";
    label.style.fontWeight = "600";
    
    label.appendChild(checkbox);
    const labelText = document.createElement("span");
    labelText.textContent = `Column ${col}`;
    label.appendChild(labelText);
    
    div.appendChild(label);
    
    const details = document.createElement("div");
    details.style.fontSize = "0.75rem";
    details.style.marginLeft = "24px";
    details.style.marginTop = "4px";
    
    if (colData.info.genotype) {
      const genDiv = document.createElement("div");
      genDiv.style.marginBottom = "2px";
      genDiv.innerHTML = `<strong>Genotype:</strong> ${colData.info.genotype}`;
      details.appendChild(genDiv);
    }
    if (colData.info.buffer) {
      const bufDiv = document.createElement("div");
      bufDiv.style.marginBottom = "2px";
      bufDiv.innerHTML = `<strong>Buffer:</strong> ${colData.info.buffer}`;
      details.appendChild(bufDiv);
    }
    if (colData.info.scientist) {
      const sciDiv = document.createElement("div");
      sciDiv.innerHTML = `<strong>Scientist:</strong> ${colData.info.scientist}`;
      details.appendChild(sciDiv);
    }
    
    if (colDatasets.length > 0) {
      const countDiv = document.createElement("div");
      countDiv.style.marginTop = "4px";
      countDiv.style.fontSize = "0.7rem";
      countDiv.style.color = "#666";
      countDiv.textContent = `${colDatasets.length} dataset${colDatasets.length !== 1 ? "s" : ""}`;
      details.appendChild(countDiv);
    }
    
    div.appendChild(details);
    legendDiv.appendChild(div);
  });
}

function toggleWellSelection(wellId, plateIds) {
  // Toggle selection for all selected datasets
  if (!selectedDatasetId) return;
  
  const selectedDatasets = getSelectedDatasets();
  const column = getColumnFromWellId(wellId);
  const isControl = isControlWell(wellId);
  
  // Control wells: allow individual selection (can select each well separately)
  if (isControl) {
    // Check if this specific control well is already selected
    let isSelected = false;
    selectedDatasets.forEach(plateId => {
      if (plateIds.includes(plateId)) {
        const key = `${plateId}:${wellId}`;
        if (selectedWells.has(key)) {
          isSelected = true;
        }
      }
    });
    
    // Toggle just this individual control well across all selected datasets
    let anyChanged = false;
    selectedDatasets.forEach(plateId => {
      // Only toggle if this plate has data for this well
      if (plateIds.includes(plateId)) {
        // Filter by enabledWells if any are set
        if (enabledWells.size > 0 && !enabledWells.has(wellId)) {
          return; // Skip excluded wells
        }
        const key = `${plateId}:${wellId}`;
        if (isSelected) {
          selectedWells.delete(key);
        } else {
          selectedWells.add(key);
        }
        anyChanged = true;
      }
    });
    
    if (anyChanged) {
      renderPlateGrid();
      updateWellInfo();
      updateChart();
    }
    return;
  }
  
  // For duplicate wells, automatically select all datasets for comparison
  if (isDuplicateWell(wellId)) {
    // Find all datasets that have this well
    const allGroupDatasets = getAllDatasetsInGroup();
    const datasetsWithWell = allGroupDatasets.filter(plateId => {
      return allWellsData[plateId] && allWellsData[plateId][wellId];
    });
    
    // Check if any dataset for this well is already selected
    let isAnySelected = datasetsWithWell.some(plateId => {
      const key = `${plateId}:${wellId}`;
      return selectedWells.has(key);
    });
    
    // Toggle all datasets for this duplicate well
    let anyChanged = false;
    datasetsWithWell.forEach(plateId => {
      const key = `${plateId}:${wellId}`;
      if (isAnySelected) {
        selectedWells.delete(key);
      } else {
        selectedWells.add(key);
      }
      anyChanged = true;
    });
    
    if (anyChanged) {
      renderPlateGrid();
      updateWellInfo();
      updateChart();
    }
    return;
  }
  
  // For non-control wells, check triplicate groups
  const triplicateGroup = getTriplicateGroup(wellId, column);
  
  let anyChanged = false;
  
  // Check if triplicate group is selected (check if ALL wells in the group are selected)
  let isGroupSelected = false;
  if (triplicateGroup) {
    const groupRowLetters = getRowLettersFromWells(triplicateGroup.wells);
    // Count how many wells in the group are selected
    let selectedCount = 0;
    let totalCount = 0;
    
    selectedDatasets.forEach(plateId => {
      groupRowLetters.forEach(rowLetter => {
        const wId = rowLetter + column;
        const key = `${plateId}:${wId}`;
        const hasData = allWellsData[plateId] && allWellsData[plateId][wId];
        if (hasData) {
          totalCount++;
          if (selectedWells.has(key)) {
            selectedCount++;
          }
        }
      });
    });
    
    // Consider group selected if at least one well is selected (for toggle behavior)
    // This allows clicking to toggle the whole group
    isGroupSelected = selectedCount > 0;
  } else {
    // For non-triplicate wells, check if this specific well is selected
    selectedDatasets.forEach(plateId => {
      if (plateIds.includes(plateId)) {
        const key = `${plateId}:${wellId}`;
        if (selectedWells.has(key)) {
          isGroupSelected = true;
        }
      }
    });
  }
  
  // Toggle selection: if group/well is selected, remove all; otherwise add all
  if (triplicateGroup) {
    // Get row letters from the triplicate group pattern
    const groupRowLetters = getRowLettersFromWells(triplicateGroup.wells);
    
    // Toggle all wells in the triplicate group (same row pattern, same column) across all selected datasets
    selectedDatasets.forEach(plateId => {
      groupRowLetters.forEach(rowLetter => {
        const wId = rowLetter + column;
        const key = `${plateId}:${wId}`;
        const hasData = allWellsData[plateId] && allWellsData[plateId][wId];
        if (hasData) {
          // Filter by enabledWells if any are set
          if (enabledWells.size > 0 && !enabledWells.has(wId)) {
            return; // Skip excluded wells
          }
          if (isGroupSelected) {
            selectedWells.delete(key);
          } else {
            selectedWells.add(key);
          }
          anyChanged = true;
        }
      });
    });
  } else {
    // Toggle individual well across all selected datasets
    selectedDatasets.forEach(plateId => {
      // Only toggle if this plate has data for this well
      if (plateIds.includes(plateId)) {
        // Filter by enabledWells if any are set
        if (enabledWells.size > 0 && !enabledWells.has(wellId)) {
          return; // Skip excluded wells
        }
        const key = `${plateId}:${wellId}`;
        if (isGroupSelected) {
          selectedWells.delete(key);
        } else {
          selectedWells.add(key);
        }
        anyChanged = true;
      }
    });
  }
  
  if (anyChanged) {
    renderPlateGrid();
    updateWellInfo();
  }
}

function clearSelection() {
  selectedWells.clear();
  renderPlateGrid();
  updateChart();
}

function updateWellInfo() {
  const el = document.getElementById("well-info");
  if (selectedWells.size === 0) {
    const selectedDatasets = getSelectedDatasets();
    const datasetCount = selectedDatasets.length;
    if (datasetCount > 1) {
      el.textContent = `Selected group with ${datasetCount} datasets. Click wells on the plate to add them to the chart.`;
    } else {
      el.textContent = "Click wells on the plate to add them to the chart.";
    }
    return;
  }
  const selectedDatasets = getSelectedDatasets();
  const datasetCount = selectedDatasets.length;
  const uniqueWells = new Set(Array.from(selectedWells).map(key => key.split(":")[1]));
  const wellCount = uniqueWells.size;
  if (datasetCount > 1) {
    el.textContent = `${wellCount} well(s) selected across ${datasetCount} datasets. Click "Update Chart" to display.`;
  } else {
    el.textContent = `${wellCount} well(s) selected. Click "Update Chart" to display.`;
  }
}

function isControlWell(wellId) {
  // Control wells are configurable via control_rows in config (default: ["N", "O"])
  const config = viewerData.config || {};
  const controlRows = config.control_rows || ["N", "O"];
  const rowLetter = wellId.replace(/[0-9]/g, '');
  return controlRows.includes(rowLetter);
}

function formatLabel(label) {
  // Format labels: remove "uM" and convert to "number-number-number" format
  // Examples: "250 uM 10-2" -> "250-10-2", "500 uM 5-2" -> "500-5-2"
  if (!label) return label;
  
  // Remove "uM" (case insensitive) and extra spaces
  let formatted = label.replace(/\s*uM\s*/gi, ' ').replace(/\s+/g, ' ').trim();
  
  // Extract numbers and hyphens, replace spaces with hyphens
  // Match pattern: number(s) followed by optional space/hyphen and more numbers
  formatted = formatted.replace(/(\d+)\s*-\s*(\d+)/g, '$1-$2'); // Handle existing hyphens
  formatted = formatted.replace(/(\d+)\s+(\d+)/g, '$1-$2'); // Replace spaces between numbers with hyphens
  
  return formatted;
}

function getRowLetter(wellId) {
  // Extract row letter from well ID (e.g., "B13" -> "B")
  return wellId.replace(/[0-9]/g, '');
}

function getTriplicateGroup(wellId, column) {
  // Control wells should never be part of triplicate groups
  if (isControlWell(wellId)) {
    return null;
  }
  
  const config = viewerData.config || {};
  const triplicateGroups = config.triplicate_groups || [];
  const rowLetter = getRowLetter(wellId);
  
  // Find group that matches the row pattern (BCD, EFG, etc.) - applies to all columns
  return triplicateGroups.find(group => {
    if (!group.wells || group.wells.length === 0) {
      return false;
    }
    
    // Extract row letters from the group's wells
    const groupRowLetters = getRowLettersFromWells(group.wells);
    
    // Check if this well's row letter matches any row in the group pattern
    // Groups apply to all columns, so no column check needed
    return groupRowLetters.includes(rowLetter);
  });
}

function getColumnFromWellId(wellId) {
  // Extract column number from well ID (e.g., "B13" -> 13)
  return parseInt(wellId.replace(/[A-Z]/g, ''));
}

function getRowLettersFromWells(wells) {
  // Extract row letters from wells array
  // Wells can be either row letters (e.g., ["B", "C", "D"]) or full well IDs (e.g., ["B13", "C13", "D13"])
  return wells.map(w => {
    // If it's a full well ID, extract row letter; otherwise assume it's already a row letter
    if (w.length > 1 && /^[A-Z]/.test(w) && /[0-9]/.test(w)) {
      return getRowLetter(w);
    }
    return w; // Already a row letter
  });
}

// Color palette for triplicate groups - each group gets a distinct color
const TRIPLICATE_COLORS = [
  { border: "#ff6b6b", bg: "#ffe0e0", bgSelected: "#ffcccc", bgHasData: "#ffe8e8" }, // Red
  { border: "#4ecdc4", bg: "#e0f7f5", bgSelected: "#ccf0ed", bgHasData: "#e8f9f7" }, // Teal
  { border: "#45b7d1", bg: "#e0f2f7", bgSelected: "#cce8f0", bgHasData: "#e8f5fa" }, // Blue
  { border: "#f9ca24", bg: "#fff9e0", bgSelected: "#fff3cc", bgHasData: "#fffbe8" }, // Yellow
  { border: "#6c5ce7", bg: "#edeaf7", bgSelected: "#ddd4f0", bgHasData: "#f0ecf9" }, // Purple
  { border: "#a29bfe", bg: "#edeaf7", bgSelected: "#ddd4f0", bgHasData: "#f0ecf9" }, // Light Purple
  { border: "#fd79a8", bg: "#ffe0ed", bgSelected: "#ffcce0", bgHasData: "#ffe8f0" }, // Pink
  { border: "#00b894", bg: "#e0f5f0", bgSelected: "#ccede5", bgHasData: "#e8f7f2" }, // Green
  { border: "#e17055", bg: "#ffe0d9", bgSelected: "#ffccbb", bgHasData: "#ffe8e0" }, // Orange
  { border: "#74b9ff", bg: "#e0f0ff", bgSelected: "#cce0ff", bgHasData: "#e8f5ff" }, // Light Blue
  { border: "#55efc4", bg: "#e0faf5", bgSelected: "#ccf5eb", bgHasData: "#e8fbf7" }, // Mint
  { border: "#fdcb6e", bg: "#fff5e0", bgSelected: "#ffebcc", bgHasData: "#fff8e8" }  // Light Orange
];

// Map triplicate group names to colors (consistent across the app)
let triplicateGroupColorMap = {};

function getTriplicateGroupColor(groupName) {
  if (!groupName) return null;
  
  // Ensure color map is initialized
  initializeTriplicateGroupColors();
  
  // Normalize group name to ensure consistent lookup
  const normalizedName = String(groupName).trim();
  
  // Return the color for this group name (same color for all wells in the group)
  // All wells with the same group name will get the same color
  // For example, all wells in rows B, C, D with group "250-10-2" get the same color
  return triplicateGroupColorMap[normalizedName] || TRIPLICATE_COLORS[0];
}

function isDuplicateWell(wellId) {
  // Check if this well ID is in the duplicate warnings
  if (!viewerData || !viewerData.config || !viewerData.config.warnings) {
    return false;
  }
  const warnings = viewerData.config.warnings.duplicate_columns || [];
  return warnings.some(warning => {
    const wellIds = warning.well_ids || [];
    return wellIds.includes(wellId);
  });
}

function getDuplicateFiles(wellId) {
  // Get the list of files that have this duplicate well
  if (!viewerData || !viewerData.config || !viewerData.config.warnings) {
    return [];
  }
  const warnings = viewerData.config.warnings.duplicate_columns || [];
  for (const warning of warnings) {
    const wellIds = warning.well_ids || [];
    if (wellIds.includes(wellId)) {
      return warning.files || [];
    }
  }
  return [];
}

function getAllUniqueTimePoints(wells) {
  // Collect all unique time points from all wells
  const allTimePoints = new Set();
  wells.forEach(well => {
    if (well.time_s) {
      well.time_s.forEach(t => {
        if (t !== null && t !== undefined) {
          allTimePoints.add(t);
        }
      });
    }
  });
  return Array.from(allTimePoints).sort((a, b) => a - b);
}

function calculateMeanAndError(wells, timePoints) {
  // Calculate mean and standard error for triplicates at each time point
  const means = [];
  const errors = [];
  
  timePoints.forEach((time, idx) => {
    const values = wells.map(w => {
      const timeIdx = w.time_s.indexOf(time);
      return timeIdx >= 0 ? w.values[timeIdx] : null;
    }).filter(v => v !== null);
    
    if (values.length === 0) {
      means.push(null);
      errors.push(null);
      return;
    }
    
    const mean = values.reduce((a, b) => a + b, 0) / values.length;
    // Use sample variance (divide by n-1) for proper SEM calculation
    const n = values.length;
    const variance = n > 1 
      ? values.reduce((sum, val) => sum + Math.pow(val - mean, 2), 0) / (n - 1)
      : 0;
    const stdDev = Math.sqrt(variance);
    const stdError = stdDev / Math.sqrt(n);
    
    means.push(mean);
    errors.push(stdError);
  });
  
  return { means, errors };
}

function processWellsAndCalculateMean(wells, isNormalized = false, baselineEndTime = 24, cutoffTime = 360) {
  // WORKFLOW: Always process each well individually first, then average
  // If normalized: normalize each well to its OWN baseline first, then calculate mean and error
  // If not normalized: process each well individually (no normalization), then calculate mean and error
  // This ensures correct normalization even when combining wells from different columns with different baselines
  // Returns { filteredTimePoints, means, errors, originalMeans }
  
  // Step 1: Process each well individually (normalize if enabled, otherwise use raw data)
  const processedWells = wells.map(well => {
    const filteredTimePoints = [];
    const processedValues = [];
    const originalValues = [];
    
    if (!well.time_s || !well.values) {
      return { timePoints: [], values: [], originalValues: [] };
    }
    
    if (isNormalized) {
      // Normalize this well by its own baseline
      const baselineValues = [];
      for (let i = 0; i < well.time_s.length; i++) {
        const time = well.time_s[i];
        const val = well.values[i];
        if (time !== null && val !== null && time <= baselineEndTime) {
          baselineValues.push(val);
        }
      }
      
      if (baselineValues.length === 0) {
        // No baseline data for this well, skip it
        return { timePoints: [], values: [], originalValues: [] };
      }
      
      const wellBaseline = baselineValues.reduce((a, b) => a + b, 0) / baselineValues.length;
      
      if (wellBaseline === 0 || !isFinite(wellBaseline)) {
        // Invalid baseline for this well, skip it
        return { timePoints: [], values: [], originalValues: [] };
      }
      
      // Normalize and filter: only include points after baselineEndTime and up to cutoffTime
      for (let i = 0; i < well.time_s.length; i++) {
        const time = well.time_s[i];
        const val = well.values[i];
        
        if (time === null || val === null) continue;
        
        // Only include points after baselineEndTime and up to cutoffTime
        if (time > baselineEndTime && time <= cutoffTime) {
          // Subtract baselineEndTime to bring first point to x=0 (time after injection)
          filteredTimePoints.push(time - baselineEndTime);
          // Divide by this well's own baseline
          processedValues.push(val / wellBaseline);
          originalValues.push(val); // Store original value
        }
      }
    } else {
      // Not normalized: use all raw data points
      for (let i = 0; i < well.time_s.length; i++) {
        const time = well.time_s[i];
        const val = well.values[i];
        
        if (time === null || val === null) continue;
        
        filteredTimePoints.push(time);
        processedValues.push(val);
        originalValues.push(val); // Store original value (same as processed when not normalized)
      }
    }
    
    return { timePoints: filteredTimePoints, values: processedValues, originalValues: originalValues };
  }).filter(well => well.timePoints.length > 0); // Filter out wells with no data
  
  // Step 2: Find all unique time points across all processed wells
  const allTimePoints = new Set();
  processedWells.forEach(well => {
    well.timePoints.forEach(t => {
      if (t !== null) allTimePoints.add(t);
    });
  });
  
  const sortedTimePoints = Array.from(allTimePoints).sort((a, b) => a - b);
  
  // Step 3: Calculate mean and error at each time point from processed values
  // Also track original values for tooltip display
  const means = [];
  const errors = [];
  const originalMeans = [];
  
  sortedTimePoints.forEach(time => {
    const values = [];
    const origValues = [];
    processedWells.forEach(well => {
      const timeIdx = well.timePoints.indexOf(time);
      if (timeIdx >= 0 && well.values[timeIdx] !== null) {
        values.push(well.values[timeIdx]);
        if (well.originalValues && well.originalValues[timeIdx] !== undefined) {
          origValues.push(well.originalValues[timeIdx]);
        }
      }
    });
    
    if (values.length === 0) {
      means.push(null);
      errors.push(null);
      originalMeans.push(null);
      return;
    }
    
    const mean = values.reduce((a, b) => a + b, 0) / values.length;
    // Use sample variance (divide by n-1) for proper SEM calculation
    const n = values.length;
    const variance = n > 1 
      ? values.reduce((sum, val) => sum + Math.pow(val - mean, 2), 0) / (n - 1)
      : 0;
    const stdDev = Math.sqrt(variance);
    const stdError = stdDev / Math.sqrt(n);
    
    means.push(mean);
    errors.push(stdError);
    
    // Calculate mean of original values
    if (origValues.length > 0) {
      const origMean = origValues.reduce((a, b) => a + b, 0) / origValues.length;
      originalMeans.push(origMean);
    } else {
      originalMeans.push(null);
    }
  });
  
  return { filteredTimePoints: sortedTimePoints, normalizedMeans: means, normalizedErrors: errors, originalMeans: originalMeans };
}

// Keep old function name for backwards compatibility, but it now uses the unified workflow
function calculateMeanAndErrorFromNormalized(wells, baselineEndTime = 24, cutoffTime = 360) {
  return processWellsAndCalculateMean(wells, true, baselineEndTime, cutoffTime);
}

function normalizeData(timePoints, values, baselineEndTime = 24, cutoffTime = 360) {
  // Normalize data: calculate baseline from 0-baselineEndTime, then plot only points after baselineEndTime
  // divided by baseline, up to cutoffTime. Baseline period is NOT included in the plot.
  // Returns { filteredTimePoints, normalizedValues, originalValues } - filtered arrays
  if (!timePoints || !values || timePoints.length !== values.length) {
    // If invalid input, return empty filtered arrays
    return { filteredTimePoints: [], normalizedValues: [], originalValues: [] };
  }
  
  // Find baseline: average of all values where time <= baselineEndTime
  const baselineValues = [];
  for (let i = 0; i < timePoints.length; i++) {
    if (timePoints[i] !== null && values[i] !== null && timePoints[i] <= baselineEndTime) {
      baselineValues.push(values[i]);
    }
  }
  
  if (baselineValues.length === 0) {
    // No baseline data found, return empty filtered arrays (can't normalize without baseline)
    return { filteredTimePoints: [], normalizedValues: [], originalValues: [] };
  }
  
  const baseline = baselineValues.reduce((a, b) => a + b, 0) / baselineValues.length;
  
  if (baseline === 0 || !isFinite(baseline)) {
    // Avoid division by zero or invalid baseline
    return { filteredTimePoints: [], normalizedValues: [], originalValues: [] };
  }
  
  // Filter and normalize: only include points where baselineEndTime < time <= cutoffTime
  // EVERY point after baselineEndTime is divided by the SAME baseline average
  const filteredTimePoints = [];
  const normalizedValues = [];
  const originalValues = [];
  
  for (let i = 0; i < timePoints.length; i++) {
    const time = timePoints[i];
    const val = values[i];
    
    if (time === null || val === null) {
      continue;
    }
    
    // Only include points after baselineEndTime and up to cutoffTime (baseline NOT included)
    if (time > baselineEndTime && time <= cutoffTime) {
      // Subtract baselineEndTime to bring first point to x=0 (time after injection)
      filteredTimePoints.push(time - baselineEndTime);
      // Divide EVERY time point by the SAME baseline average
      normalizedValues.push(val / baseline);
      originalValues.push(val); // Store original value for tooltip
    }
  }
  
  return { filteredTimePoints, normalizedValues, originalValues };
}

// Fit baseline function using specified method
function fitBaselineFunction(timePoints, values, method = 'lowess', frac = 0.5, polyOrder = 1) {
  // Extract baseline window (time <= baselineEndTime will be handled by caller)
  // This function fits on the provided points and returns a function g(t) that can be evaluated at any time
  
  if (!timePoints || !values || timePoints.length === 0) {
    return null;
  }
  
  // Filter out NaN values
  const validData = [];
  for (let i = 0; i < timePoints.length; i++) {
    if (timePoints[i] !== null && values[i] !== null && !isNaN(timePoints[i]) && !isNaN(values[i])) {
      validData.push({ t: timePoints[i], v: values[i] });
    }
  }
  
  if (validData.length === 0) {
    return null;
  }
  
  validData.sort((a, b) => a.t - b.t);
  const times = validData.map(d => d.t);
  const vals = validData.map(d => d.v);
  
  if (method === 'constant') {
    // Constant baseline: g(t) = mean(F_baseline)
    const baselineValue = vals.reduce((a, b) => a + b, 0) / vals.length;
    return (t) => baselineValue;
  } else if (method === 'polynomial') {
    // Polynomial fit: g(t) = c₀ + c₁t + c₂t² + ...
    // Simple polynomial regression
    const n = times.length;
    const X = [];
    const y = vals;
    
    // Build design matrix
    for (let i = 0; i < n; i++) {
      const row = [];
      for (let p = polyOrder; p >= 0; p--) {
        row.push(Math.pow(times[i], p));
      }
      X.push(row);
    }
    
    // Solve using normal equations: (X'X)^(-1)X'y
    // Simplified for small orders
    let coeffs;
    if (polyOrder === 1) {
      // Linear: y = a + b*t
      const sumT = times.reduce((a, b) => a + b, 0);
      const sumT2 = times.reduce((a, b) => a + b * b, 0);
      const sumY = y.reduce((a, b) => a + b, 0);
      const sumTY = times.reduce((sum, t, i) => sum + t * y[i], 0);
      const n = times.length;
      
      const det = n * sumT2 - sumT * sumT;
      if (Math.abs(det) < 1e-10) {
        // Fallback to constant
        const mean = sumY / n;
        return (t) => mean;
      }
      
      const a = (sumY * sumT2 - sumT * sumTY) / det;
      const b = (n * sumTY - sumT * sumY) / det;
      coeffs = [b, a]; // [c1, c0] for poly1d format
    } else {
      // Higher order - use simple least squares
      // For simplicity, use a basic implementation
      const meanT = sumT / n;
      const meanY = sumY / n;
      // Fallback to constant for now if order > 1
      const mean = meanY;
      return (t) => mean;
    }
    
    return (t) => {
      let result = 0;
      for (let i = 0; i < coeffs.length; i++) {
        result += coeffs[i] * Math.pow(t, coeffs.length - 1 - i);
      }
      return result;
    };
  } else if (method === 'lowess') {
    // LOWESS smoothing - simplified implementation
    // For a proper LOWESS, we'd need a library, but we can approximate with moving average
    // This is a simplified version - for production, consider using a library like regression.js
    
    // Use a simple moving average with weighted window as approximation
    const windowSize = Math.max(1, Math.floor(frac * times.length));
    const smoothed = [];
    
    for (let i = 0; i < times.length; i++) {
      let sum = 0;
      let weightSum = 0;
      const center = times[i];
      
      for (let j = 0; j < times.length; j++) {
        const dist = Math.abs(times[j] - center);
        const maxDist = Math.max(...times.map(t => Math.abs(t - center)));
        const weight = dist < maxDist * frac ? Math.pow(1 - (dist / (maxDist * frac)), 3) : 0; // Tricube weight
        
        sum += vals[j] * weight;
        weightSum += weight;
      }
      
      smoothed.push(weightSum > 0 ? sum / weightSum : vals[i]);
    }
    
    // Return interpolation function
    return (t) => {
      if (t <= times[0]) return smoothed[0];
      if (t >= times[times.length - 1]) return smoothed[smoothed.length - 1];
      
      // Linear interpolation
      for (let i = 0; i < times.length - 1; i++) {
        if (t >= times[i] && t <= times[i + 1]) {
          const ratio = (t - times[i]) / (times[i + 1] - times[i]);
          return smoothed[i] * (1 - ratio) + smoothed[i + 1] * ratio;
        }
      }
      return smoothed[smoothed.length - 1];
    };
  }
  
  return null;
}

// Normalize well data using fitted baseline (Step 4)
function normalizeWithFittedBaseline(wellData, baselineEndTime, method, frac, polyOrder, normalizationMode) {
  if (!wellData || !wellData.time_s || !wellData.values) {
    return wellData;
  }
  
  // Extract baseline window
  const baselineTimes = [];
  const baselineValues = [];
  for (let i = 0; i < wellData.time_s.length; i++) {
    const time = wellData.time_s[i];
    const val = wellData.values[i];
    if (time !== null && val !== null && time <= baselineEndTime) {
      baselineTimes.push(time);
      baselineValues.push(val);
    }
  }
  
  if (baselineTimes.length === 0) {
    return wellData; // No baseline data
  }
  
  // Fit baseline function
  const baselineFunc = fitBaselineFunction(baselineTimes, baselineValues, method, frac, polyOrder);
  if (!baselineFunc) {
    return wellData;
  }
  
  // Normalize all timepoints using fitted baseline
  const normalizedTimes = [];
  const normalizedValues = [];
  const originalValues = [];
  
  for (let i = 0; i < wellData.time_s.length; i++) {
    const time = wellData.time_s[i];
    const value = wellData.values[i];
    
    if (time === null || value === null) continue;
    
    const baselineValue = baselineFunc(time);
    if (baselineValue === 0 || !isFinite(baselineValue)) continue;
    
    let normValue;
    if (normalizationMode === 'delta_f_over_f') {
      // ΔF/F = (F(t) - g(t)) / g(t)
      normValue = (value - baselineValue) / baselineValue;
    } else {
      // Multiplicative: F_norm(t) = F(t) / g(t)
      normValue = value / baselineValue;
    }
    
    normalizedTimes.push(time);
    normalizedValues.push(normValue);
    originalValues.push(value);
  }
  
  if (normalizedTimes.length === 0) {
    return null;
  }
  
  return {
    time_s: normalizedTimes,
    values: normalizedValues,
    originalValues: originalValues,
    label: wellData.label,
    content: wellData.content
  };
}

function normalizeWellData(wellData, isNormalized, baselineEndTime = 24, cutoffTime = 360) {
  // STEP 3: Normalize a single well's data object (if normalization is enabled)
  // This happens AFTER baseline alignment and grouping into triplicates
  // Returns a new well data object with normalized time_s and values (or original if not normalized)
  if (!wellData || !wellData.time_s || !wellData.values) {
    return wellData; // Return as-is if invalid
  }
  
  if (!isNormalized) {
    // Not normalizing - return original data
    return {
      time_s: wellData.time_s,
      values: wellData.values,
      label: wellData.label,
      content: wellData.content,
      originalValues: wellData.originalValues || wellData.values
    };
  }
  
  // Normalize this well
  const normalized = normalizeData(wellData.time_s, wellData.values, baselineEndTime, cutoffTime);
  
  if (normalized.filteredTimePoints.length === 0) {
    // No data after normalization - return null to indicate this well should be skipped
    return null;
  }
  
  // Return normalized well data
  return {
    time_s: normalized.filteredTimePoints,
    values: normalized.normalizedValues,
    originalValues: normalized.originalValues, // Store for tooltip display
    label: wellData.label,
    content: wellData.content
  };
}


// Fit regression curve to data points (Step 4)
function fitRegressionCurve(dataPoints, regressionType = 'mono-exponential', polynomialOrder = 2) {
  // dataPoints is an array of {x, y} objects
  if (!dataPoints || dataPoints.length === 0) {
    return null;
  }
  
  // Filter valid points
  const validPoints = dataPoints.filter(p => p.x !== null && p.y !== null && isFinite(p.x) && isFinite(p.y));
  if (validPoints.length < 2) {
    return null;
  }
  
  // Sort by x
  validPoints.sort((a, b) => a.x - b.x);
  const x = validPoints.map(p => p.x);
  const y = validPoints.map(p => p.y);
  
  let fitFunction;
  
  if (regressionType === 'linear') {
    // Linear: y = a + b*x
    const n = x.length;
    const sumX = x.reduce((a, b) => a + b, 0);
    const sumY = y.reduce((a, b) => a + b, 0);
    const sumXY = x.reduce((sum, xi, i) => sum + xi * y[i], 0);
    const sumX2 = x.reduce((sum, xi) => sum + xi * xi, 0);
    
    const det = n * sumX2 - sumX * sumX;
    if (Math.abs(det) < 1e-10) {
      // Fallback to constant
      const mean = sumY / n;
      fitFunction = (t) => mean;
    } else {
      const a = (sumY * sumX2 - sumX * sumXY) / det;
      const b = (n * sumXY - sumX * sumY) / det;
      fitFunction = (t) => a + b * t;
    }
  } else if (regressionType === 'polynomial') {
    // Polynomial: y = c0 + c1*x + c2*x^2 + ... + cn*x^n
    const n = x.length;
    const order = Math.min(polynomialOrder, n - 1); // Can't fit order >= n points
    
    // Build design matrix
    const X = [];
    for (let i = 0; i < n; i++) {
      const row = [];
      for (let p = order; p >= 0; p--) {
        row.push(Math.pow(x[i], p));
      }
      X.push(row);
    }
    
    // Solve using normal equations: (X'X)^(-1)X'y
    // Simplified for small orders
    let coeffs;
    if (order === 1) {
      // Linear case
      const sumX = x.reduce((a, b) => a + b, 0);
      const sumX2 = x.reduce((sum, xi) => sum + xi * xi, 0);
      const sumY = y.reduce((a, b) => a + b, 0);
      const sumXY = x.reduce((sum, xi, i) => sum + xi * y[i], 0);
      const det = n * sumX2 - sumX * sumX;
      if (Math.abs(det) < 1e-10) {
        coeffs = [sumY / n, 0];
      } else {
        const a = (sumY * sumX2 - sumX * sumXY) / det;
        const b = (n * sumXY - sumX * sumY) / det;
        coeffs = [a, b]; // [c0, c1]
      }
    } else {
      // Higher order - use simple least squares approximation
      // For simplicity, use a basic implementation
      const meanY = y.reduce((a, b) => a + b, 0) / n;
      coeffs = [meanY];
      for (let p = 1; p <= order; p++) {
        coeffs.push(0);
      }
    }
    
    fitFunction = (t) => {
      let result = 0;
      for (let i = 0; i < coeffs.length; i++) {
        result += coeffs[i] * Math.pow(t, coeffs.length - 1 - i);
      }
      return result;
    };
  } else if (regressionType === 'exponential') {
    // Exponential: y = a * e^(b*x) or y = a * b^x
    // Transform to linear: ln(y) = ln(a) + b*x
    const logY = y.map(yi => yi > 0 ? Math.log(yi) : null).filter(ly => ly !== null);
    if (logY.length < 2) {
      return null;
    }
    
    const n = x.length;
    const sumX = x.reduce((a, b) => a + b, 0);
    const sumLogY = logY.reduce((a, b) => a + b, 0);
    const sumXLogY = x.slice(0, logY.length).reduce((sum, xi, i) => sum + xi * logY[i], 0);
    const sumX2 = x.slice(0, logY.length).reduce((sum, xi) => sum + xi * xi, 0);
    
    const det = logY.length * sumX2 - sumX * sumX;
    if (Math.abs(det) < 1e-10) {
      return null;
    }
    
    const lnA = (sumLogY * sumX2 - sumX * sumXLogY) / det;
    const b = (logY.length * sumXLogY - sumX * sumLogY) / det;
    const a = Math.exp(lnA);
    
    fitFunction = (t) => a * Math.exp(b * t);
  } else if (regressionType === 'logarithmic') {
    // Logarithmic: y = a + b*ln(x)
    // Filter out non-positive x values
    const validIndices = [];
    for (let i = 0; i < x.length; i++) {
      if (x[i] > 0) {
        validIndices.push(i);
      }
    }
    
    if (validIndices.length < 2) {
      return null;
    }
    
    const logX = validIndices.map(i => Math.log(x[i]));
    const validY = validIndices.map(i => y[i]);
    const n = validIndices.length;
    
    const sumLogX = logX.reduce((a, b) => a + b, 0);
    const sumY = validY.reduce((a, b) => a + b, 0);
    const sumLogXY = logX.reduce((sum, lxi, i) => sum + lxi * validY[i], 0);
    const sumLogX2 = logX.reduce((sum, lxi) => sum + lxi * lxi, 0);
    
    const det = n * sumLogX2 - sumLogX * sumLogX;
    if (Math.abs(det) < 1e-10) {
      return null;
    }
    
    const a = (sumY * sumLogX2 - sumLogX * sumLogXY) / det;
    const b = (n * sumLogXY - sumLogX * sumY) / det;
    
    fitFunction = (t) => t > 0 ? a + b * Math.log(t) : null;
  } else if (regressionType === 'mono-exponential') {
    // Mono-exponential rise to plateau: y(t) = 1 + A(1 - e^(-k(t - t₀)))
    // For normalized data, baseline ≈ 1, so we fix t₀ = 0 (or first timepoint)
    // y(t) = 1 + A(1 - e^(-k*t))
    
    // Find t₀ as the first timepoint (or 0 if data starts at 0)
    const t0 = Math.min(...x);
    const xShifted = x.map(xi => xi - t0);
    
    // Initial parameter estimates
    // A: amplitude = max(y) - 1 (since baseline is normalized to 1)
    const yMax = Math.max(...y);
    const yMin = Math.min(...y);
    let A = Math.max(0, yMax - 1); // Amplitude should be positive
    
    // k: rate constant - estimate from rise time
    // Find timepoint where signal reaches ~63%% of amplitude (1 - 1/e ≈ 0.632)
    const timeSpan = Math.max(...xShifted) - Math.min(...xShifted);
    const targetY = 1 + A * 0.632;
    let k = 0.01; // Default initial guess
    if (A > 0.01 && timeSpan > 0) {
      // Find closest point to target
      let closestIdx = 0;
      let minDiff = Math.abs(y[0] - targetY);
      for (let i = 1; i < y.length; i++) {
        const diff = Math.abs(y[i] - targetY);
        if (diff < minDiff) {
          minDiff = diff;
          closestIdx = i;
        }
      }
      // Estimate k from: 1 - e^(-k*t) = 0.632 => k ≈ 1/t
      if (xShifted[closestIdx] > 0) {
        k = 1 / xShifted[closestIdx];
        k = Math.max(1e-6, Math.min(10, k)); // Bound between 1e-6 and 10
      }
    }
    
    // Use simple iterative nonlinear least squares (Gauss-Newton)
    // Minimize sum of squared residuals: Σ(y - (1 + A(1 - e^(-k*t))))²
    const maxIterations = 100;
    const tolerance = 1e-6;
    
    for (let iter = 0; iter < maxIterations; iter++) {
      // Calculate residuals and Jacobian
      let sumResidual = 0;
      let sumResidualA = 0;
      let sumResidualK = 0;
      let sumJAA = 0;
      let sumJAK = 0;
      let sumJKK = 0;
      
      for (let i = 0; i < xShifted.length; i++) {
        const t = xShifted[i];
        const yObs = y[i];
        const yPred = 1 + A * (1 - Math.exp(-k * t));
        const residual = yObs - yPred;
        
        // Partial derivatives
        const dA = 1 - Math.exp(-k * t);
        const dK = A * t * Math.exp(-k * t);
        
        sumResidual += residual;
        sumResidualA += residual * dA;
        sumResidualK += residual * dK;
        sumJAA += dA * dA;
        sumJAK += dA * dK;
        sumJKK += dK * dK;
      }
      
      // Solve normal equations: J^T * J * delta = J^T * residual
      const det = sumJAA * sumJKK - sumJAK * sumJAK;
      if (Math.abs(det) < 1e-10) break;
      
      const deltaA = (sumResidualA * sumJKK - sumResidualK * sumJAK) / det;
      const deltaK = (sumResidualK * sumJAA - sumResidualA * sumJAK) / det;
      
      // Update parameters with damping
      const damping = 0.5;
      A = Math.max(0, A + damping * deltaA);
      k = Math.max(1e-6, Math.min(10, k + damping * deltaK));
      
      // If parameters didn't change much, we're done
      if (Math.abs(deltaA) < tolerance && Math.abs(deltaK) < tolerance) {
        break;
      }
    }
    
    // Store parameters for info panel
    const tau = 1 / k; // time constant
    const tHalf = Math.log(2) / k; // half-time
    const plateau = 1 + A; // plateau value
    
    // Create fit function
    fitFunction = (t) => {
      const tShifted = t - t0;
      if (tShifted < 0) return 1; // Before activation
      return 1 + A * (1 - Math.exp(-k * tShifted));
    };
    
    // Attach parameters to function for info panel
    fitFunction.params = {
      A: A,
      k: k,
      t0: t0,
      tau: tau,
      tHalf: tHalf,
      plateau: plateau
    };
  } else {
    return null;
  }
  
  return fitFunction;
}

function normalizeMeanAndError(timePoints, means, errors, baselineEndTime = 24, cutoffTime = 360) {
  // Normalize mean and error data: calculate baseline from 0-baselineEndTime, then plot only points after baselineEndTime
  // divided by baseline, up to cutoffTime. Baseline period is NOT included in the plot.
  // Returns { filteredTimePoints, normalizedMeans, normalizedErrors } - filtered arrays
  // First calculate baseline from means
  const baselineValues = [];
  for (let i = 0; i < timePoints.length; i++) {
    if (timePoints[i] !== null && means[i] !== null && timePoints[i] <= baselineEndTime) {
      baselineValues.push(means[i]);
    }
  }
  
  if (baselineValues.length === 0) {
    return { filteredTimePoints: timePoints, normalizedMeans: means, normalizedErrors: errors };
  }
  
  const baseline = baselineValues.reduce((a, b) => a + b, 0) / baselineValues.length;
  
  if (baseline === 0) {
    return { filteredTimePoints: timePoints, normalizedMeans: means, normalizedErrors: errors };
  }
  
  // Filter and normalize: only include points where baselineEndTime < time <= cutoffTime
  const filteredTimePoints = [];
  const normalizedMeans = [];
  const normalizedErrors = [];
  
  for (let i = 0; i < timePoints.length; i++) {
    const time = timePoints[i];
    const mean = means[i];
    const error = errors[i];
    
    if (time === null || mean === null) {
      continue;
    }
    
    // Only include points after baselineEndTime and up to cutoffTime (baseline NOT included)
    if (time > baselineEndTime && time <= cutoffTime) {
      // Subtract baselineEndTime to bring first point to x=0 (time after injection)
      filteredTimePoints.push(time - baselineEndTime);
      normalizedMeans.push(mean / baseline);
      // Errors scale proportionally
      normalizedErrors.push(error !== null ? error / baseline : null);
    }
  }
  
  return { filteredTimePoints, normalizedMeans, normalizedErrors };
}

function updateChart() {
  // Guard: ensure viewerData is loaded
  if (!viewerData || !allWellsData || Object.keys(allWellsData).length === 0) {
    console.warn("viewerData not loaded yet, skipping chart update");
    return;
  }
  
  const ctx = document.getElementById("well-chart").getContext("2d");

  // STEP 1: Check if normalization using fitted baseline is enabled (must be first)
  const normalizeBaselineToggle = document.getElementById("normalize-baseline-toggle");
  const normalizeBaseline = normalizeBaselineToggle ? normalizeBaselineToggle.checked : false;
  
  // STEP 2: Check if triplicate grouping is enabled (requires Step 1)
  const groupTriplicatesToggle = document.getElementById("group-triplicates-toggle");
  const groupTriplicates = groupTriplicatesToggle ? groupTriplicatesToggle.checked : false;
  
  // Enforce workflow: Step 2 (grouping) requires Step 1 (normalization)
  if (groupTriplicates && !normalizeBaseline) {
    // Disable grouping if normalization is not enabled
    if (groupTriplicatesToggle) {
      groupTriplicatesToggle.checked = false;
    }
    groupTriplicates = false;
  }
  
  // STEP 0: Check if individual wells should be shown
  // Step 0 only shows when NO other processing steps (1-2) are enabled
  const showIndividualWellsToggle = document.getElementById("show-individual-wells-toggle");
  const step0Requested = showIndividualWellsToggle ? showIndividualWellsToggle.checked : true;
  const anyProcessingStepEnabled = normalizeBaseline || groupTriplicates;
  // Show Step 0 if requested AND no other steps enabled
  // If no steps are enabled at all, default to showing Step 0 (raw data) as fallback
  const showIndividualWells = (step0Requested && !anyProcessingStepEnabled) || (!anyProcessingStepEnabled);
  
  // Show/hide Step 0 warning
  const step0Warning = document.getElementById("step0-warning");
  if (step0Warning) {
    step0Warning.style.display = (step0Requested && anyProcessingStepEnabled) ? "block" : "none";
  }
  
  // Disable Step 0 checkbox if other steps (1-2) are enabled (visual feedback)
  if (showIndividualWellsToggle) {
    showIndividualWellsToggle.disabled = anyProcessingStepEnabled;
  }
  
  // Get normalization baseline parameters
  const normalizeBaselineEndInput = document.getElementById("normalize-baseline-end-time");
  const normalizeBaselineEndTime = normalizeBaselineEndInput ? parseFloat(normalizeBaselineEndInput.value) || 24 : 24;
  
  // Get baseline fitting method and parameters
  const baselineMethodSelect = document.getElementById("baseline-method-select");
  const baselineMethod = baselineMethodSelect ? baselineMethodSelect.value : "lowess";
  const baselineFracInput = document.getElementById("baseline-frac-input");
  const baselineFrac = baselineFracInput ? parseFloat(baselineFracInput.value) || 0.5 : 0.5;
  const baselinePolyOrderInput = document.getElementById("baseline-poly-order-input");
  const baselinePolyOrder = baselinePolyOrderInput ? parseInt(baselinePolyOrderInput.value) || 1 : 1;
  
  // Get normalization mode
  const normalizationModeSelect = document.getElementById("normalization-mode-select");
  const normalizationMode = normalizationModeSelect ? normalizationModeSelect.value : "delta_f_over_f";
  
  // Show/hide normalization baseline parameters
  const normalizeBaselineParametersDiv = document.getElementById("normalize-baseline-parameters");
  if (normalizeBaselineParametersDiv && normalizeBaselineToggle) {
    normalizeBaselineParametersDiv.style.display = normalizeBaseline ? "block" : "none";
  }
  
  // STEP 3: Check if cutoff is enabled
  const cutoffBaselineToggle = document.getElementById("cutoff-baseline-toggle");
  const cutoffBaseline = cutoffBaselineToggle ? cutoffBaselineToggle.checked : false;
  const cutoffBaselineTimeInput = document.getElementById("cutoff-baseline-time");
  const cutoffBaselineTime = cutoffBaselineTimeInput ? parseFloat(cutoffBaselineTimeInput.value) || 360 : 360;
  
  // Show/hide cutoff parameters
  const cutoffBaselineParametersDiv = document.getElementById("cutoff-baseline-parameters");
  if (cutoffBaselineParametersDiv && cutoffBaselineToggle) {
    cutoffBaselineParametersDiv.style.display = cutoffBaseline ? "block" : "none";
  }
  
  // STEP 4: Check if regression is enabled
  const fitRegressionToggle = document.getElementById("fit-regression-toggle");
  const fitRegression = fitRegressionToggle ? fitRegressionToggle.checked : false;
  const regressionTypeSelect = document.getElementById("regression-type-select");
  const regressionType = regressionTypeSelect ? regressionTypeSelect.value : "mono-exponential";
  const polynomialOrderInput = document.getElementById("polynomial-order-input");
  const polynomialOrder = polynomialOrderInput ? parseInt(polynomialOrderInput.value) || 2 : 2;
  
  // Show/hide regression parameters
  const fitRegressionParametersDiv = document.getElementById("fit-regression-parameters");
  if (fitRegressionParametersDiv && fitRegressionToggle) {
    fitRegressionParametersDiv.style.display = fitRegression ? "block" : "none";
  }
  
  // Update regression type params to show/hide info panel when regression is toggled
  if (fitRegression) {
    updateRegressionTypeParams(true); // Skip chart update since we're already updating
  } else {
    // Hide info panel when regression is disabled
    const infoPanel = document.getElementById("regression-info-panel");
    if (infoPanel) {
      infoPanel.style.display = "none";
    }
  }

  if (selectedWells.size === 0) {
    // Clear the chart if no wells selected
    // Determine y-axis label based on processing steps
    let yAxisLabel = "Fluorescence (FI)";
    if (normalizeBaseline) {
      if (normalizationMode === 'delta_f_over_f') {
        yAxisLabel = "ΔF/F (Normalized)";
      } else {
        yAxisLabel = "F/F₀ (Normalized)";
      }
    }
    if (chart) {
      // Transform existing chart
      chart.data.datasets = [];
      chart.options.scales.x.title.text = "Time (s)";
      chart.options.scales.y.title.text = yAxisLabel;
      chart.options.scales.y.min = undefined;
      chart.options.plugins.title.text = "Select wells to display data";
      chart.update('none'); // Update without animation
    } else {
      // Create new chart only if it doesn't exist
      chart = new Chart(ctx, {
        type: "line",
        data: { datasets: [] },
        options: {
          responsive: true,
          maintainAspectRatio: false,
          scales: {
            x: { 
              type: "linear", 
              title: { display: true, text: "Time (s)" } 
            },
            y: { 
              title: { display: true, text: yAxisLabel },
              min: undefined
            }
          },
          plugins: {
            legend: { display: false },
            title: { display: true, text: "Select wells to display data" }
          }
        }
      });
    }
    return;
  }

  const datasets = [];
  const colors = [
    "rgba(54, 162, 235, 1.0)",
    "rgba(255, 99, 132, 1.0)",
    "rgba(255, 206, 86, 1.0)",
    "rgba(75, 192, 192, 1.0)",
    "rgba(153, 102, 255, 1.0)",
    "rgba(255, 159, 64, 1.0)",
    "rgba(199, 199, 199, 1.0)",
    "rgba(83, 102, 255, 1.0)"
  ];
  
  let colorIdx = 0;
  
  // STEP 0: Add individual wells (raw data, before any processing) if enabled
  if (showIndividualWells) {
    Array.from(selectedWells).forEach(key => {
      const [plateId, wellId] = key.split(":");
      const column = getColumnFromWellId(wellId);
      
      // Filter by enabledColumns if any are set
      if (enabledColumns.size > 0 && !enabledColumns.has(column)) {
        return;
      }
      
      // Filter by enabledWells if any are set
      if (enabledWells.size > 0 && !enabledWells.has(wellId)) {
        return;
      }
      
      if (allWellsData[plateId] && allWellsData[plateId][wellId]) {
        const wellData = allWellsData[plateId][wellId];
        const plate = viewerData.plates.find(p => p.id === plateId);
        const filename = plate ? (plate.file || plate.id) : plateId;
        const config = viewerData.config || {};
        const wellLabels = config.well_labels || {};
        const rowLetter = getRowLetter(wellId);
        const isControl = isControlWell(wellId);
        const label = wellLabels[rowLetter] || wellData.label || wellData.content || `Well ${wellId}`;
        const formattedLabel = formatLabel(label);
        
        // Use raw data (not processed) for Step 0
        const color = colors[colorIdx %% colors.length];
        colorIdx++;
        
        // Add "Control" suffix for control wells
        const displayLabel = isControl 
          ? `${wellId} - ${formattedLabel} Control Col ${column}`
          : `${wellId} - ${formattedLabel} Col ${column}`;
        
        datasets.push({
          label: displayLabel,
          data: wellData.time_s.map((t, i) => ({ x: t, y: wellData.values[i] })).filter(d => d.y !== null && d.x !== null),
          borderColor: color.replace("1.0", "0.5"), // Lighter color for raw data
          backgroundColor: "transparent",
          tension: 0.15,
          pointRadius: 1.5,
          borderWidth: 1,
          borderDash: [3, 3], // Dashed line for raw data
          fill: false
        });
      }
    });
  }
  

  // Group selected wells by triplicate group or individual well
  const groupedSelections = {};
  const selectedDatasets = getSelectedDatasets();
  
  // groupTriplicates already defined above
  
  // First pass: identify all triplicate groups from selected wells
  // Only group into triplicates if Step 2 is enabled
  const triplicateGroupsFound = new Map(); // groupKey -> { name, column, wells: Set of wellIds }
  
  Array.from(selectedWells).forEach(key => {
    const [plateId, wellId] = key.split(":");
    const column = getColumnFromWellId(wellId);
    
    // Filter by enabledColumns if any are set
    if (enabledColumns.size > 0 && !enabledColumns.has(column)) {
      return; // Skip wells from disabled columns
    }
    
    // Filter by enabledWells if any are set (exclude wells)
    if (enabledWells.size > 0 && !enabledWells.has(wellId)) {
      return; // Skip excluded wells
    }
    
    const isControl = isControlWell(wellId);
    
    // Controls are always handled as individual wells (no grouping)
    if (isControl) {
      return; // Skip controls - they'll be handled as individual wells in the third pass
    }
    
    // STEP 2: Only group into triplicates if grouping is enabled
    if (groupTriplicates) {
      const triplicateGroup = getTriplicateGroup(wellId, column);
      if (triplicateGroup) {
        // Create group key that includes column to separate same pattern across columns
        const groupKey = `${triplicateGroup.name}_col${column}`;
        if (!triplicateGroupsFound.has(groupKey)) {
          // Get row letters from the triplicate group pattern
          const groupRowLetters = getRowLettersFromWells(triplicateGroup.wells);
          // Create well IDs for this column using the row pattern
          const wellIdsForColumn = groupRowLetters.map(row => row + column);
          
          triplicateGroupsFound.set(groupKey, {
            name: triplicateGroup.name,
            column: column,
            wellIds: new Set(wellIdsForColumn)
          });
        }
      }
    }
    // If grouping is disabled, wells will be handled as individual wells in the third pass
  });
  
  // Second pass: build grouped selections for triplicate groups
  triplicateGroupsFound.forEach((groupInfo, groupKey) => {
    groupedSelections[groupKey] = {
      type: "triplicate",
      wells: [],
      name: groupInfo.name
    };
    
    // Collect wells from datasets that have selected wells in this triplicate group
    Array.from(selectedWells).forEach(key => {
      const [plateId, wellId] = key.split(":");
      const column = getColumnFromWellId(wellId);
      if (groupInfo.wellIds.has(wellId) && column === groupInfo.column) {
        const wells = allWellsData[plateId] || {};
        if (wells[wellId]) {
          // Check if already added
          const exists = groupedSelections[groupKey].wells.some(
            w => w.wellId === wellId && w.plateId === plateId
          );
          if (!exists) {
            // STEP 1: Apply normalization using fitted baseline if enabled (must happen before grouping)
            let processedWellData = wells[wellId];
            if (normalizeBaseline) {
              processedWellData = normalizeWithFittedBaseline(
                processedWellData,
                normalizeBaselineEndTime,
                baselineMethod,
                baselineFrac,
                baselinePolyOrder,
                normalizationMode
              );
              if (!processedWellData) return; // Skip if normalization failed
            }
            // Add to group (Step 2 happens here - grouping)
            // Step 3 (averaging) will be applied later
            groupedSelections[groupKey].wells.push({
              wellId: wellId,
              plateId: plateId,
              data: processedWellData
            });
          }
        }
      }
    });
  });
  
  // Third pass: handle individual wells (skip if already in triplicate group)
  Array.from(selectedWells).forEach(key => {
    // Verify this well is still selected (may have been deselected)
    if (!selectedWells.has(key)) {
      return; // Skip wells that are no longer selected
    }
    
    const [plateId, wellId] = key.split(":");
    const column = getColumnFromWellId(wellId);
    const isControl = isControlWell(wellId);
    
    // Controls are always handled as individual wells (no grouping)
    if (isControl) {
      // Continue to process as individual well
    } else {
      // STEP 2: Check if this well is part of a triplicate group we've already processed
      // Only skip if grouping is enabled and this well was grouped
      if (groupTriplicates) {
        const triplicateGroup = getTriplicateGroup(wellId, column);
        if (triplicateGroup) {
          // Already handled in triplicateGroupsFound - skip
          return;
        }
      }
    }
    
    // Filter by enabledWells if any are set
    if (enabledWells.size > 0 && !enabledWells.has(wellId)) {
      return; // Skip excluded wells
    }
    
    // Individual well (not part of a triplicate group)
    const groupKey = `well_${wellId}`;
    if (!groupedSelections[groupKey]) {
      const config = viewerData.config || {};
      const wellLabels = config.well_labels || {};
      // Get label from row-based config or well data
      const rowLetter = getRowLetter(wellId);
      let label = wellLabels[rowLetter];
      if (!label) {
        for (const pId of selectedDatasets) {
          if (allWellsData[pId] && allWellsData[pId][wellId]) {
            label = allWellsData[pId][wellId].label || allWellsData[pId][wellId].content || `Well ${wellId}`;
            break;
          }
        }
      }
      groupedSelections[groupKey] = {
        type: "single",
        wells: [],
        label: formatLabel(label) || `Well ${wellId}`
      };
    }
    
    // Add this well's data if it exists
    if (allWellsData[plateId] && allWellsData[plateId][wellId]) {
      const existing = groupedSelections[groupKey].wells.find(w => w.wellId === wellId && w.plateId === plateId);
      if (!existing) {
          // STEP 1: Apply normalization using fitted baseline if enabled (must happen before grouping)
          let processedWellData = allWellsData[plateId][wellId];
          if (normalizeBaseline) {
            processedWellData = normalizeWithFittedBaseline(
              processedWellData,
              normalizeBaselineEndTime,
              baselineMethod,
              baselineFrac,
              baselinePolyOrder,
              normalizationMode
            );
            if (!processedWellData) return; // Skip if normalization failed
          }
          // STEP 4: Apply normalization using fitted baseline if enabled
          if (normalizeBaseline) {
            processedWellData = normalizeWithFittedBaseline(
              processedWellData,
              normalizeBaselineEndTime,
              baselineMethod,
              baselineFrac,
              baselinePolyOrder,
              normalizationMode
            );
            if (!processedWellData) return; // Skip if normalization failed
          }
          // Add to group (Step 2 happens here - grouping)
          // Step 3 (averaging) will be applied later
          groupedSelections[groupKey].wells.push({
          wellId, 
          plateId, 
          data: processedWellData
        });
      }
    }
  });
  
  // Continue with grouped selections (colorIdx already initialized for Step 0)
  // Process grouped selections when:
  // 1. Step 0 is not showing (because other steps are enabled), OR
  // 2. Processing steps are enabled (which automatically disables Step 0)
  // This ensures data always shows - either Step 0 individual wells OR grouped selections
  if (!showIndividualWells) {
    Object.keys(groupedSelections).forEach(groupKey => {
      const group = groupedSelections[groupKey];
      
      // Filter out any wells that are no longer selected or are excluded
      group.wells = group.wells.filter(w => {
        const key = `${w.plateId}:${w.wellId}`;
        // Check if still selected
        if (!selectedWells.has(key)) {
          return false;
        }
        // Check if excluded via enabledWells
        if (enabledWells.size > 0 && !enabledWells.has(w.wellId)) {
          return false;
        }
        return true;
      });
      
      // Skip groups with no wells after filtering
      if (group.wells.length === 0) {
        return; // Skip to next group
      }
      
      const color = colors[colorIdx %% colors.length];
      colorIdx++;

    if (group.type === "triplicate" && group.wells.length > 0) {
      // STEP 3: Average within triplicate group
      const result = processWellsAndCalculateMean(
        group.wells.map(w => w.data),
        false, // No normalization
        0, // baselineEndTime not used
        999999 // cutoffTime not used
      );
      const timePoints = result.filteredTimePoints;
      const means = result.normalizedMeans;
      const errors = result.normalizedErrors;
      
      // Create dataset with mean and error information
      const lineData = timePoints.map((t, i) => {
        const point = {
          x: t,
          y: means[i],
          error: errors[i]
        };
        return point;
      }).filter(d => d.y !== null);
      
      // Get group name from first well's plate
      const firstPlate = viewerData.plates.find(p => p.id === group.wells[0].plateId);
      const columnGroupName = firstPlate ? firstPlate.column_group : "";
      const datasetCount = new Set(group.wells.map(w => w.plateId)).size;
      const wellCount = group.wells.length;
      
      // Create appropriate label based on number of replicates
      let replicateLabel = "";
      if (wellCount === 3) {
        replicateLabel = "triplicate";
      } else if (wellCount === 6) {
        replicateLabel = "hexplicate";
      } else if (wellCount === 9) {
        replicateLabel = "nonuplicate";
      } else {
        replicateLabel = `${wellCount}-plicate`;
      }
      
      // Get column number from first well
      const firstWellId = group.wells[0].wellId;
      const column = getColumnFromWellId(firstWellId);
      
      // Add well IDs when Step 1 (normalization) is enabled
      let wellIdsPrefix = "";
      if (normalizeBaseline && group.wells.length > 0) {
        const wellIds = group.wells.map(w => w.wellId).join(", ");
        wellIdsPrefix = `${wellIds} - `;
      }
      
      const labelSuffix = datasetCount > 1 
        ? ` (${replicateLabel}, ${datasetCount} datasets, mean ± SE)` 
        : ` (${replicateLabel}, mean ± SE)`;
      const formattedName = formatLabel(group.name);
      const label = columnGroupName 
        ? `${wellIdsPrefix}${formattedName} (${columnGroupName}) Col ${column}${labelSuffix}`
        : `${wellIdsPrefix}${formattedName} Col ${column}${labelSuffix}`;
      
      datasets.push({
        label: label,
        data: lineData,
        borderColor: color,
        backgroundColor: color.replace("1.0", "0.2"),
        tension: 0.15,
        pointRadius: 3,
        borderWidth: 2,
        fill: false
      });
    } else if (group.type === "single") {
      // Single well - check if it's a duplicate
      const wellId = group.wells[0].wellId;
      const isDuplicate = isDuplicateWell(wellId);
      
      if (group.wells.length > 1 && isDuplicate) {
        // Duplicate well - show each dataset separately for comparison
        const lineStyles = ["solid", "dashed", "dotted"];
        const formattedLabel = formatLabel(group.label);
        
        group.wells.forEach((well, idx) => {
          const plate = viewerData.plates.find(p => p.id === well.plateId);
          const filename = plate ? (plate.file || plate.id) : well.plateId;
          const lineStyle = lineStyles[idx %% lineStyles.length];
          
          // Use slightly different shades of the same color for duplicates
          const baseColor = color;
          const alpha = 1.0 - (idx * 0.15); // Slightly fade each subsequent line
          const adjustedColor = baseColor.replace("1.0", alpha.toFixed(1));
          
          let timePoints = well.data.time_s;
          let values = well.data.values;
          
          const column = getColumnFromWellId(well.wellId);
          // Add well ID when Step 1 (normalization) is enabled
          const wellIdPrefix = normalizeBaseline ? `${well.wellId} - ` : "";
          datasets.push({
            label: `${wellIdPrefix}${formattedLabel} Col ${column} (${filename})`,
            data: timePoints.map((t, i) => {
              const point = { x: t, y: values[i] };
              return point;
            }).filter(d => d.y !== null),
            borderColor: adjustedColor,
            backgroundColor: "transparent",
            borderDash: lineStyle === "dashed" ? [5, 5] : lineStyle === "dotted" ? [2, 2] : [],
            tension: 0.15,
            pointRadius: 2,
            borderWidth: 1.5,
            fill: false,
          });
        });
      } else if (group.wells.length > 1) {
        // STEP 3: Average multiple datasets for same well (not duplicate) - calculate mean and error
        const result = processWellsAndCalculateMean(
          group.wells.map(w => w.data),
          false, // No normalization
          0, // baselineEndTime not used
          999999 // cutoffTime not used
        );
        const timePoints = result.filteredTimePoints;
        const means = result.normalizedMeans;
        const errors = result.normalizedErrors;
        
        const lineData = timePoints.map((t, i) => {
          const point = {
            x: t,
            y: means[i],
            error: errors[i]
          };
          return point;
        }).filter(d => d.y !== null);
        
        const datasetCount = group.wells.length;
        const column = getColumnFromWellId(wellId);
        const formattedLabel = formatLabel(group.label);
        // Add well ID when Step 1 (normalization) is enabled
        const wellIdPrefix = normalizeBaseline ? `${wellId} - ` : "";
        const label = `${wellIdPrefix}${formattedLabel} Col ${column} (${datasetCount} datasets, mean ± SE)`;
        
        datasets.push({
          label: label,
          data: lineData,
          borderColor: color,
          backgroundColor: color.replace("1.0", "0.2"),
          tension: 0.15,
          pointRadius: 3,
          borderWidth: 2,
          fill: false
        });
      } else {
        // Single dataset for this well
        // Skip showing processed individual wells if Step 0 is already showing the raw data
        // But only skip if Step 0 is actually enabled (not just requested)
        if (showIndividualWells && !anyProcessingStepEnabled) {
          // Step 0 already shows this well, skip the processed version to avoid duplicates
          return; // Skip to next group in forEach
        }
        
        const well = group.wells[0];
        const column = getColumnFromWellId(well.wellId);
        const formattedLabel = formatLabel(group.label);
        
        let timePoints = well.data.time_s;
        let values = well.data.values;
        
        // Add well ID when Step 1 (normalization) is enabled
        const wellIdPrefix = normalizeBaseline ? `${well.wellId} - ` : "";
        const label = formattedLabel 
          ? `${wellIdPrefix}${formattedLabel} Col ${column}` 
          : `${wellIdPrefix}Well ${well.wellId} Col ${column}`;
        
        datasets.push({
          label: label,
          data: timePoints.map((t, i) => {
            const point = { x: t, y: values[i] };
            return point;
          }).filter(d => d.y !== null),
          borderColor: color,
          backgroundColor: "transparent",
          tension: 0.15,
          pointRadius: 2,
          borderWidth: 1.5,
        });
      }
    }
    });
  }

  // STEP 3: Apply cutoff filtering to all datasets
  let xAxisMin = undefined;
  let xAxisMax = undefined;
  
  if (cutoffBaseline) {
    datasets.forEach(dataset => {
      if (dataset.data && dataset.data.length > 0) {
        // Find first real timepoint (first non-null time)
        let firstTimepoint = null;
        for (let i = 0; i < dataset.data.length; i++) {
          if (dataset.data[i].x !== null && isFinite(dataset.data[i].x)) {
            firstTimepoint = dataset.data[i].x;
            break;
          }
        }
        
        // Determine minimum time to include
        // If normalization is enabled, exclude baseline period (points before baselineEndTime)
        // Otherwise, exclude points before first real timepoint
        const minTime = normalizeBaseline ? normalizeBaselineEndTime : (firstTimepoint !== null ? firstTimepoint : -Infinity);
        
        // Filter data points: exclude before minTime and after cutoff time
        dataset.data = dataset.data.filter(point => {
          if (point.x === null || !isFinite(point.x)) return false;
          if (point.x < minTime) return false;
          if (point.x > cutoffBaselineTime) return false;
          return true;
        });
        
        // Track min/max x values for axis adjustment
        if (dataset.data.length > 0) {
          const times = dataset.data.map(p => p.x).filter(t => t !== null && isFinite(t));
          if (times.length > 0) {
            const datasetMin = Math.min(...times);
            const datasetMax = Math.max(...times);
            if (xAxisMin === undefined || datasetMin < xAxisMin) {
              xAxisMin = datasetMin;
            }
            if (xAxisMax === undefined || datasetMax > xAxisMax) {
              xAxisMax = datasetMax;
            }
          }
        }
      }
    });
  } else {
    // Even if cutoff is disabled, calculate x-axis range from all datasets for proper display
    datasets.forEach(dataset => {
      if (dataset.data && dataset.data.length > 0) {
        const times = dataset.data.map(p => p.x).filter(t => t !== null && isFinite(t));
        if (times.length > 0) {
          const datasetMin = Math.min(...times);
          const datasetMax = Math.max(...times);
          if (xAxisMin === undefined || datasetMin < xAxisMin) {
            xAxisMin = datasetMin;
          }
          if (xAxisMax === undefined || datasetMax > xAxisMax) {
            xAxisMax = datasetMax;
          }
        }
      }
    });
  }
  
  // STEP 4: Add regression curves to datasets
  if (fitRegression) {
    const regressionDatasets = [];
    const monoExpParams = []; // Store parameters for mono-exponential fits
    
    datasets.forEach(dataset => {
      if (dataset.data && dataset.data.length >= 2) {
        // Check if this dataset represents a control well
        // Extract wellId from label (format: "N1 - Label" or "N1 - Label Control Col 1")
        let isControl = false;
        const label = dataset.label || "";
        
        // Method 1: Check if label contains "Control" (case-insensitive)
        if (label.toLowerCase().includes("control")) {
          isControl = true;
        } else {
          // Method 2: Try to extract wellId from label and check
          // Label format is typically: "N1 - Label Col 1" or "N1 - Label"
          const wellIdMatch = label.match(/^([A-Z]\d+)\s*-/);
          if (wellIdMatch) {
            const wellId = wellIdMatch[1];
            isControl = isControlWell(wellId);
          }
        }
        
        // Use linear regression for controls, otherwise use selected regression type
        const regressionTypeToUse = isControl ? 'linear' : regressionType;
        
        // Fit regression curve
        const fitFunction = fitRegressionCurve(dataset.data, regressionTypeToUse, polynomialOrder);
        if (fitFunction) {
          // Get time range from data
          const times = dataset.data.map(p => p.x).filter(t => t !== null && isFinite(t));
          if (times.length === 0) return;
          
          const minTime = Math.min(...times);
          const maxTime = Math.max(...times);
          
          // Generate regression curve points
          const numPoints = 100;
          const regressionData = [];
          for (let i = 0; i <= numPoints; i++) {
            const t = minTime + (maxTime - minTime) * (i / numPoints);
            const y = fitFunction(t);
            if (y !== null && isFinite(y)) {
              regressionData.push({ x: t, y: y });
            }
          }
          
          if (regressionData.length > 0) {
            // Create regression dataset with same color but different style
            const regressionColor = dataset.borderColor || dataset.backgroundColor || "rgba(128, 128, 128, 1.0)";
            regressionDatasets.push({
              label: `${dataset.label} (${regressionTypeToUse} fit)`,
              data: regressionData,
              borderColor: regressionColor,
              backgroundColor: "transparent",
              borderWidth: 2,
              borderDash: [5, 5], // Dashed line for regression
              pointRadius: 0, // No points on regression line
              tension: 0,
              fill: false
            });
            
            // Store mono-exponential parameters for info panel
            // Only store if using mono-exponential (not for controls which use linear)
            if (regressionTypeToUse === 'mono-exponential' && fitFunction.params) {
              // Remove replicate and dataset info from label for regression lines
              // Pattern: remove "(X-plicate, Y datasets, mean ± SE)" or similar patterns
              let cleanLabel = dataset.label;
              // Remove patterns containing plicate, datasets, or mean ± SE in parentheses
              cleanLabel = cleanLabel.replace(/\s*\([^)]*(?:plicate|datasets|mean\s*±\s*SE)[^)]*\)/gi, '');
              // Clean up any double spaces or trailing/leading spaces
              cleanLabel = cleanLabel.replace(/\s+/g, ' ').trim();
              
              monoExpParams.push({
                label: cleanLabel,
                params: fitFunction.params,
                color: regressionColor // Use the same color as the regression curve
              });
            }
          }
        }
      }
    });
    
    // Add regression datasets to the main datasets array
    datasets.push(...regressionDatasets);
    
    // Update info panel for mono-exponential
    if (regressionType === 'mono-exponential') {
      updateMonoExponentialInfoPanel(monoExpParams);
    } else {
      const infoPanel = document.getElementById("regression-info-panel");
      if (infoPanel) {
        infoPanel.style.display = "none";
      }
    }
    
    // Recalculate x-axis range to include regression curves
    if (cutoffBaseline) {
      datasets.forEach(dataset => {
        if (dataset.data && dataset.data.length > 0) {
          const times = dataset.data.map(p => p.x).filter(t => t !== null && isFinite(t));
          if (times.length > 0) {
            const datasetMin = Math.min(...times);
            const datasetMax = Math.max(...times);
            if (xAxisMin === undefined || datasetMin < xAxisMin) {
              xAxisMin = datasetMin;
            }
            if (xAxisMax === undefined || datasetMax > xAxisMax) {
              xAxisMax = datasetMax;
            }
          }
        }
      });
    }
  }

  // Determine y-axis label based on processing steps
  let yAxisLabel = "Fluorescence (FI)";
  if (normalizeBaseline) {
    if (normalizationMode === 'delta_f_over_f') {
      yAxisLabel = "ΔF/F (Normalized)";
    } else {
      yAxisLabel = "F/F₀ (Normalized)";
    }
  }

  // Transform existing chart instead of recreating
  if (chart) {
    // Clear existing datasets
    while (chart.data.datasets.length > 0) {
      chart.data.datasets.pop();
    }
    
    // Add new datasets
    datasets.forEach(dataset => {
      chart.data.datasets.push(dataset);
    });
    
    // Update chart options
    chart.options.scales.x.title.text = "Time (s)";
    chart.options.scales.y.title.text = yAxisLabel;
    chart.options.scales.y.min = undefined;
    
    // Adjust x-axis range based on filtered data (Step 3)
    if (cutoffBaseline && xAxisMin !== undefined && xAxisMax !== undefined) {
      // Add small padding (2%% of range) for better visualization
      const padding = (xAxisMax - xAxisMin) * 0.02;
      chart.options.scales.x.min = Math.max(0, xAxisMin - padding);
      chart.options.scales.x.max = xAxisMax + padding;
    } else {
      // Reset to auto if cutoff is disabled
      chart.options.scales.x.min = undefined;
      chart.options.scales.x.max = undefined;
    }
    
    // Count only wells that are actually displayed (accounting for exclusions)
    // Count unique well groups/datasets that are actually displayed
    const displayedWellCount = datasets.length;
    chart.options.plugins.title.text = `Time Course - ${displayedWellCount} well(s) selected`;
    chart.options.plugins.legend.display = datasets.length > 0;
    
    // Force Chart.js to update - use update() without mode to ensure proper redraw
    chart.update();
  } else {
    // Create new chart only if it doesn't exist
    chart = new Chart(ctx, {
      type: "line",
      data: {
        datasets: datasets
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        parsing: false,
        scales: {
            x: { 
              type: "linear", 
              title: { 
                display: true,
                text: "Time (s)"
              },
              min: cutoffBaseline && xAxisMin !== undefined ? Math.max(0, xAxisMin - (xAxisMax - xAxisMin) * 0.02) : undefined,
              max: cutoffBaseline && xAxisMax !== undefined ? xAxisMax + (xAxisMax - xAxisMin) * 0.02 : undefined
            },
            y: { 
              title: { 
                display: true,
                text: yAxisLabel
            },
            min: undefined
          }
        },
        plugins: {
          legend: {
            display: datasets.length > 0, // Show legend when there are datasets
            position: "top",
            onClick: null,
            labels: {
              usePointStyle: true,
              padding: 15,
              font: {
                size: 11
              },
              generateLabels: function(chart) {
                // Custom label generation to show color and well name clearly
                const original = Chart.defaults.plugins.legend.labels.generateLabels;
                const labels = original.call(this, chart);
                // Labels already have the correct text from dataset.label
                return labels;
              }
            }
          },
          title: {
            display: true,
            text: `Time Course - ${datasets.length} well(s) selected`
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                let label = context.dataset.label || "";
                if (label) {
                  label += ": ";
                }
                if (context.parsed.y !== null) {
                  label += context.parsed.y.toFixed(2);
                  // Show error if available
                  if (context.raw && context.raw.error !== undefined && context.raw.error !== null) {
                    label += " ± " + context.raw.error.toFixed(2);
                  }
                }
                return label;
              }
            }
          }
        }
      }
    });
  }
}

// Initialize settings display on page load
updateCurrentSettingsDisplay();
// Initialize baseline method params to show/hide LOWESS/polynomial parameters based on default selection
updateBaselineMethodParams();
// Initialize regression type params to show/hide info panel based on default selection
// Skip chart update during initialization
updateRegressionTypeParams(true);

loadData().catch(err => {
  console.error("Failed to load viewer_data.json:", err);
  alert("Could not load viewer_data.json. Did you run the Python script in the same folder?");
});
</script>
</body>
</html>
"""


def write_html(html_path: str):
    html = HTML_TEMPLATE % {
        "rows": PLATE_ROWS,
        "cols": PLATE_COLS,
    }
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)


def write_web_command(script_dir: str, web_dir: str):
    """
    Create a double-clickable 'web.command' file in the main folder.
    On macOS, .command files can be double-clicked to run in Terminal.
    """
    web_command_path = os.path.join(script_dir, "web.command")
    
    web_command_content = f"""#!/bin/bash
# Double-clickable web server launcher for Plate Viewer
# This file is auto-generated by plate_viewer.py

PORT=8000
WEB_DIR="{web_dir}"

# Change to web directory
cd "$WEB_DIR" || {{
    echo "Error: Could not find web directory: $WEB_DIR"
    echo "Please run plate_viewer.py first to generate the web files."
    read -p "Press Enter to exit..."
    exit 1
}}

# Check if viewer_data.json exists
if [ ! -f "viewer_data.json" ]; then
    echo "Warning: viewer_data.json not found!"
    echo "Please run plate_viewer.py first to generate the data files."
    echo ""
fi

URL="http://localhost:${{PORT}}/index.html"

echo "============================================================"
echo "Plate Viewer Web Server"
echo "============================================================"
echo "Server running at: $URL"
echo "Serving directory: $WEB_DIR"
echo ""
echo "Press Ctrl+C to stop the server"
echo "============================================================"
echo ""

# Try to open browser automatically after a short delay
(sleep 2 && if command -v open &> /dev/null; then
    open "$URL" 2>/dev/null && echo "Opened $URL in your default browser" || echo "Please open $URL manually in your browser"
elif command -v xdg-open &> /dev/null; then
    xdg-open "$URL" 2>/dev/null && echo "Opened $URL in your default browser" || echo "Please open $URL manually in your browser"
else
    echo "Please open $URL manually in your browser"
fi) &

# Try Python 3 first, then Python 2, then exit with error
if command -v python3 &> /dev/null; then
    python3 -m http.server "$PORT"
elif command -v python &> /dev/null; then
    python -m SimpleHTTPServer "$PORT"
else
    echo "Error: Python not found. Please install Python to run a local server."
    read -p "Press Enter to exit..."
    exit 1
fi
"""
    
    with open(web_command_path, "w", encoding="utf-8") as f:
        f.write(web_command_content)
    
    # Make it executable
    os.chmod(web_command_path, 0o755)
    
    return web_command_path


# ============================================================================
# Baseline Identification and Normalization Functions
# ============================================================================

def identify_baseline_window(df: pd.DataFrame, baseline_end_time: float = 24.0) -> Dict[str, pd.DataFrame]:
    """
    For each well column, define the baseline window as all rows with Time_s <= baseline_end_time.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Long-format DataFrame with columns ['plate_id', 'well', 'content', 'time_s', 'value']
    baseline_end_time : float
        Maximum time (in seconds) for baseline window (default: 24.0)
    
    Returns:
    --------
    Dict[str, pd.DataFrame]
        Dictionary mapping well_id to DataFrame containing only baseline points for that well
    """
    baseline_data = {}
    
    for well_id, well_df in df.groupby('well'):
        baseline_mask = well_df['time_s'] <= baseline_end_time
        baseline_data[well_id] = well_df[baseline_mask].copy()
    
    return baseline_data


def fit_baseline(
    time: pd.Series,
    fluorescence: pd.Series,
    method: str = 'lowess',
    frac: float = 0.5,
    poly_order: int = 1,
    well_id: str = None
) -> pd.Series:
    """
    Fit a smooth baseline function on the baseline window and extrapolate to full timecourse.
    
    Parameters:
    -----------
    time : pd.Series
        Time values for baseline window
    fluorescence : pd.Series
        Fluorescence values for baseline window
    method : str
        Baseline fitting method: 'lowess', 'constant', or 'polynomial' (default: 'constant')
        - 'lowess': Locally Weighted Scatterplot Smoothing. Use when baseline has gradual drift.
          Best for non-linear baselines with smooth trends. Formula: g(t) = LOWESS(F_baseline, t_baseline, frac)
        - 'constant': Constant mean baseline. Use when baseline is stable/flat.
          Best for stable baselines with minimal drift. Formula: g(t) = mean(F_baseline)
        - 'polynomial': Polynomial fit. Use when baseline has linear or polynomial drift.
          Best for baselines with known polynomial trends. Formula: g(t) = c₀ + c₁t + c₂t² + ...
    frac : float
        Fraction of data used for LOWESS smoothing (default: 0.5). Higher = smoother, lower = more local.
    poly_order : int
        Polynomial order for polynomial fit (default: 1). 1 = linear, 2 = quadratic, etc.
    well_id : str, optional
        Well identifier for output
    
    Returns:
    --------
    pd.Series
        Fitted baseline values for all timepoints (same index as input time)
    """
    import numpy as np
    
    # Remove NaN values
    valid_mask = ~(pd.isna(time) | pd.isna(fluorescence))
    time_clean = time[valid_mask].values
    fluo_clean = fluorescence[valid_mask].values
    
    if len(time_clean) == 0:
        # No valid data, return NaN series
        return pd.Series(index=time.index, dtype=float)
    
    if method == 'constant':
        # 0th-order (constant mean baseline)
        # Formula: g(t) = mean(F_baseline) for all t
        baseline_value = np.mean(fluo_clean)
        baseline_fit = pd.Series(baseline_value, index=time.index)
    
    elif method == 'polynomial':
        # Polynomial fit using np.polyfit
        # Formula: g(t) = c₀ + c₁t + c₂t² + ... (depending on order)
        coeffs = np.polyfit(time_clean, fluo_clean, poly_order)
        poly_func = np.poly1d(coeffs)
        baseline_fit = pd.Series(poly_func(time.values), index=time.index)
    
    elif method == 'lowess':
        # LOWESS smoothing
        # Formula: g(t) = LOWESS(F_baseline, t_baseline, frac) interpolated/extrapolated
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            # LOWESS returns array of [x, y] pairs
            smoothed = lowess(fluo_clean, time_clean, frac=frac, return_sorted=False)
            # Interpolate/extrapolate to full timecourse
            from scipy.interpolate import interp1d
            interp_func = interp1d(
                time_clean, smoothed, 
                kind='linear', 
                bounds_error=False, 
                fill_value='extrapolate'
            )
            baseline_fit = pd.Series(interp_func(time.values), index=time.index)
        except ImportError:
            raise ImportError(
                "statsmodels is required for LOWESS fitting. "
                "Install with: pip install statsmodels"
            )
    
    else:
        raise ValueError(f"Unknown baseline method: {method}. Choose from 'lowess', 'constant', or 'polynomial'")
    
    return baseline_fit


def fit_well_baselines(
    df: pd.DataFrame,
    baseline_end_time: float = 24.0,
    method: str = 'lowess',
    frac: float = 0.5,
    poly_order: int = 1
) -> Dict[str, pd.Series]:
    """
    For each well, fit a smooth baseline function on the baseline window ONLY,
    then extrapolate over the full timecourse.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Long-format DataFrame with columns ['plate_id', 'well', 'content', 'time_s', 'value']
    baseline_end_time : float
        Maximum time (in seconds) for baseline window (default: 24.0)
    method : str
        Baseline fitting method: 'lowess', 'constant', or 'polynomial' (default: 'constant')
    frac : float
        Fraction of data used for LOWESS smoothing (default: 0.5)
    poly_order : int
        Polynomial order for polynomial fit (default: 1)
    
    Returns:
    --------
    Dict[str, pd.Series]
        Dictionary mapping well_id to Series of fitted baseline values (indexed by time_s)
    """
    baseline_fits = {}
    
    # Get baseline windows for each well
    baseline_windows = identify_baseline_window(df, baseline_end_time)
    
    # Fit baseline for each well
    for well_id, well_df in df.groupby('well'):
        # Get full timecourse for this well
        full_timecourse = well_df.sort_values('time_s')
        
        # Get baseline window
        baseline_window = baseline_windows.get(well_id)
        
        if baseline_window is None or len(baseline_window) == 0:
            # No baseline data, create NaN series
            baseline_fits[well_id] = pd.Series(
                index=full_timecourse['time_s'],
                dtype=float
            )
            continue
        
        # Sort baseline window by time
        baseline_window = baseline_window.sort_values('time_s')
        
        # Fit baseline on baseline window only
        baseline_fit = fit_baseline(
            baseline_window['time_s'],
            baseline_window['value'],
            method=method,
            frac=frac,
            poly_order=poly_order,
            well_id=well_id
        )
        
        # Extrapolate to full timecourse
        # Interpolate baseline fit to all timepoints in full timecourse
        from scipy.interpolate import interp1d
        
        # Get unique timepoints from baseline fit
        baseline_times = baseline_fit.index.values
        baseline_values = baseline_fit.values
        
        # Create interpolation function
        interp_func = interp1d(
            baseline_times,
            baseline_values,
            kind='linear',
            bounds_error=False,
            fill_value='extrapolate'
        )
        
        # Interpolate to full timecourse timepoints
        full_times = full_timecourse['time_s'].values
        extrapolated_baseline = pd.Series(
            interp_func(full_times),
            index=full_timecourse['time_s']
        )
        
        baseline_fits[well_id] = extrapolated_baseline
    
    return baseline_fits


def normalize_wells(
    df: pd.DataFrame,
    baseline_fits: Dict[str, pd.Series],
    normalization_mode: str = 'multiplicative'
) -> pd.DataFrame:
    """
    For each well and each timepoint, compute normalized fluorescence using the fitted baseline.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Long-format DataFrame with columns ['plate_id', 'well', 'content', 'time_s', 'value']
    baseline_fits : Dict[str, pd.Series]
        Dictionary mapping well_id to Series of fitted baseline values
    normalization_mode : str
        Normalization mode (default: 'multiplicative')
        - 'delta_f_over_f': ΔF/F normalization. Formula: F'(t) = (F(t) - g(t)) / g(t)
          Use when you want to measure relative change from baseline. Best for detecting
          increases/decreases relative to baseline. Values can be negative (below baseline)
          or positive (above baseline).
        - 'multiplicative': Multiplicative normalization. Formula: F_norm(t) = F(t) / g(t)
          Use when you want to normalize by baseline without subtracting. Best for fold-change
          analysis. Values are always positive (ratio to baseline).
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with same shape as input but with normalized 'value' column
    """
    norm_df = df.copy()
    norm_values = []
    
    for idx, row in df.iterrows():
        well_id = row['well']
        time_s = row['time_s']
        value = row['value']
        
        # Get baseline fit for this well
        baseline_series = baseline_fits.get(well_id)
        
        if baseline_series is None or len(baseline_series) == 0:
            norm_values.append(np.nan)
            continue
        
        # Find closest timepoint in baseline series
        baseline_times = baseline_series.index.values
        closest_idx = np.argmin(np.abs(baseline_times - time_s))
        baseline_value = baseline_series.iloc[closest_idx]
        
        if pd.isna(baseline_value) or baseline_value == 0:
            norm_values.append(np.nan)
            continue
        
        # Apply normalization
        if normalization_mode == 'delta_f_over_f':
            # ΔF/F = (F(t) - g(t)) / g(t)
            norm_value = (value - baseline_value) / baseline_value
        elif normalization_mode == 'multiplicative':
            # F_norm(t) = F(t) / g(t)
            norm_value = value / baseline_value
        else:
            raise ValueError(
                f"Unknown normalization_mode: {normalization_mode}. "
                "Choose from 'delta_f_over_f' or 'multiplicative'"
            )
        
        norm_values.append(norm_value)
    
    norm_df['value'] = norm_values
    return norm_df


def get_triplicate_group_color(
    group_name: str,
    config: Dict[str, Any] = None
) -> str:
    """
    Get the color for a triplicate group name.
    Matches the color assignment from the web interface.
    
    Parameters:
    -----------
    group_name : str
        Name of the triplicate group
    config : Dict[str, Any], optional
        Configuration dictionary with 'triplicate_groups' key
    
    Returns:
    --------
    str
        Hex color code for the group
    """
    if config is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "plate_config.json")
        config = load_config(config_path)
    
    triplicate_groups = config.get("triplicate_groups", DEFAULT_TRIPLICATE_GROUPS)
    
    # Build color map: group name -> color index
    color_map = {}
    for idx, group in enumerate(triplicate_groups):
        name = group.get("name", "").strip()
        if name and name not in color_map:
            color_map[name] = idx
    
    # Get color index for this group
    normalized_name = group_name.strip() if group_name else ""
    color_idx = color_map.get(normalized_name, 0)
    
    # Return color from palette (cycle if needed)
    return TRIPLICATE_COLORS[color_idx % len(TRIPLICATE_COLORS)]


def build_condition_map(
    df: pd.DataFrame,
    config: Dict[str, Any] = None
) -> Dict[str, str]:
    """
    Build a mapping from well_id to condition name.
    
    Uses triplicate group names from config to assign conditions.
    Wells in the same triplicate group (same row pattern) get the same condition.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Long-format DataFrame with columns ['plate_id', 'well', 'content', 'time_s', 'value']
    config : Dict[str, Any]
        Configuration dictionary with 'triplicate_groups' key
    
    Returns:
    --------
    Dict[str, str]
        Dictionary mapping well_id to condition name
    """
    if config is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "plate_config.json")
        config = load_config(config_path)
    
    condition_map = {}
    triplicate_groups = config.get("triplicate_groups", DEFAULT_TRIPLICATE_GROUPS)
    
    # Build mapping: row_letter -> condition name
    row_to_condition = {}
    for group in triplicate_groups:
        group_name = group.get("name", "")
        group_wells = group.get("wells", [])
        for well_pattern in group_wells:
            # Extract row letter from well pattern (e.g., "B" from "B" or "B13")
            row_letter = well_pattern[0] if well_pattern and well_pattern[0].isalpha() else ""
            if row_letter:
                row_to_condition[row_letter] = group_name
    
    # Map each well to its condition based on row letter
    for well_id in df['well'].unique():
        if well_id and len(well_id) > 0:
            row_letter = well_id[0] if well_id[0].isalpha() else ""
            condition = row_to_condition.get(row_letter, well_id)  # Default to well_id if not found
            condition_map[well_id] = condition
    
    return condition_map


def compute_condition_statistics(
    norm_df: pd.DataFrame,
    condition_map: Dict[str, str]
) -> pd.DataFrame:
    """
    Group wells by condition and compute mean ± SEM for each condition and timepoint.
    
    Parameters:
    -----------
    norm_df : pd.DataFrame
        Normalized DataFrame with columns ['plate_id', 'well', 'content', 'time_s', 'value']
    condition_map : Dict[str, str]
        Dictionary mapping well_id to condition name
    
    Returns:
    --------
    pd.DataFrame
        Tidy DataFrame with columns: ['Time_s', 'condition', 'mean', 'sd', 'sem', 'n']
    """
    import numpy as np
    
    # Add condition column
    norm_df = norm_df.copy()
    norm_df['condition'] = norm_df['well'].map(condition_map)
    
    # Group by condition and time
    results = []
    
    for (condition, time_s), group in norm_df.groupby(['condition', 'time_s']):
        values = group['value'].dropna().values
        
        if len(values) == 0:
            continue
        
        n = len(values)
        mean = np.mean(values)
        sd = np.std(values, ddof=1)  # Sample standard deviation
        sem = sd / np.sqrt(n) if n > 0 else np.nan
        
        results.append({
            'Time_s': time_s,
            'condition': condition,
            'mean': mean,
            'sd': sd,
            'sem': sem,
            'n': n
        })
    
    stats_df = pd.DataFrame(results)
    return stats_df.sort_values(['condition', 'Time_s']).reset_index(drop=True)


def plot_timecourse(
    stats_df: pd.DataFrame,
    conditions: List[str] = None,
    figsize: tuple = (10, 6),
    save_path: str = None,
    config: Dict[str, Any] = None
):
    """
    Generate time-course plots with mean ± SEM for each condition.
    Uses triplicate group colors to match the web interface.
    
    Parameters:
    -----------
    stats_df : pd.DataFrame
        Statistics DataFrame with columns ['Time_s', 'condition', 'mean', 'sem']
    conditions : List[str], optional
        Subset of conditions to plot. If None, plots all conditions.
    figsize : tuple
        Figure size (width, height) in inches (default: (10, 6))
    save_path : str, optional
        Path to save the figure. If None, displays the plot.
    config : Dict[str, Any], optional
        Configuration dictionary with 'triplicate_groups' key
    """
    import matplotlib.pyplot as plt
    
    if config is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "plate_config.json")
        config = load_config(config_path)
    
    if conditions is None:
        conditions = stats_df['condition'].unique()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    for idx, condition in enumerate(conditions):
        cond_data = stats_df[stats_df['condition'] == condition].sort_values('Time_s')
        
        if len(cond_data) == 0:
            continue
        
        times = cond_data['Time_s'].values
        means = cond_data['mean'].values
        sems = cond_data['sem'].values
        
        # Get color for this condition (triplicate group name)
        # Use triplicate group color if condition matches a group name, otherwise use default palette
        color = get_triplicate_group_color(condition, config)
        
        # Plot line
        ax.plot(times, means, label=condition, color=color, linewidth=2)
        
        # Plot error bars
        ax.errorbar(
            times, means, yerr=sems,
            color=color,
            alpha=0.5,
            capsize=3,
            capthick=1,
            elinewidth=1
        )
    
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('ΔF/F', fontsize=12)
    ax.set_title('Time Course by Condition', fontsize=14)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    else:
        plt.show()
    
    plt.close()


def print_baseline_options():
    """
    Print descriptions of baseline fitting and normalization options.
    """
    print("=" * 70)
    print("BASELINE FITTING OPTIONS")
    print("=" * 70)
    print("\n1. LOWESS (Locally Weighted Scatterplot Smoothing)")
    print("   Formula: g(t) = LOWESS(F_baseline, t_baseline, frac)")
    print("   When to use: Baseline has gradual drift or non-linear trends")
    print("   Best for: Non-linear baselines with smooth trends")
    print("   Parameters:")
    print("     - frac: Fraction of data used for smoothing (default: 0.5)")
    print("       * Higher frac (0.7-0.9) = smoother, more global fit")
    print("       * Lower frac (0.2-0.4) = more local, follows data closely")
    print()
    print("2. CONSTANT (0th-order polynomial)")
    print("   Formula: g(t) = mean(F_baseline) for all t")
    print("   When to use: Baseline is stable/flat with minimal drift")
    print("   Best for: Stable baselines with no time-dependent trends")
    print("   Parameters: None")
    print()
    print("3. POLYNOMIAL (1st-order or higher)")
    print("   Formula: g(t) = c₀ + c₁t + c₂t² + ...")
    print("   When to use: Baseline has linear or polynomial drift")
    print("   Best for: Baselines with known polynomial trends")
    print("   Parameters:")
    print("     - poly_order: Polynomial order (default: 1 = linear)")
    print("       * Order 1: Linear drift (g(t) = c₀ + c₁t)")
    print("       * Order 2: Quadratic drift (g(t) = c₀ + c₁t + c₂t²)")
    print("       * Higher orders: More complex trends")
    print()
    print("=" * 70)
    print("NORMALIZATION OPTIONS")
    print("=" * 70)
    print("\n1. DELTA_F_OVER_F (ΔF/F)")
    print("   Formula: F'(t) = (F(t) - g(t)) / g(t)")
    print("   When to use: Measure relative change from baseline")
    print("   Best for: Detecting increases/decreases relative to baseline")
    print("   Interpretation:")
    print("     * F'(t) = 0: No change from baseline")
    print("     * F'(t) > 0: Increase above baseline (e.g., 0.5 = 50% increase)")
    print("     * F'(t) < 0: Decrease below baseline (e.g., -0.3 = 30% decrease)")
    print("   Values can be negative (below baseline) or positive (above baseline)")
    print()
    print("2. MULTIPLICATIVE")
    print("   Formula: F_norm(t) = F(t) / g(t)")
    print("   When to use: Normalize by baseline without subtracting")
    print("   Best for: Fold-change analysis")
    print("   Interpretation:")
    print("     * F_norm(t) = 1.0: No change from baseline")
    print("     * F_norm(t) = 2.0: 2-fold increase (100% increase)")
    print("     * F_norm(t) = 0.5: 2-fold decrease (50% of baseline)")
    print("   Values are always positive (ratio to baseline)")
    print()
    print("=" * 70)
    print("STATISTICS COMPUTATION")
    print("=" * 70)
    print("\nFor each condition and timepoint, the following statistics are computed:")
    print("  Mean: μ(t) = (1/n) * Σ F'(t)")
    print("    Average normalized value across all wells in the condition")
    print()
    print("  SD: σ(t) = sqrt(Σ(F'(t) - μ(t))² / (n-1))")
    print("    Sample standard deviation (ddof=1) across wells")
    print()
    print("  SEM: SEM(t) = σ(t) / sqrt(n)")
    print("    Standard error of the mean = SD / sqrt(n_wells)")
    print("    Note: SEM is computed across wells, not propagated from baseline error")
    print()
    print("=" * 70)
    print()


def process_baseline_normalization(
    df: pd.DataFrame,
    baseline_end_time: float = 24.0,
    baseline_method: str = 'constant',
    baseline_frac: float = 0.5,
    baseline_poly_order: int = 1,
    normalization_mode: str = 'multiplicative',
    config: Dict[str, Any] = None,
    output_dir: str = None
) -> tuple:
    """
    Complete pipeline: baseline identification, fitting, normalization, and statistics.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Long-format DataFrame with columns ['plate_id', 'well', 'content', 'time_s', 'value']
    baseline_end_time : float
        Maximum time (in seconds) for baseline window (default: 24.0)
    baseline_method : str
        Baseline fitting method: 'lowess', 'constant', or 'polynomial' (default: 'constant')
    baseline_frac : float
        Fraction of data used for LOWESS smoothing (default: 0.5)
    baseline_poly_order : int
        Polynomial order for polynomial fit (default: 1)
    normalization_mode : str
        Normalization mode: 'delta_f_over_f' or 'multiplicative' (default: 'multiplicative')
    config : Dict[str, Any], optional
        Configuration dictionary. If None, loads from plate_config.json
    output_dir : str, optional
        Directory to save outputs. If None, uses current directory.
    
    Returns:
    --------
    tuple
        (norm_df, stats_df, baseline_fits)
        - norm_df: Normalized DataFrame
        - stats_df: Condition statistics DataFrame
        - baseline_fits: Dictionary of baseline fits per well
    """
    print("=" * 70)
    print("BASELINE IDENTIFICATION AND NORMALIZATION PIPELINE")
    print("=" * 70)
    print(f"\nCurrent Settings:")
    print(f"  Baseline window: Time_s <= {baseline_end_time} s")
    print(f"  Baseline method: {baseline_method}")
    if baseline_method == 'lowess':
        print(f"    LOWESS smoothing fraction (frac): {baseline_frac}")
    elif baseline_method == 'polynomial':
        print(f"    Polynomial order: {baseline_poly_order}")
    print(f"  Normalization mode: {normalization_mode}")
    print()
    print("For detailed option descriptions, call: print_baseline_options()")
    print("=" * 70)
    print()
    
    print("Step 1: Identifying baseline windows...", flush=True)
    baseline_windows = identify_baseline_window(df, baseline_end_time)
    print(f"  Identified baseline windows for {len(baseline_windows)} wells", flush=True)
    total_baseline_points = sum(len(bw) for bw in baseline_windows.values())
    print(f"  Total baseline points: {total_baseline_points}", flush=True)
    
    print("Step 2: Fitting baselines...", flush=True)
    baseline_fits = fit_well_baselines(
        df,
        baseline_end_time=baseline_end_time,
        method=baseline_method,
        frac=baseline_frac,
        poly_order=baseline_poly_order
    )
    print(f"  Fitted baselines for {len(baseline_fits)} wells", flush=True)
    
    print("Step 3: Normalizing wells...", flush=True)
    norm_df = normalize_wells(df, baseline_fits, normalization_mode=normalization_mode)
    print(f"  Normalized {len(norm_df)} data points", flush=True)
    
    print("Step 4: Building condition map...", flush=True)
    condition_map = build_condition_map(df, config=config)
    print(f"  Mapped {len(condition_map)} wells to conditions", flush=True)
    condition_counts = {}
    for well_id, condition in condition_map.items():
        condition_counts[condition] = condition_counts.get(condition, 0) + 1
    print(f"  Condition distribution:", flush=True)
    for condition, count in sorted(condition_counts.items()):
        print(f"    {condition}: {count} wells", flush=True)
    
    print("Step 5: Computing condition statistics...", flush=True)
    stats_df = compute_condition_statistics(norm_df, condition_map)
    print(f"  Computed statistics for {len(stats_df)} condition-timepoint combinations", flush=True)
    
    # Save outputs if output_dir is provided
    if output_dir:
        ensure_dir(output_dir)
        
        print("Step 6: Saving outputs...", flush=True)
        
        # Save normalized data
        norm_csv_path = os.path.join(output_dir, "normalized_data.csv")
        norm_df.to_csv(norm_csv_path, index=False)
        print(f"  Saved normalized data to {norm_csv_path}", flush=True)
        
        # Save condition statistics
        stats_csv_path = os.path.join(output_dir, "condition_statistics.csv")
        stats_df.to_csv(stats_csv_path, index=False)
        print(f"  Saved condition statistics to {stats_csv_path}", flush=True)
        
        # Save plot
        plot_path = os.path.join(output_dir, "timecourse_plot.png")
        plot_timecourse(stats_df, save_path=plot_path, config=config)
        print(f"  Saved plot to {plot_path}", flush=True)
    
    return norm_df, stats_df, baseline_fits


def main():
    # Use the directory where this script lives as the working directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    # Find Excel files in this directory
    excel_files = [
        f for f in os.listdir(".")
        if f.lower().endswith(".xlsx") and not f.startswith("~$")
    ]

    if not excel_files:
        print("No .xlsx files found in this folder.", file=sys.stderr)
        sys.exit(1)
    
    # Validate filename formats before processing
    print("Starting plate viewer script...", flush=True)
    print("Validating filename formats...", flush=True)
    invalid_files = []
    for fname in excel_files:
        try:
            validate_filename_format(fname)
        except ValueError as e:
            invalid_files.append((fname, str(e)))
    
    if invalid_files:
        print("\n❌ ERROR: Invalid filename format(s) detected:", file=sys.stderr)
        for fname, error_msg in invalid_files:
            print(f"\nFile: {fname}", file=sys.stderr)
            print(error_msg, file=sys.stderr)
        print("\nPlease rename the file(s) to match the expected format and rerun the script.", file=sys.stderr)
        sys.exit(1)
    
    print(f"✓ All {len(excel_files)} file(s) have valid filename formats", flush=True)

    all_plates: Dict[str, pd.DataFrame] = {}
    csv_root = os.path.join(script_dir, "csv")
    
    # Load configuration file early so it can be used for CSV generation
    config_path = os.path.join(script_dir, "plate_config.json")
    config = load_config(config_path)
    
    # Track plate_ids and their source files to detect duplicates
    plate_id_to_files: Dict[str, List[str]] = {}
    # Track plate_id to original filename mapping
    plate_id_to_filename: Dict[str, str] = {}

    for idx, fname in enumerate(excel_files, 1):
        print(f"\n[{idx}/{len(excel_files)}] Processing {fname}...", flush=True)
        path = os.path.join(script_dir, fname)
        try:
            long_df = parse_plate_file(path)
        except Exception as e:
            print(f"Error parsing {fname}: {e}", file=sys.stderr, flush=True)
            continue

        plate_id = long_df["plate_id"].iloc[0]
        
        # Track which files produce which plate_id
        if plate_id not in plate_id_to_files:
            plate_id_to_files[plate_id] = []
        plate_id_to_files[plate_id].append(fname)
        
        # Track filename for this plate_id (use first occurrence)
        if plate_id not in plate_id_to_filename:
            plate_id_to_filename[plate_id] = fname
        
        # If this plate_id already exists, make it unique by appending filename
        if plate_id in all_plates:
            # This is a duplicate - we'll handle it below
            pass
        else:
            all_plates[plate_id] = long_df
            write_csvs_for_plate(long_df, csv_root, config, fname)

    # Check for duplicates and warn
    duplicates_found = False
    for plate_id, files in plate_id_to_files.items():
        if len(files) > 1:
            duplicates_found = True
            print(f"\n⚠️  WARNING: Duplicate dataset detected!", file=sys.stderr)
            print(f"   Plate ID '{plate_id}' is produced by multiple files:", file=sys.stderr)
            for f in files:
                print(f"     - {f}", file=sys.stderr)
            print(f"   Please delete the duplicate file(s) and rerun the script.", file=sys.stderr)
    
    if duplicates_found:
        print("\n⚠️  ERROR: Duplicate datasets found. Please delete duplicates and rerun.", file=sys.stderr)
        print("   Only the first occurrence of each duplicate will be included.", file=sys.stderr)
        sys.exit(1)

    if not all_plates:
        print("No plates were successfully parsed.", file=sys.stderr)
        sys.exit(1)
    
    print(f"\n✓ Successfully loaded {len(all_plates)} dataset(s):", flush=True)
    for plate_id in sorted(all_plates.keys()):
        print(f"   - {plate_id}", flush=True)

    # Config already loaded above

    # Build JSON for viewer
    web_dir = os.path.join(script_dir, "web")
    ensure_dir(web_dir)
    json_path = os.path.join(web_dir, "viewer_data.json")
    build_viewer_json(all_plates, json_path, config, plate_id_to_filename)

    # Write HTML
    html_path = os.path.join(web_dir, "index.html")
    write_html(html_path)

    # Create double-clickable web.command file
    web_command_path = write_web_command(script_dir, web_dir)

    print("Done.", flush=True)
    print("Generated:", flush=True)
    print(f"  CSVs in: {csv_root}", flush=True)
    print(f"  Web viewer: {html_path}", flush=True)
    print(f"  Web launcher: {web_command_path}", flush=True)
    print("\nTo view the plates:", flush=True)
    print(f"  Double-click 'web.command' in Finder, or", flush=True)
    print(f"  Run: ./web.command", flush=True)
    print("The browser will open automatically.", flush=True)


if __name__ == "__main__":
    main()
