"""
Parsers for Delphi Lab XPonent and Accession files.
"""

import pandas as pd
import io
import re
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, field


@dataclass
class XPonentData:
    """Parsed XPonent file data."""
    net_mfi: pd.DataFrame  # Sample x Gene matrix of Net MFI values
    counts: pd.DataFrame   # Sample x Gene matrix of bead counts
    metadata: Dict[str, str]  # File metadata (batch, date, etc.)
    samples: List[str]  # List of sample IDs
    negative_controls: List[str]  # List of NC sample IDs


@dataclass
class AccessionData:
    """Parsed Accession file data."""
    samples: pd.DataFrame  # Sample info with staging data
    sample_id_column: str  # Name of sample ID column


def parse_xponent_csv(file_content: bytes) -> XPonentData:
    """
    Parse XPonent CSV file exported from Luminex instrument.

    The file has multiple sections:
    - Header metadata
    - DataType: Median
    - DataType: Net MFI (what we need)
    - DataType: Count (for QC)
    """
    # Decode content
    try:
        content = file_content.decode('utf-8')
    except UnicodeDecodeError:
        content = file_content.decode('latin-1')

    lines = content.strip().split('\n')

    # Parse metadata from header
    metadata = {}
    for line in lines[:30]:
        if ',' in line:
            parts = line.split(',')
            if len(parts) >= 2:
                key = parts[0].strip().strip('"')
                val = parts[1].strip().strip('"')
                if key and val and key not in ['DataType:', 'Location']:
                    metadata[key] = val

    # Find section boundaries using flexible pattern matching
    net_mfi_start = None
    count_start = None

    for i, line in enumerate(lines):
        line_lower = line.lower().replace('"', '').replace("'", "")
        # Match variations like "DataType:,Net MFI", "DataType:, Net MFI", etc.
        if re.search(r'datatype:\s*,?\s*net\s*mfi', line_lower):
            net_mfi_start = i
        elif re.search(r'datatype:\s*,?\s*count', line_lower):
            count_start = i

    if net_mfi_start is None:
        raise ValueError("Could not find 'Net MFI' section in XPonent file. "
                        "Make sure you're uploading the correct file format. "
                        "The file should contain a section marked 'DataType:,Net MFI'.")

    # Parse Net MFI section
    net_mfi_df = _parse_xponent_section(lines, net_mfi_start)

    # Parse Count section if available
    if count_start is not None:
        count_df = _parse_xponent_section(lines, count_start)
    else:
        # Create dummy counts with all passing values
        count_df = net_mfi_df.copy()
        for col in count_df.columns:
            if col != 'Sample':
                count_df[col] = 100

    # Identify negative controls using flexible pattern matching
    negative_controls = []
    nc_patterns = [
        r'^nc$',                          # Exact "nc"
        r'^nc[-_\s]?\d*$',                # nc, nc1, nc-1, nc_1, nc 1
        r'^nc[-_\s]',                     # nc-anything, nc_buffer, etc.
        r'^neg(ative)?[-_\s]?(ctrl|control)?$',  # neg, negative, negctrl, neg control
        r'^negative[-_\s]?control',       # negative control, negative-control
        r'^background',                   # background control
        r'^blank',                        # blank control
    ]

    for sample in net_mfi_df['Sample'].tolist():
        sample_lower = str(sample).lower().strip()
        if any(re.match(p, sample_lower) for p in nc_patterns):
            negative_controls.append(sample)

    return XPonentData(
        net_mfi=net_mfi_df,
        counts=count_df,
        metadata=metadata,
        samples=net_mfi_df['Sample'].tolist(),
        negative_controls=negative_controls
    )


def _parse_xponent_section(lines: List[str], start_idx: int) -> pd.DataFrame:
    """Parse a data section from XPonent file."""
    # Find header row (Location,Sample,gene1,gene2,...)
    header_idx = None
    for i in range(start_idx + 1, min(start_idx + 5, len(lines))):
        line = lines[i]
        if 'Location' in line and 'Sample' in line:
            header_idx = i
            break
        elif 'Sample' in line and ('SLC39A6' in line or 'ESR1' in line):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"Could not find header row in section starting at line {start_idx}")

    # Read data from header to next empty line or next DataType
    data_lines = [lines[header_idx]]
    for i in range(header_idx + 1, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith('"DataType:') or line.startswith('DataType:'):
            break
        if line.startswith('"Avg') or line.startswith('Avg'):
            break
        data_lines.append(line)

    # Parse as CSV
    csv_content = '\n'.join(data_lines)
    df = pd.read_csv(io.StringIO(csv_content))

    # Clean up column names
    df.columns = [str(c).strip().strip('"') for c in df.columns]

    # Drop Location column if present
    if 'Location' in df.columns:
        df = df.drop(columns=['Location'])

    # Drop Total Events if present (we'll keep it separate)
    if 'Total Events' in df.columns:
        df = df.drop(columns=['Total Events'])

    # Convert Sample column to string
    df['Sample'] = df['Sample'].astype(str).str.strip()

    # Convert gene columns to numeric
    for col in df.columns:
        if col != 'Sample':
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

    return df


def parse_accession_xlsx(file_content: bytes) -> AccessionData:
    """
    Parse Accession XLSX file with specimen information.

    The file has:
    - Header rows with form info
    - Column headers row with: Sample Id, Clinical or Pathology, Tumor Stage, etc.
    - Data rows
    """
    df = pd.read_excel(io.BytesIO(file_content), header=None)

    # Find the header row by looking for multiple expected column names
    header_row = None
    sample_id_col = None

    # Keywords that indicate a header row
    header_keywords = ['sample id', 'sample_id', 'sampleid', 'clinical', 'patholog',
                       'tumor stage', 'tumor node', 'tumor size', 'accession']

    for idx, row in df.iterrows():
        # Count how many header-like values are in this row
        keyword_matches = 0
        has_sample_id = False

        for val in row.values:
            if pd.notna(val):
                val_lower = str(val).lower()
                if any(kw in val_lower for kw in header_keywords):
                    keyword_matches += 1
                if 'sample id' in val_lower or 'sample_id' in val_lower or 'sampleid' in val_lower:
                    has_sample_id = True

        # Need at least 2 keyword matches AND must have sample id
        if keyword_matches >= 2 and has_sample_id:
            header_row = idx
            # Find which column has sample id
            for col_idx, val in enumerate(row.values):
                if pd.notna(val):
                    val_str = str(val).lower()
                    if 'sample id' in val_str or 'sample_id' in val_str or 'sampleid' in val_str:
                        sample_id_col = col_idx
                        break
            break

    if header_row is None:
        raise ValueError("Could not find header row with 'Sample Id' column in Accession file. "
                        "Make sure your file has columns like 'Sample Id', 'Clinical or Pathology', "
                        "'Tumor Stage', 'Tumor Node', etc.")

    # Set the header row
    df.columns = df.iloc[header_row]
    df = df.iloc[header_row + 1:].reset_index(drop=True)

    # Clean column names
    column_mapping = {}
    for col in df.columns:
        col_str = str(col).strip()
        col_lower = col_str.lower()

        # Normalize column names
        if 'sample id' in col_lower or 'sample_id' in col_lower:
            column_mapping[col] = 'sample_id'
        elif 'clinical or path' in col_lower:
            column_mapping[col] = 'sample_type'
        elif 'tumor size' in col_lower:
            column_mapping[col] = 'tumor_size_mm'
        elif 'tumor stage' in col_lower:
            column_mapping[col] = 'tumor_stage'
        elif 'positive lymph' in col_lower or 'number of positive' in col_lower:
            column_mapping[col] = 'positive_lymph_nodes'
        elif 'tumor node' in col_lower:
            column_mapping[col] = 'tumor_node'
        elif 'accession' in col_lower or 'case number' in col_lower:
            column_mapping[col] = 'accession_number'

    df = df.rename(columns=column_mapping)

    # Clean up sample_id
    if 'sample_id' in df.columns:
        df['sample_id'] = df['sample_id'].astype(str).str.strip()
        # Remove rows where sample_id is empty or NaN
        df = df[df['sample_id'].notna() & (df['sample_id'] != '') & (df['sample_id'] != 'nan')]

    # Clean up sample_type
    if 'sample_type' in df.columns:
        df['sample_type'] = df['sample_type'].astype(str).str.strip()
        # Normalize to Clinical/Pathological
        df['sample_type'] = df['sample_type'].apply(_normalize_sample_type)

    # Clean up staging columns
    if 'tumor_stage' in df.columns:
        df['tumor_stage'] = df['tumor_stage'].astype(str).str.strip()

    if 'tumor_node' in df.columns:
        df['tumor_node'] = df['tumor_node'].astype(str).str.strip()

    if 'tumor_size_mm' in df.columns:
        df['tumor_size_mm'] = pd.to_numeric(df['tumor_size_mm'], errors='coerce')

    if 'positive_lymph_nodes' in df.columns:
        df['positive_lymph_nodes'] = pd.to_numeric(df['positive_lymph_nodes'], errors='coerce')

    return AccessionData(
        samples=df,
        sample_id_column='sample_id'
    )


def _normalize_sample_type(val: str) -> str:
    """Normalize sample type to Clinical or Pathological."""
    val_lower = val.lower().strip()
    if val_lower in ['pathological', 'pathologic', 'pathology', 'path', 'p']:
        return 'Pathologic'
    elif val_lower in ['clinical', 'clin', 'c']:
        return 'Clinical'
    elif val_lower in ['control', 'ctrl']:
        return 'Control'
    return val


def merge_xponent_accession(xponent: XPonentData, accession: AccessionData) -> pd.DataFrame:
    """
    Merge XPonent gene expression data with Accession sample info.

    Returns a DataFrame with all gene values and staging info per sample.
    """
    # Get accession data
    acc_df = accession.samples.copy()

    # Create merged dataframe
    merged_rows = []

    for _, row in xponent.net_mfi.iterrows():
        sample_id = str(row['Sample']).strip()

        # Find matching accession record
        acc_match = acc_df[acc_df['sample_id'].str.strip() == sample_id]

        if len(acc_match) > 0:
            acc_row = acc_match.iloc[0]
            merged_row = {
                'sample_id': sample_id,
                'sample_type': acc_row.get('sample_type', 'Clinical'),
                'tumor_stage': acc_row.get('tumor_stage', ''),
                'tumor_node': acc_row.get('tumor_node', ''),
                'tumor_size_mm': acc_row.get('tumor_size_mm', 0),
                'positive_lymph_nodes': acc_row.get('positive_lymph_nodes', 0),
            }
        else:
            # No accession match - use defaults
            merged_row = {
                'sample_id': sample_id,
                'sample_type': 'Clinical',
                'tumor_stage': '',
                'tumor_node': '',
                'tumor_size_mm': 0,
                'positive_lymph_nodes': 0,
            }

        # Add gene values
        for col in xponent.net_mfi.columns:
            if col != 'Sample':
                merged_row[col] = row[col]

        merged_rows.append(merged_row)

    return pd.DataFrame(merged_rows)
