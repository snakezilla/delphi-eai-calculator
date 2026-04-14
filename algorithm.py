"""
Delphi EAI Algorithm v2.0 - SET2,3 Score Calculator

This implements the complete SET2,3 algorithm for HR+/HER2- breast cancer
risk assessment, supporting both Clinical and Pathological sample types.
"""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from enum import Enum
import pandas as pd

ALGO_VERSION = "2.0.0"

# Gene classifications
REFERENCE_GENES = [
    'ATP5J2', 'VDAC2', 'DARS', 'UGP2', 'UBE2Z',
    'AK2', 'WIPF2', 'APPBP2', 'TRIM2', 'LDHA'
]

TARGET_GENES = [
    'SLC39A6', 'STC2', 'CA12', 'ESR1', 'PDZK1', 'NPY1R', 'CD2', 'MAPT',
    'QDPR', 'AZGP1', 'ABAT', 'ADCY1', 'CD3D', 'NAT1', 'MRPS30', 'DNAJC12',
    'SCUBE2', 'KCNE4'
]

RNA4_GENES = ['ESR1', 'ERBB2', 'PGR', 'AURKA']

ALL_GENES = list(set(REFERENCE_GENES + TARGET_GENES + RNA4_GENES))

# Calibration coefficients
CALIBRATION = {
    'SET_ERPR_SLOPE': 0.7548,
    'SET_ERPR_INTERCEPT': 0.2912,
    'ESR1_SLOPE': 0.7178,
    'ESR1_INTERCEPT': 4.0114,
    'ERBB2_SLOPE': 0.7738,
    'ERBB2_INTERCEPT': 1.8455,
    'PGR_SLOPE': 0.768,
    'PGR_INTERCEPT': 0.854,
    'AURKA_SLOPE': 0.8991,
    'AURKA_INTERCEPT': -0.1311
}

# Risk thresholds
THRESHOLDS = {
    'ESR1': 8.93,
    'ERBB2': 11.97,
    'PGR_HIGH': 5.1,
    'PGR_BORDERLINE': 4.5,
    'SET23_CLINICAL': 1.95,
    'SET23_PATHOLOGICAL': 2.1
}

# QC thresholds
QC_MIN_BEAD_COUNT = 30
QC_MIN_MEAN_REFERENCE = 3.0
QC_MIN_BACKGROUND = 3.5


class PGRRiskLevel(Enum):
    HIGH = 0
    BORDERLINE = 1
    LOW = 2


@dataclass
class EAIResult:
    """Complete result from EAI algorithm calculation."""
    sample_id: str = ""
    sample_type: str = ""  # Clinical or Pathological

    # QC Results
    qc_passed: bool = True
    qc_bead_count_passed: bool = True
    qc_mean_ref_passed: bool = True
    qc_background_passed: bool = True
    qc_genes_present: bool = True
    qc_accession_matched: bool = True
    qc_fail_reason: str = ""
    min_bead_count: int = 0
    missing_genes: str = ""

    # Intermediate calculations
    avg_reference_genes: float = 0.0
    avg_target_genes: float = 0.0
    background: float = 0.0  # Mean of NC values

    # SET ER/PR
    set_erpr: float = 0.0
    calibrated_set_erpr: float = 0.0

    # Calibrated RNA4 values
    calibrated_esr1: float = 0.0
    calibrated_erbb2: float = 0.0
    calibrated_pgr: float = 0.0
    calibrated_aurka: float = 0.0

    # Risk scores
    esr1_risk: float = 0.0
    erbb2_risk: float = 0.0
    pgr_risk: float = 0.0
    pgr_risk_level: str = ""
    aurka_risk: float = 0.0
    rna4_score: float = 0.0

    # Staging info
    tumor_stage: str = ""
    tumor_node: str = ""
    tumor_size_mm: float = 0.0
    positive_lymph_nodes: int = 0

    # Staging risk
    stage_risk: float = 0.0
    node_risk: float = 0.0

    # Final scores
    bpi: float = 0.0
    set23: float = 0.0
    set23_category: str = ""

    algo_version: str = ALGO_VERSION


def log2_transform(value: float) -> float:
    """Log2 transform with floor at 0."""
    if value <= 0:
        return 0.0
    result = math.log2(value)
    return max(0.0, result)


def calculate_mean(values: List[float]) -> float:
    """Calculate mean of non-NaN values. Returns NaN if no valid values."""
    valid = []
    for v in values:
        if v is not None:
            try:
                if not math.isnan(v):
                    valid.append(v)
            except (TypeError, ValueError):
                pass
    if not valid:
        return float('nan')  # Return NaN to signal invalid data
    return sum(valid) / len(valid)


def calculate_pgr_risk(calibrated_pgr: float) -> Tuple[float, PGRRiskLevel, str]:
    """Calculate PGR risk with 3 levels."""
    if calibrated_pgr > THRESHOLDS['PGR_HIGH']:
        return 0.0, PGRRiskLevel.HIGH, "High"
    elif calibrated_pgr > THRESHOLDS['PGR_BORDERLINE']:
        return 1.0, PGRRiskLevel.BORDERLINE, "Borderline"
    else:
        return 2.0, PGRRiskLevel.LOW, "Low"


def calculate_aurka_risk(calibrated_aurka: float, pgr_level: PGRRiskLevel) -> float:
    """Calculate AURKA risk based on PGR risk level."""
    if pgr_level in [PGRRiskLevel.HIGH, PGRRiskLevel.BORDERLINE]:
        if calibrated_aurka < 7.0:
            return 0.0
        elif calibrated_aurka > 9.0:
            return 2.0
        else:
            return calibrated_aurka - 7.0
    else:  # LOW
        if calibrated_aurka < 7.5:
            return 0.0
        elif calibrated_aurka > 8.5:
            return 1.0
        else:
            return calibrated_aurka - 7.5


def calculate_clinical_staging(tumor_stage: str, tumor_node: str) -> Tuple[float, float]:
    """Calculate clinical staging risk values."""
    stage_map = {'ct0': 0, 'ct1': 0, 'ct2': 1, 'ct3': 2, 'ct4': 3}
    node_map = {'cn0': 0, 'cn1': 1, 'cn2': 2, 'cn3': 3}

    tumor_stage_clean = str(tumor_stage).lower().strip()
    tumor_node_clean = str(tumor_node).lower().strip()

    stage_risk = stage_map.get(tumor_stage_clean, float('nan'))
    node_risk = node_map.get(tumor_node_clean, float('nan'))

    return stage_risk, node_risk


def calculate_pathological_staging(
    tumor_stage: str,
    tumor_node: str,
    tumor_size: float,
    positive_nodes: int
) -> Tuple[float, float]:
    """Calculate pathological staging risk values."""
    tumor_stage_clean = str(tumor_stage).lower().strip()
    tumor_node_clean = str(tumor_node).lower().strip()

    # Handle NaN values
    if pd.isna(tumor_size):
        tumor_size = 0.0
    if pd.isna(positive_nodes):
        positive_nodes = 0

    # pStageRisk
    if tumor_size < 10.0:
        stage_risk = 0.0
    elif tumor_size >= 40.0 or tumor_stage_clean == 'pt4':
        stage_risk = 3.0
    else:
        stage_risk = (tumor_size / 10.0) - 1.0

    # pNodeRisk
    if positive_nodes >= 10 or tumor_node_clean == 'pn3':
        node_risk = 5.0
    elif tumor_node_clean == 'pn0':
        node_risk = 0.0
    else:
        node_risk = 0.5 * positive_nodes

    return stage_risk, node_risk


def calculate_negative_control_adjustment(
    net_mfi_df: pd.DataFrame,
    nc_sample_ids: List[str]
) -> Dict[str, float]:
    """
    Calculate mean adjustment factor from negative controls.

    Returns dict of gene -> mean NC value.
    """
    if not nc_sample_ids:
        return {gene: 0.0 for gene in ALL_GENES}

    adjustment = {}
    nc_rows = net_mfi_df[net_mfi_df['Sample'].isin(nc_sample_ids)]

    for gene in ALL_GENES:
        if gene in nc_rows.columns:
            values = nc_rows[gene].tolist()
            adjustment[gene] = calculate_mean(values)
        else:
            adjustment[gene] = 0.0

    return adjustment


def calculate_eai(
    gene_mfi: Dict[str, float],
    bead_counts: Optional[Dict[str, int]],
    nc_adjustment: Dict[str, float],
    sample_id: str,
    sample_type: str,
    tumor_stage: str,
    tumor_node: str,
    tumor_size_mm: float,
    positive_lymph_nodes: int
) -> EAIResult:
    """
    Run the complete EAI algorithm for a single sample.

    Args:
        gene_mfi: Dict of gene name -> Net MFI value
        bead_counts: Dict of gene name -> bead count (for QC)
        nc_adjustment: Dict of gene name -> NC adjustment factor
        sample_id: Sample identifier
        sample_type: "Clinical" or "Pathological"
        tumor_stage: cT0-cT4 or pT0-pT4
        tumor_node: cN0-cN3 or pN0-pN3
        tumor_size_mm: Tumor size in mm (pathological only)
        positive_lymph_nodes: Number of positive nodes (pathological only)

    Returns:
        EAIResult with all calculated values
    """
    result = EAIResult(
        sample_id=sample_id,
        sample_type=sample_type,
        tumor_stage=tumor_stage,
        tumor_node=tumor_node,
        tumor_size_mm=tumor_size_mm if not pd.isna(tumor_size_mm) else 0.0,
        positive_lymph_nodes=int(positive_lymph_nodes) if not pd.isna(positive_lymph_nodes) else 0
    )

    is_pathological = sample_type.lower().strip() == "pathologic"

    # Check for missing critical genes
    missing_genes = []
    for gene in ALL_GENES:
        if gene not in gene_mfi or gene_mfi.get(gene) is None:
            missing_genes.append(gene)

    if missing_genes:
        result.qc_genes_present = False
        result.qc_passed = False
        result.missing_genes = ', '.join(missing_genes[:5])
        if len(missing_genes) > 5:
            result.missing_genes += f'... (+{len(missing_genes) - 5} more)'
        result.qc_fail_reason = f"Missing genes: {result.missing_genes}"

    # Step 1: Subtract adjustment factor (negative control)
    adjusted_mfi = {}
    for gene in ALL_GENES:
        raw = gene_mfi.get(gene, 0.0) or 0.0
        # Handle NaN in NC adjustment
        nc_val = nc_adjustment.get(gene, 0.0)
        nc = nc_val if nc_val is not None and not (isinstance(nc_val, float) and math.isnan(nc_val)) else 0.0
        adjusted_mfi[gene] = max(0, raw - nc)

    # Calculate background (mean of NC adjustments used)
    nc_values = [nc_adjustment.get(g, 0) for g in REFERENCE_GENES]
    result.background = calculate_mean(nc_values)

    # Step 1b-c: Log2 transform with floor at 0
    log2_values = {gene: log2_transform(val) for gene, val in adjusted_mfi.items()}

    # Step 2-3: Calculate means
    ref_vals = [log2_values.get(g, 0.0) for g in REFERENCE_GENES]
    target_vals = [log2_values.get(g, 0.0) for g in TARGET_GENES]

    ref_mean = calculate_mean(ref_vals)
    target_mean = calculate_mean(target_vals)

    result.avg_reference_genes = ref_mean
    result.avg_target_genes = target_mean

    # Step 4: QC checks
    # 4a: Bead count check
    if bead_counts:
        min_count = min(bead_counts.values()) if bead_counts else 100
        result.min_bead_count = min_count
        if min_count < QC_MIN_BEAD_COUNT:
            result.qc_bead_count_passed = False
            result.qc_passed = False
            result.qc_fail_reason = f"Bead count {min_count} < {QC_MIN_BEAD_COUNT}"
    else:
        result.qc_bead_count_passed = True
        result.min_bead_count = 100  # Assume passing if no counts

    # 4b: Mean reference genes check (also fails if ref_mean is NaN)
    if math.isnan(ref_mean) or ref_mean < QC_MIN_MEAN_REFERENCE:
        result.qc_mean_ref_passed = False
        result.qc_passed = False
        if result.qc_fail_reason:
            result.qc_fail_reason += f"; Mean ref genes {ref_mean:.2f} < {QC_MIN_MEAN_REFERENCE}"
        else:
            result.qc_fail_reason = f"Mean ref genes {ref_mean:.2f} < {QC_MIN_MEAN_REFERENCE}"

    # 4c: Background check (mean of NC adjustment values for reference genes)
    bg = result.background
    if not math.isnan(bg) and bg < QC_MIN_BACKGROUND:
        result.qc_background_passed = False
        result.qc_passed = False
        if result.qc_fail_reason:
            result.qc_fail_reason += f"; Background {bg:.2f} < {QC_MIN_BACKGROUND}"
        else:
            result.qc_fail_reason = f"Background {bg:.2f} < {QC_MIN_BACKGROUND}"

    # Step 5: If QC fails for controls, we still calculate but flag it
    # For actual samples, we continue calculation even if QC fails but report it

    # Step 6-8: SET ER/PR
    set_erpr = target_mean - ref_mean + 2
    calibrated_set_erpr = (set_erpr * CALIBRATION['SET_ERPR_SLOPE']) + CALIBRATION['SET_ERPR_INTERCEPT']
    calibrated_set_erpr = max(0.0, calibrated_set_erpr)

    result.set_erpr = set_erpr
    result.calibrated_set_erpr = calibrated_set_erpr

    # Step 9-13: RNA4 normalized and calibrated values
    norm_esr1 = log2_values.get('ESR1', 0.0) - ref_mean + 10
    norm_erbb2 = log2_values.get('ERBB2', 0.0) - ref_mean + 10
    norm_pgr = log2_values.get('PGR', 0.0) - ref_mean + 10
    norm_aurka = log2_values.get('AURKA', 0.0) - ref_mean + 10

    cal_esr1 = (norm_esr1 * CALIBRATION['ESR1_SLOPE']) + CALIBRATION['ESR1_INTERCEPT']
    cal_erbb2 = (norm_erbb2 * CALIBRATION['ERBB2_SLOPE']) + CALIBRATION['ERBB2_INTERCEPT']
    cal_pgr = (norm_pgr * CALIBRATION['PGR_SLOPE']) + CALIBRATION['PGR_INTERCEPT']
    cal_aurka = (norm_aurka * CALIBRATION['AURKA_SLOPE']) + CALIBRATION['AURKA_INTERCEPT']

    result.calibrated_esr1 = cal_esr1
    result.calibrated_erbb2 = cal_erbb2
    result.calibrated_pgr = cal_pgr
    result.calibrated_aurka = cal_aurka

    # Step 14: ESR1 Risk
    esr1_risk = 0.0 if cal_esr1 > THRESHOLDS['ESR1'] else 0.5
    result.esr1_risk = esr1_risk

    # Step 15: ERBB2 Risk (threshold 11.97)
    erbb2_risk = 0.5 if cal_erbb2 >= THRESHOLDS['ERBB2'] else 0.0
    result.erbb2_risk = erbb2_risk

    # Step 16: PGR Risk (3 levels)
    pgr_risk, pgr_level, pgr_level_str = calculate_pgr_risk(cal_pgr)
    result.pgr_risk = pgr_risk
    result.pgr_risk_level = pgr_level_str

    # Step 17: AURKA Risk (depends on PGR)
    aurka_risk = calculate_aurka_risk(cal_aurka, pgr_level)
    result.aurka_risk = aurka_risk

    # Step 18: RNA4 Score
    rna4_score = esr1_risk + erbb2_risk + pgr_risk + aurka_risk
    result.rna4_score = rna4_score

    # Steps 19-28: Staging and BPI
    if is_pathological:
        stage_risk, node_risk = calculate_pathological_staging(
            tumor_stage, tumor_node, tumor_size_mm, positive_lymph_nodes
        )
        # Pathological BPI: [11 - (pStageRisk + pNodeRisk + RNA4)] * 4/11
        bpi = (11 - (stage_risk + node_risk + rna4_score)) * (4.0 / 11.0)
        set23_threshold = THRESHOLDS['SET23_PATHOLOGICAL']
    else:
        stage_risk, node_risk = calculate_clinical_staging(tumor_stage, tumor_node)
        # Clinical BPI: [8 - (cStageRisk + cNodeRisk + RNA4*2/3)] * 3/8
        if not math.isnan(stage_risk) and not math.isnan(node_risk):
            bpi = (8 - (stage_risk + node_risk + (rna4_score * 2.0 / 3.0))) * (3.0 / 8.0)
        else:
            bpi = float('nan')
        set23_threshold = THRESHOLDS['SET23_CLINICAL']

    if not math.isnan(bpi):
        bpi = max(0.0, bpi)

    result.stage_risk = stage_risk
    result.node_risk = node_risk
    result.bpi = bpi

    # SET2,3 calculation
    if not math.isnan(bpi):
        set23 = (0.5056 * bpi) + (0.7542 * calibrated_set_erpr)
        result.set23 = set23
        result.set23_category = "High" if set23 >= set23_threshold else "Low"
    else:
        result.set23 = float('nan')
        result.set23_category = "N/A"

    return result


def process_batch(
    net_mfi_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    accession_df: pd.DataFrame,
    nc_sample_ids: List[str]
) -> Tuple[List[EAIResult], List[str]]:
    """
    Process a batch of samples through the EAI algorithm.

    Args:
        net_mfi_df: DataFrame with Sample column and gene columns (Net MFI values)
        counts_df: DataFrame with Sample column and gene columns (bead counts)
        accession_df: DataFrame with sample info (sample_id, sample_type, staging)
        nc_sample_ids: List of negative control sample IDs

    Returns:
        Tuple of (List of EAIResult objects, List of skipped sample IDs)
    """
    # Calculate NC adjustment
    nc_adjustment = calculate_negative_control_adjustment(net_mfi_df, nc_sample_ids)

    results = []
    skipped_samples = []

    for _, row in net_mfi_df.iterrows():
        sample_id = str(row['Sample']).strip()

        # Skip negative controls in results (they're used for adjustment)
        if sample_id in nc_sample_ids:
            continue

        # Get gene MFI values
        gene_mfi = {}
        for gene in ALL_GENES:
            if gene in row.index:
                gene_mfi[gene] = float(row[gene]) if pd.notna(row[gene]) else 0.0
            else:
                gene_mfi[gene] = 0.0

        # Get bead counts
        bead_counts = None
        if counts_df is not None:
            count_row = counts_df[counts_df['Sample'] == sample_id]
            if len(count_row) > 0:
                count_row = count_row.iloc[0]
                bead_counts = {}
                for gene in ALL_GENES:
                    if gene in count_row.index:
                        val = count_row[gene]
                        bead_counts[gene] = int(val) if pd.notna(val) else 100
                    else:
                        bead_counts[gene] = 100

        # Get accession info — exact match required, skip if no match
        acc_match = accession_df[accession_df['sample_id'].str.strip() == sample_id]

        if len(acc_match) == 0:
            # Per Federico (Apr 8, 2026): skip samples not in accession file
            skipped_samples.append(sample_id)
            continue

        acc_row = acc_match.iloc[0]
        sample_type = str(acc_row.get('sample_type', 'Clinical'))
        tumor_stage = str(acc_row.get('tumor_stage', ''))
        tumor_node = str(acc_row.get('tumor_node', ''))
        tumor_size_mm = acc_row.get('tumor_size_mm', 0.0)
        positive_lymph_nodes = acc_row.get('positive_lymph_nodes', 0)

        # Handle NaN values
        if pd.isna(tumor_size_mm):
            tumor_size_mm = 0.0
        if pd.isna(positive_lymph_nodes):
            positive_lymph_nodes = 0

        # Calculate EAI
        result = calculate_eai(
            gene_mfi=gene_mfi,
            bead_counts=bead_counts,
            nc_adjustment=nc_adjustment,
            sample_id=sample_id,
            sample_type=sample_type,
            tumor_stage=tumor_stage,
            tumor_node=tumor_node,
            tumor_size_mm=tumor_size_mm,
            positive_lymph_nodes=int(positive_lymph_nodes)
        )

        results.append(result)

    return results, skipped_samples


def results_to_dataframe(results: List[EAIResult]) -> pd.DataFrame:
    """Convert list of EAIResult to pandas DataFrame for display/export."""
    rows = []
    for r in results:
        # Determine bead count display
        if r.qc_bead_count_passed:
            bead_display = 'PASS'
        else:
            bead_display = f'FAIL ({r.min_bead_count})'

        rows.append({
            'Sample': r.sample_id,
            'Bead Count': bead_display,
            'Min Bead Count': r.min_bead_count,
            'Mean Ref Genes': round(r.avg_reference_genes, 2) if not math.isnan(r.avg_reference_genes) else 'N/A',
            'Background': round(r.background, 2) if not math.isnan(r.background) else 0.0,
            'PASSED (Y/N)': 'PASS' if r.qc_passed else 'FAIL',
            'Accession Match': 'Yes' if r.qc_accession_matched else 'No',
            'SET ER/PR': round(r.calibrated_set_erpr, 2),
            'ESR1': round(r.calibrated_esr1, 2),
            'ERBB2': round(r.calibrated_erbb2, 2),
            'PGR': round(r.calibrated_pgr, 2),
            'AURKA': round(r.calibrated_aurka, 2),
            'Sample Type': r.sample_type,
            'cT': r.tumor_stage if not r.sample_type.lower() == 'pathologic' else '',
            'cN': r.tumor_node if not r.sample_type.lower() == 'pathologic' else '',
            'Tumor size (mm)': r.tumor_size_mm if r.sample_type.lower() == 'pathologic' else '',
            'pT': r.tumor_stage if r.sample_type.lower() == 'pathologic' else '',
            '# Nodes': r.positive_lymph_nodes if r.sample_type.lower() == 'pathologic' else '',
            'pN': r.tumor_node if r.sample_type.lower() == 'pathologic' else '',
            'SET2,3': round(r.set23, 2) if not math.isnan(r.set23) else 'N/A',
            'SET2,3 Category': r.set23_category,
            'ESR1 Risk': r.esr1_risk,
            'ERBB2 Risk': r.erbb2_risk,
            'PGR Risk': r.pgr_risk,
            'PGR Level': r.pgr_risk_level,
            'AURKA Risk': round(r.aurka_risk, 2),
            'RNA4': round(r.rna4_score, 2),
            'Stage Risk': round(r.stage_risk, 2) if not math.isnan(r.stage_risk) else 'N/A',
            'Node Risk': round(r.node_risk, 2) if not math.isnan(r.node_risk) else 'N/A',
            'BPI': round(r.bpi, 2) if not math.isnan(r.bpi) else 'N/A',
            'QC Fail Reason': r.qc_fail_reason if r.qc_fail_reason else '',
            'Algorithm Version': r.algo_version
        })

    return pd.DataFrame(rows)
