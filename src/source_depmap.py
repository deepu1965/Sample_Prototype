from __future__ import annotations

import logging
import requests

from evidence_schema import (
    ConfidenceLevel,
    EvidenceCategory,
    EvidenceItem,
    EvidenceSource,
    TargetProfile,
)

logger = logging.getLogger(__name__)

BUILTIN_DEPMAP = {
    "EGFR": {"mean_dep": -0.35, "n_dependent_lines": 42, "total_lines": 1095, "common_essential": False},
    "TP53": {"mean_dep": 0.05, "n_dependent_lines": 3, "total_lines": 1095, "common_essential": False},
    "BRCA1": {"mean_dep": -0.18, "n_dependent_lines": 15, "total_lines": 1095, "common_essential": False},
    "KRAS": {"mean_dep": -0.52, "n_dependent_lines": 98, "total_lines": 1095, "common_essential": False},
    "MYC": {"mean_dep": -0.78, "n_dependent_lines": 412, "total_lines": 1095, "common_essential": True},
    "BRAF": {"mean_dep": -0.22, "n_dependent_lines": 58, "total_lines": 1095, "common_essential": False},
    "PIK3CA": {"mean_dep": -0.25, "n_dependent_lines": 35, "total_lines": 1095, "common_essential": False},
    "CDK4": {"mean_dep": -0.30, "n_dependent_lines": 45, "total_lines": 1095, "common_essential": False},
    "AKT1": {"mean_dep": -0.15, "n_dependent_lines": 20, "total_lines": 1095, "common_essential": False},
    "MTOR": {"mean_dep": -0.65, "n_dependent_lines": 320, "total_lines": 1095, "common_essential": True},
    "RPS6": {"mean_dep": -1.10, "n_dependent_lines": 950, "total_lines": 1095, "common_essential": True},
    "PCNA": {"mean_dep": -1.05, "n_dependent_lines": 980, "total_lines": 1095, "common_essential": True},
}


def dep_score_to_essentiality(mean_dep):
    if mean_dep >= 0:
        return 0.0
    return min(abs(mean_dep), 1.0)


def selectivity_score(n_dependent, total):
    if total == 0:
        return 0.0
    return 1.0 - (n_dependent / total)


def fetch_depmap_evidence(gene_symbol, use_builtin=True):
    items = []
    gene_data = None

    if use_builtin:
        gene_data = BUILTIN_DEPMAP.get(gene_symbol.upper())
    else:
        gene_data = try_fetch_depmap_api(gene_symbol)
        if gene_data is None:
            gene_data = BUILTIN_DEPMAP.get(gene_symbol.upper())

    if gene_data is None:
        logger.info("No DepMap data available for %s", gene_symbol)
        return items

    mean_dep = gene_data["mean_dep"]
    n_dep = gene_data["n_dependent_lines"]
    total = gene_data["total_lines"]
    common_essential = gene_data.get("common_essential", False)

    items.append(EvidenceItem(
        source=EvidenceSource.DEPMAP,
        category=EvidenceCategory.ESSENTIALITY,
        feature_name="crispr_mean_dependency",
        value=mean_dep,
        score=dep_score_to_essentiality(mean_dep),
        confidence=ConfidenceLevel.HIGH,
        description=f"Mean CRISPR dependency score across {total} cell lines: {mean_dep:.3f} (more negative = more essential)",
        metadata={"mean_dep": mean_dep, "total_lines": total},
    ))

    dep_frac = n_dep / total if total > 0 else 0
    items.append(EvidenceItem(
        source=EvidenceSource.DEPMAP,
        category=EvidenceCategory.ESSENTIALITY,
        feature_name="dependent_cell_line_fraction",
        value=dep_frac,
        score=dep_frac,
        confidence=ConfidenceLevel.HIGH,
        description=f"{n_dep}/{total} cell lines ({dep_frac:.1%}) depend on {gene_symbol}",
        metadata={"n_dependent": n_dep, "total_lines": total},
    ))

    sel = selectivity_score(n_dep, total)
    items.append(EvidenceItem(
        source=EvidenceSource.DEPMAP,
        category=EvidenceCategory.ESSENTIALITY,
        feature_name="selectivity_score",
        value=sel,
        score=sel,
        confidence=ConfidenceLevel.MEDIUM,
        description=f"Selectivity score (1 = highly selective dependency): {sel:.3f}",
        metadata={"selectivity": sel},
    ))

    items.append(EvidenceItem(
        source=EvidenceSource.DEPMAP,
        category=EvidenceCategory.ESSENTIALITY,
        feature_name="common_essential",
        value=common_essential,
        score=0.0 if common_essential else 0.7,
        confidence=ConfidenceLevel.HIGH,
        description=f"{'IS' if common_essential else 'NOT'} a common essential gene. Common essentials are poor drug targets (toxic on inhibition).",
        metadata={"common_essential": common_essential},
    ))

    logger.info("Fetched %d DepMap evidence items for %s", len(items), gene_symbol)
    return items


def try_fetch_depmap_api(gene_symbol):
    try:
        url = f"https://depmap.org/portal/api/gene/{gene_symbol}/dep_summary"
        resp = requests.get(url, timeout=15)
        if resp.status_code != 200:
            return None
        data = resp.json()
        return {
            "mean_dep": data.get("mean_chronos", 0.0),
            "n_dependent_lines": data.get("n_dependent", 0),
            "total_lines": data.get("n_total", 1095),
            "common_essential": data.get("is_common_essential", False),
        }
    except Exception as e:
        logger.debug("DepMap API request failed: %s", e)
        return None


def populate_profile_from_depmap(profile, use_builtin=True):
    items = fetch_depmap_evidence(profile.gene_symbol, use_builtin=use_builtin)
    for item in items:
        profile.add_evidence(item)
    return profile
