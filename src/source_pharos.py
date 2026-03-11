from __future__ import annotations

import logging
import requests

from evidence_schema import (
    ConfidenceLevel,
    EvidenceCategory,
    EvidenceItem,
    EvidenceSource,
    TargetDevelopmentLevel,
    TargetProfile,
)

logger = logging.getLogger(__name__)

PHAROS_URL = "https://pharos-api.ncats.io/graphql"

QUERY = """
query TargetDetails($gene: String!) {
  targets(filter: { associatedTarget: $gene }) {
    targets {
      name
      sym
      tdl
      novelty
      fam
      uniprot
      description
    }
  }
}
"""

TDL_MAP = {
    "Tclin": TargetDevelopmentLevel.TCLIN,
    "Tchem": TargetDevelopmentLevel.TCHEM,
    "Tbio": TargetDevelopmentLevel.TBIO,
    "Tdark": TargetDevelopmentLevel.TDARK,
}

TDL_SCORES = {"Tclin": 1.0, "Tchem": 0.75, "Tbio": 0.4, "Tdark": 0.1}

FAMILY_SCORES = {
    "Kinase": 0.9,
    "GPCR": 0.9,
    "Ion Channel": 0.85,
    "Nuclear Receptor": 0.8,
    "Enzyme": 0.7,
    "Transporter": 0.65,
    "Transcription Factor": 0.4,
}


def run_graphql(query, variables, timeout=30):
    try:
        resp = requests.post(
            PHAROS_URL,
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        resp.raise_for_status()
        data = resp.json()
        if "errors" in data:
            logger.warning("PHAROS GraphQL errors: %s", data["errors"])
        return data.get("data")
    except requests.RequestException as e:
        logger.error("PHAROS request failed: %s", e)
        return None


def fetch_pharos_evidence(gene_symbol):
    items = []
    data = run_graphql(QUERY, {"gene": gene_symbol})

    target_data = None
    if data and data.get("targets") and data["targets"].get("targets"):
        target_list = data["targets"]["targets"]
        for t in target_list:
            if t.get("sym", "").upper() == gene_symbol.upper():
                target_data = t
                break
        if not target_data and target_list:
            target_data = target_list[0]

    if not target_data:
        logger.warning("No PHAROS data for gene: %s", gene_symbol)
        return items

    tdl = target_data.get("tdl", "Unknown")
    items.append(EvidenceItem(
        source=EvidenceSource.PHAROS,
        category=EvidenceCategory.DRUGGABILITY,
        feature_name="target_development_level",
        value=tdl,
        score=TDL_SCORES.get(tdl, 0.0),
        confidence=ConfidenceLevel.HIGH,
        description=f"PHAROS Target Development Level: {tdl}",
        metadata={"tdl": tdl},
    ))

    novelty = target_data.get("novelty")
    if novelty is not None:
        novelty_score = min(novelty / 10.0, 1.0) if novelty else 0.5
        items.append(EvidenceItem(
            source=EvidenceSource.PHAROS,
            category=EvidenceCategory.DRUGGABILITY,
            feature_name="pharos_novelty",
            value=novelty,
            score=novelty_score,
            confidence=ConfidenceLevel.MEDIUM,
            description=f"PHAROS novelty score: {novelty:.3f}",
            metadata={"novelty": novelty},
        ))

    fam = target_data.get("fam")
    if fam:
        items.append(EvidenceItem(
            source=EvidenceSource.PHAROS,
            category=EvidenceCategory.PROTEIN_STRUCTURE,
            feature_name="target_family",
            value=fam,
            score=FAMILY_SCORES.get(fam, 0.5),
            confidence=ConfidenceLevel.HIGH,
            description=f"Target family: {fam}",
            metadata={"family": fam},
        ))

    desc = target_data.get("description")
    if desc:
        items.append(EvidenceItem(
            source=EvidenceSource.PHAROS,
            category=EvidenceCategory.PHENOTYPE,
            feature_name="function_description",
            value=desc[:300],
            score=0.5,
            confidence=ConfidenceLevel.HIGH,
            description="PHAROS function description",
            metadata={"full_description": desc},
        ))

    logger.info("Fetched %d PHAROS evidence items for %s", len(items), gene_symbol)
    return items


def populate_profile_from_pharos(profile):
    items = fetch_pharos_evidence(profile.gene_symbol)
    for item in items:
        profile.add_evidence(item)

    for item in items:
        if item.feature_name == "target_development_level":
            profile.tdl = TDL_MAP.get(str(item.value), TargetDevelopmentLevel.UNKNOWN)
            break

    if not profile.uniprot_id:
        data = run_graphql(QUERY, {"gene": profile.gene_symbol})
        if data and data.get("targets") and data["targets"].get("targets"):
            for t in data["targets"]["targets"]:
                if t.get("sym", "").upper() == profile.gene_symbol.upper():
                    uid = t.get("uniprot")
                    if uid:
                        profile.uniprot_id = uid
                    break

    return profile
