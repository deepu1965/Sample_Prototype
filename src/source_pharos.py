from __future__ import annotations

import logging
import requests

from evidence_schema import (
    ConfidenceLevel,
    EvidenceCategory,
    EvidenceCollection,
    EvidenceContext,
    EvidenceItem,
    EvidenceMapping,
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
        target_id=gene_symbol.upper(),
        source=EvidenceSource.PHAROS,
        data_type="tractability",
        category=EvidenceCategory.DRUGGABILITY,
        feature_name="target_development_level",
        value=tdl,
        score=TDL_SCORES.get(tdl, 0.0),
        disease_context="pan_disease",
        confidence=ConfidenceLevel.HIGH,
        confidence_score=0.92,
        description=f"PHAROS Target Development Level: {tdl}",
        context_id=f"pharos_target_{gene_symbol.upper()}",
        collection_id="pharos_tractability",
        metadata={"tdl": tdl},
        provenance={
            "quality_tier": "high",
            "replicated": True,
            "sample_size": 1,
            "freshness_days": 30,
        },
    ))

    novelty = target_data.get("novelty")
    if novelty is not None:
        novelty_score = min(novelty / 10.0, 1.0) if novelty else 0.5
        items.append(EvidenceItem(
            target_id=gene_symbol.upper(),
            source=EvidenceSource.PHAROS,
            data_type="novelty",
            category=EvidenceCategory.DRUGGABILITY,
            feature_name="pharos_novelty",
            value=novelty,
            score=novelty_score,
            disease_context="pan_disease",
            confidence=ConfidenceLevel.MEDIUM,
            confidence_score=0.78,
            description=f"PHAROS novelty score: {novelty:.3f}",
            context_id=f"pharos_target_{gene_symbol.upper()}",
            collection_id="pharos_tractability",
            metadata={"novelty": novelty},
            provenance={
                "quality_tier": "medium",
                "replicated": True,
                "sample_size": 1,
                "freshness_days": 45,
            },
        ))

    fam = target_data.get("fam")
    if fam:
        items.append(EvidenceItem(
            target_id=gene_symbol.upper(),
            source=EvidenceSource.PHAROS,
            data_type="protein_family",
            category=EvidenceCategory.PROTEIN_STRUCTURE,
            feature_name="target_family",
            value=fam,
            score=FAMILY_SCORES.get(fam, 0.5),
            disease_context="pan_disease",
            confidence=ConfidenceLevel.HIGH,
            confidence_score=0.88,
            description=f"Target family: {fam}",
            context_id=f"pharos_target_{gene_symbol.upper()}",
            collection_id="pharos_structure",
            metadata={"family": fam},
            provenance={
                "quality_tier": "high",
                "replicated": True,
                "sample_size": 1,
                "freshness_days": 30,
            },
        ))

    desc = target_data.get("description")
    if desc:
        items.append(EvidenceItem(
            target_id=gene_symbol.upper(),
            source=EvidenceSource.PHAROS,
            data_type="functional_annotation",
            category=EvidenceCategory.PHENOTYPE,
            feature_name="function_description",
            value=desc[:300],
            score=0.5,
            disease_context="pan_disease",
            confidence=ConfidenceLevel.HIGH,
            confidence_score=0.85,
            description="PHAROS function description",
            context_id=f"pharos_target_{gene_symbol.upper()}",
            collection_id="pharos_annotation",
            metadata={"full_description": desc},
            provenance={
                "quality_tier": "high",
                "replicated": True,
                "sample_size": 1,
                "freshness_days": 30,
            },
        ))

    logger.info("Fetched %d PHAROS evidence items for %s", len(items), gene_symbol)
    return items


def populate_profile_from_pharos(profile):
    profile.add_context(
        EvidenceContext(
            context_id=f"pharos_target_{profile.gene_symbol}",
            source=EvidenceSource.PHAROS,
            context_type="target_profile",
            disease_context="pan_disease",
            metadata={"provider": "pharos"},
        )
    )
    profile.add_collection(
        EvidenceCollection(
            collection_id="pharos_tractability",
            collection_type="knowledge_base",
            source=EvidenceSource.PHAROS,
            metadata={"provider": "pharos"},
        )
    )
    profile.add_collection(
        EvidenceCollection(
            collection_id="pharos_structure",
            collection_type="knowledge_base",
            source=EvidenceSource.PHAROS,
            metadata={"provider": "pharos"},
        )
    )
    profile.add_collection(
        EvidenceCollection(
            collection_id="pharos_annotation",
            collection_type="knowledge_base",
            source=EvidenceSource.PHAROS,
            metadata={"provider": "pharos"},
        )
    )
    profile.add_mapping(
        EvidenceMapping(
            source_key=profile.gene_symbol,
            canonical_key=profile.gene_symbol,
            mapping_type="symbol_identity",
            metadata={"source": "pharos"},
        )
    )

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
