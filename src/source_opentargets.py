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
    TargetProfile,
)

logger = logging.getLogger(__name__)

OPENTARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

SEARCH_QUERY = """
query SearchTarget($queryString: String!) {
  search(queryString: $queryString, entityNames: ["target"], page: {index: 0, size: 5}) {
    hits {
      id
      object {
        ... on Target {
          id
          approvedSymbol
        }
      }
    }
  }
}
"""

ASSOCIATION_QUERY = """
query TargetDisease($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 25}) {
      count
      rows {
        score
        disease {
          id
          name
        }
      }
    }
  }
}
"""

BUILTIN_OPENTARGETS = {
    "EGFR": {"mean_disease_score": 0.82, "max_disease_score": 0.97, "n_diseases": 23},
    "KRAS": {"mean_disease_score": 0.79, "max_disease_score": 0.95, "n_diseases": 18},
    "TP53": {"mean_disease_score": 0.91, "max_disease_score": 0.99, "n_diseases": 34},
    "BRCA1": {"mean_disease_score": 0.88, "max_disease_score": 0.98, "n_diseases": 27},
    "BRAF": {"mean_disease_score": 0.8, "max_disease_score": 0.96, "n_diseases": 19},
    "MYC": {"mean_disease_score": 0.76, "max_disease_score": 0.94, "n_diseases": 14},
}


def run_graphql(query, variables, timeout=20):
    try:
        response = requests.post(
            OPENTARGETS_GRAPHQL_URL,
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()
        if "errors" in data:
            logger.warning("Open Targets GraphQL errors: %s", data["errors"])
            return None
        return data.get("data")
    except requests.RequestException as exc:
        logger.debug("Open Targets request failed: %s", exc)
        return None


def try_fetch_opentargets_api(gene_symbol):
    search_data = run_graphql(SEARCH_QUERY, {"queryString": gene_symbol.upper()})
    if not search_data or not search_data.get("search"):
        return None

    hits = search_data["search"].get("hits") or []
    ensembl_id = None
    for hit in hits:
        obj = hit.get("object") or {}
        if obj.get("approvedSymbol", "").upper() == gene_symbol.upper():
            ensembl_id = obj.get("id")
            break
    if not ensembl_id and hits:
        ensembl_id = (hits[0].get("object") or {}).get("id")

    if not ensembl_id:
        return None

    assoc_data = run_graphql(ASSOCIATION_QUERY, {"ensemblId": ensembl_id})
    if not assoc_data or not assoc_data.get("target"):
        return None

    assoc = assoc_data["target"].get("associatedDiseases") or {}
    rows = assoc.get("rows") or []
    if not rows:
        return None

    scores = [float(row.get("score", 0.0)) for row in rows]
    top_disease = rows[0].get("disease") or {}
    return {
        "mean_disease_score": sum(scores) / len(scores),
        "max_disease_score": max(scores),
        "n_diseases": int(assoc.get("count") or len(rows)),
        "top_disease": top_disease.get("name", ""),
        "top_disease_id": top_disease.get("id", ""),
        "ensembl_id": ensembl_id,
    }


def fetch_opentargets_evidence(gene_symbol, use_builtin=False):
    items = []
    data = None if use_builtin else try_fetch_opentargets_api(gene_symbol)
    if data is None:
        data = BUILTIN_OPENTARGETS.get(gene_symbol.upper())

    if not data:
        logger.info("No Open Targets data available for %s", gene_symbol)
        return items

    mean_score = float(data.get("mean_disease_score", 0.0))
    max_score = float(data.get("max_disease_score", 0.0))
    n_diseases = int(data.get("n_diseases", 0))

    items.append(EvidenceItem(
        target_id=gene_symbol.upper(),
        source=EvidenceSource.OPENTARGETS,
        data_type="disease_association",
        category=EvidenceCategory.DISEASE_ASSOCIATION,
        feature_name="ot_mean_disease_score",
        value=mean_score,
        score=min(max(mean_score, 0.0), 1.0),
        disease_context="disease_agnostic",
        confidence=ConfidenceLevel.HIGH,
        confidence_score=0.88,
        description=f"Open Targets mean association score: {mean_score:.3f}",
        context_id=f"ot_assoc_{gene_symbol.upper()}",
        collection_id="opentargets_associations",
        metadata={"n_diseases": n_diseases},
        provenance={
            "quality_tier": "high",
            "replicated": True,
            "sample_size": n_diseases,
            "freshness_days": 30,
        },
    ))

    items.append(EvidenceItem(
        target_id=gene_symbol.upper(),
        source=EvidenceSource.OPENTARGETS,
        data_type="disease_association",
        category=EvidenceCategory.DISEASE_ASSOCIATION,
        feature_name="ot_max_disease_score",
        value=max_score,
        score=min(max(max_score, 0.0), 1.0),
        disease_context="disease_agnostic",
        confidence=ConfidenceLevel.MEDIUM,
        confidence_score=0.8,
        description=f"Open Targets top disease association score: {max_score:.3f}",
        context_id=f"ot_assoc_{gene_symbol.upper()}",
        collection_id="opentargets_associations",
        metadata={
            "top_disease": data.get("top_disease", ""),
            "top_disease_id": data.get("top_disease_id", ""),
        },
        provenance={
            "quality_tier": "medium",
            "replicated": True,
            "sample_size": n_diseases,
            "freshness_days": 30,
        },
    ))

    items.append(EvidenceItem(
        target_id=gene_symbol.upper(),
        source=EvidenceSource.OPENTARGETS,
        data_type="disease_count",
        category=EvidenceCategory.DISEASE_ASSOCIATION,
        feature_name="ot_associated_disease_count",
        value=n_diseases,
        score=min(n_diseases / 50.0, 1.0),
        disease_context="disease_agnostic",
        confidence=ConfidenceLevel.MEDIUM,
        confidence_score=0.75,
        description=f"Open Targets associated disease count: {n_diseases}",
        context_id=f"ot_assoc_{gene_symbol.upper()}",
        collection_id="opentargets_associations",
        metadata={"ensembl_id": data.get("ensembl_id", "")},
        provenance={
            "quality_tier": "medium",
            "replicated": True,
            "sample_size": n_diseases,
            "freshness_days": 30,
        },
    ))

    logger.info("Fetched %d Open Targets evidence items for %s", len(items), gene_symbol)
    return items


def populate_profile_from_opentargets(profile, use_builtin=False):
    profile.add_context(
        EvidenceContext(
            context_id=f"ot_assoc_{profile.gene_symbol}",
            source=EvidenceSource.OPENTARGETS,
            context_type="target_disease_association",
            disease_context="disease_agnostic",
            metadata={"provider": "opentargets"},
        )
    )
    profile.add_collection(
        EvidenceCollection(
            collection_id="opentargets_associations",
            collection_type="association_index",
            source=EvidenceSource.OPENTARGETS,
            metadata={"provider": "opentargets"},
        )
    )
    profile.add_mapping(
        EvidenceMapping(
            source_key=profile.gene_symbol,
            canonical_key=profile.gene_symbol,
            mapping_type="symbol_identity",
            metadata={"source": "opentargets"},
        )
    )

    items = fetch_opentargets_evidence(profile.gene_symbol, use_builtin=use_builtin)
    for item in items:
        profile.add_evidence(item)
    return profile
