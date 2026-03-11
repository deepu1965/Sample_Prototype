from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional


class EvidenceSource(Enum):
    PHAROS = "pharos"
    DEPMAP = "depmap"
    OPENTARGETS = "opentargets"
    UNIPROT = "uniprot"
    STRING = "string"
    GENE_ONTOLOGY = "gene_ontology"
    MANUAL = "manual"


class EvidenceCategory(Enum):
    DRUGGABILITY = "druggability"
    GENETIC_CONSTRAINT = "genetic_constraint"
    PROTEIN_INTERACTION = "protein_interaction"
    EXPRESSION = "expression"
    DISEASE_ASSOCIATION = "disease_association"
    ESSENTIALITY = "essentiality"
    PATHWAY = "pathway"
    PROTEIN_STRUCTURE = "protein_structure"
    PHENOTYPE = "phenotype"


class TargetDevelopmentLevel(Enum):
    TCLIN = "Tclin"
    TCHEM = "Tchem"
    TBIO = "Tbio"
    TDARK = "Tdark"
    UNKNOWN = "Unknown"


class ConfidenceLevel(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class EvidenceItem:
    source: EvidenceSource
    category: EvidenceCategory
    feature_name: str
    value: object
    score: float = 0.0
    confidence: ConfidenceLevel = ConfidenceLevel.MEDIUM
    description: str = ""
    metadata: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())

    @property
    def uid(self) -> str:
        raw = f"{self.source.value}:{self.category.value}:{self.feature_name}"
        return hashlib.md5(raw.encode()).hexdigest()[:12]

    def to_dict(self) -> dict:
        return {
            "source": self.source.value,
            "category": self.category.value,
            "feature_name": self.feature_name,
            "value": self.value,
            "score": self.score,
            "confidence": self.confidence.value,
            "description": self.description,
            "metadata": self.metadata,
            "timestamp": self.timestamp,
            "uid": self.uid,
        }


@dataclass
class TargetProfile:
    gene_symbol: str
    uniprot_id: str = ""
    tdl: TargetDevelopmentLevel = TargetDevelopmentLevel.UNKNOWN
    evidence: list[EvidenceItem] = field(default_factory=list)
    aggregate_score: Optional[float] = None

    def add_evidence(self, item):
        existing_uids = {e.uid for e in self.evidence}
        if item.uid not in existing_uids:
            self.evidence.append(item)

    def get_evidence_by_source(self, source):
        return [e for e in self.evidence if e.source == source]

    def get_evidence_by_category(self, category):
        return [e for e in self.evidence if e.category == category]

    def compute_aggregate_score(self, category_weights=None):
        if not self.evidence:
            self.aggregate_score = 0.0
            return 0.0

        if category_weights is None:
            category_weights = {cat: 1.0 for cat in EvidenceCategory}

        cat_scores = {}
        for e in self.evidence:
            cat_scores.setdefault(e.category, []).append(e.score)

        weighted_sum = 0.0
        weight_sum = 0.0
        for cat, scores in cat_scores.items():
            w = category_weights.get(cat, 1.0)
            mean_score = sum(scores) / len(scores)
            weighted_sum += w * mean_score
            weight_sum += w

        self.aggregate_score = weighted_sum / weight_sum if weight_sum else 0.0
        return self.aggregate_score

    def summary(self):
        return {
            "gene_symbol": self.gene_symbol,
            "uniprot_id": self.uniprot_id,
            "tdl": self.tdl.value,
            "n_evidence_items": len(self.evidence),
            "sources": list({e.source.value for e in self.evidence}),
            "categories": list({e.category.value for e in self.evidence}),
            "aggregate_score": self.aggregate_score,
        }

    def to_json(self, indent=2):
        payload = self.summary()
        payload["evidence"] = [e.to_dict() for e in self.evidence]
        return json.dumps(payload, indent=indent, default=str)
