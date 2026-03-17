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
class EvidenceContext:
    context_id: str
    source: EvidenceSource
    context_type: str
    disease_context: str = "pan_disease"
    metadata: dict = field(default_factory=dict)


@dataclass
class EvidenceCollection:
    collection_id: str
    collection_type: str
    source: EvidenceSource
    metadata: dict = field(default_factory=dict)


@dataclass
class EvidenceMapping:
    source_key: str
    canonical_key: str
    mapping_type: str
    metadata: dict = field(default_factory=dict)


@dataclass
class EvidenceItem:
    target_id: str
    source: EvidenceSource
    data_type: str
    category: EvidenceCategory
    feature_name: str
    value: object
    score: float = 0.0
    disease_context: str = "pan_disease"
    confidence: ConfidenceLevel = ConfidenceLevel.MEDIUM
    confidence_score: float = 0.5
    description: str = ""
    context_id: str = ""
    collection_id: str = ""
    metadata: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())

    @property
    def uid(self) -> str:
        raw = f"{self.source.value}:{self.category.value}:{self.feature_name}"
        return hashlib.md5(raw.encode()).hexdigest()[:12]

    def to_dict(self) -> dict:
        return {
            "target_id": self.target_id,
            "source": self.source.value,
            "data_type": self.data_type,
            "category": self.category.value,
            "feature_name": self.feature_name,
            "value": self.value,
            "score": self.score,
            "disease_context": self.disease_context,
            "confidence": self.confidence.value,
            "confidence_score": self.confidence_score,
            "description": self.description,
            "context_id": self.context_id,
            "collection_id": self.collection_id,
            "metadata": self.metadata,
            "provenance": self.provenance,
            "timestamp": self.timestamp,
            "uid": self.uid,
        }


@dataclass
class TargetProfile:
    gene_symbol: str
    uniprot_id: str = ""
    tdl: TargetDevelopmentLevel = TargetDevelopmentLevel.UNKNOWN
    evidence: list[EvidenceItem] = field(default_factory=list)
    contexts: list[EvidenceContext] = field(default_factory=list)
    collections: list[EvidenceCollection] = field(default_factory=list)
    mappings: list[EvidenceMapping] = field(default_factory=list)
    conflict_flags: list[dict] = field(default_factory=list)
    human_notes: list[str] = field(default_factory=list)
    layer_scores: dict = field(default_factory=dict)
    aggregate_score: Optional[float] = None

    def add_evidence(self, item):
        existing_uids = {e.uid for e in self.evidence}
        if item.uid not in existing_uids:
            self.evidence.append(item)

    def add_context(self, context):
        if not any(c.context_id == context.context_id for c in self.contexts):
            self.contexts.append(context)

    def add_collection(self, collection):
        if not any(c.collection_id == collection.collection_id for c in self.collections):
            self.collections.append(collection)

    def add_mapping(self, mapping):
        if not any(
            m.source_key == mapping.source_key and m.canonical_key == mapping.canonical_key
            for m in self.mappings
        ):
            self.mappings.append(mapping)

    def add_human_note(self, note):
        if note:
            self.human_notes.append(note)

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
            "n_contexts": len(self.contexts),
            "n_collections": len(self.collections),
            "sources": list({e.source.value for e in self.evidence}),
            "categories": list({e.category.value for e in self.evidence}),
            "n_conflicts": len(self.conflict_flags),
            "layer_scores": self.layer_scores,
            "aggregate_score": self.aggregate_score,
        }

    def to_json(self, indent=2):
        payload = self.summary()
        payload["evidence"] = [e.to_dict() for e in self.evidence]
        payload["contexts"] = [c.__dict__ | {"source": c.source.value} for c in self.contexts]
        payload["collections"] = [
            c.__dict__ | {"source": c.source.value} for c in self.collections
        ]
        payload["mappings"] = [m.__dict__ for m in self.mappings]
        payload["conflict_flags"] = self.conflict_flags
        payload["human_notes"] = self.human_notes
        return json.dumps(payload, indent=indent, default=str)
