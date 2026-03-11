from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field

from evidence_schema import (
    EvidenceCategory,
    EvidenceSource,
    TargetProfile,
)
from source_pharos import populate_profile_from_pharos
from source_depmap import populate_profile_from_depmap

logger = logging.getLogger(__name__)

DEFAULT_WEIGHTS = {
    EvidenceCategory.DRUGGABILITY: 2.0,
    EvidenceCategory.ESSENTIALITY: 1.5,
    EvidenceCategory.DISEASE_ASSOCIATION: 1.2,
    EvidenceCategory.PROTEIN_INTERACTION: 1.0,
    EvidenceCategory.PROTEIN_STRUCTURE: 1.0,
    EvidenceCategory.GENETIC_CONSTRAINT: 1.0,
    EvidenceCategory.EXPRESSION: 0.8,
    EvidenceCategory.PATHWAY: 0.8,
    EvidenceCategory.PHENOTYPE: 0.5,
}


@dataclass
class TargetPipeline:
    category_weights: dict = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))
    use_pharos: bool = True
    use_depmap: bool = True
    depmap_builtin: bool = True
    profiles: list = field(default_factory=list)

    def run(self, gene_symbols):
        self.profiles = []

        for sym in gene_symbols:
            logger.info("━━━ Processing %s ━━━", sym)
            profile = TargetProfile(gene_symbol=sym.upper())

            if self.use_pharos:
                try:
                    populate_profile_from_pharos(profile)
                    logger.info("  PHAROS: %d items", len(profile.get_evidence_by_source(EvidenceSource.PHAROS)))
                except Exception as e:
                    logger.warning("  PHAROS failed for %s: %s", sym, e)

            if self.use_depmap:
                try:
                    populate_profile_from_depmap(profile, use_builtin=self.depmap_builtin)
                    logger.info("  DepMap: %d items", len(profile.get_evidence_by_source(EvidenceSource.DEPMAP)))
                except Exception as e:
                    logger.warning("  DepMap failed for %s: %s", sym, e)

            profile.compute_aggregate_score(self.category_weights)
            logger.info("  Aggregate score: %.4f  (%d evidence items)", profile.aggregate_score or 0.0, len(profile.evidence))

            self.profiles.append(profile)

        self.profiles.sort(key=lambda p: p.aggregate_score or 0.0, reverse=True)
        return self.profiles

    def ranking_table(self):
        lines = [
            f"{'Rank':<5} {'Gene':<10} {'TDL':<8} {'Score':>8} {'#Evidence':>10} {'Sources'}",
            "─" * 75,
        ]
        for i, p in enumerate(self.profiles, 1):
            sources = ", ".join(sorted({e.source.value for e in p.evidence}))
            lines.append(f"{i:<5} {p.gene_symbol:<10} {p.tdl.value:<8} {p.aggregate_score or 0:>8.4f} {len(p.evidence):>10} {sources}")
        return "\n".join(lines)

    def to_json(self, indent=2):
        return json.dumps(
            [json.loads(p.to_json()) for p in self.profiles],
            indent=indent,
            default=str,
        )

    def evidence_breakdown(self, gene_symbol):
        profile = next((p for p in self.profiles if p.gene_symbol == gene_symbol.upper()), None)
        if profile is None:
            return f"No profile found for {gene_symbol}"

        lines = [
            f"\n{'=' * 60}",
            f"  Evidence Breakdown: {profile.gene_symbol}",
            f"  TDL: {profile.tdl.value}  |  Aggregate Score: {profile.aggregate_score or 0:.4f}",
            f"{'=' * 60}",
        ]

        by_cat = {}
        for e in profile.evidence:
            by_cat.setdefault(e.category.value, []).append(e)

        for cat, evidences in sorted(by_cat.items()):
            lines.append(f"\n  ▸ {cat.upper().replace('_', ' ')}")
            for e in evidences:
                lines.append(f"    [{e.source.value:<12}] {e.feature_name:<30} = {str(e.value):<20} score={e.score:.3f}  ({e.confidence.value})")
                if e.description:
                    lines.append(f"    {'':>12}  └─ {e.description}")

        lines.append(f"\n{'=' * 60}\n")
        return "\n".join(lines)
