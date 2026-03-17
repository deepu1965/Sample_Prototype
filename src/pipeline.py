from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

from conflict_resolution import ConflictResolver
from evidence_schema import (
    EvidenceCategory,
    EvidenceSource,
    TargetProfile,
)
from source_pharos import populate_profile_from_pharos
from source_depmap import populate_profile_from_depmap
from source_opentargets import populate_profile_from_opentargets

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
    use_opentargets: bool = True
    depmap_builtin: bool = True
    opentargets_builtin: bool = False
    human_override_path: str = "human_review_overrides.json"
    profiles: list = field(default_factory=list)
    resolver: ConflictResolver = field(default_factory=ConflictResolver)

    def run(self, gene_symbols):
        self.profiles = []
        human_overrides = self._load_human_overrides()

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

            if self.use_opentargets:
                try:
                    populate_profile_from_opentargets(profile, use_builtin=self.opentargets_builtin)
                    logger.info(
                        "  Open Targets: %d items",
                        len(profile.get_evidence_by_source(EvidenceSource.OPENTARGETS)),
                    )
                except Exception as e:
                    logger.warning("  Open Targets failed for %s: %s", sym, e)

            profile.compute_aggregate_score(self.category_weights)
            override = human_overrides.get(profile.gene_symbol, {})
            final_score = self.resolver.resolve_profile(
                profile,
                category_weights=self.category_weights,
                human_override=override,
            )
            logger.info(
                "  Final score: %.4f  (%d evidence items, %d conflict flags)",
                final_score,
                len(profile.evidence),
                len(profile.conflict_flags),
            )

            self.profiles.append(profile)

        self.profiles.sort(key=lambda p: p.aggregate_score or 0.0, reverse=True)
        return self.profiles

    def ranking_table(self):
        lines = [
            f"{'Rank':<5} {'Gene':<10} {'TDL':<8} {'Score':>8} {'#Evidence':>10} {'#Flags':>8} {'Sources'}",
            "─" * 92,
        ]
        for i, p in enumerate(self.profiles, 1):
            sources = ", ".join(sorted({e.source.value for e in p.evidence}))
            lines.append(
                f"{i:<5} {p.gene_symbol:<10} {p.tdl.value:<8} {p.aggregate_score or 0:>8.4f} "
                f"{len(p.evidence):>10} {len(p.conflict_flags):>8} {sources}"
            )
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

        if profile.layer_scores:
            lines.append(
                f"  Layer scores: L1={profile.layer_scores.get('layer1_concordance', 0):.4f} | "
                f"Penalty={profile.layer_scores.get('layer2_penalty', 0):.4f} | "
                f"Final={profile.layer_scores.get('layer3_final', 0):.4f}"
            )

        if profile.conflict_flags:
            lines.append("  Conflict Flags:")
            for flag in profile.conflict_flags:
                lines.append(f"    - [{flag['severity']}] {flag['id']}: {flag['message']}")

        if profile.human_notes:
            lines.append("  Human Notes:")
            for note in profile.human_notes:
                lines.append(f"    - {note}")

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

    def _load_human_overrides(self):
        path = Path(self.human_override_path)
        if not path.exists():
            return {}
        try:
            data = json.loads(path.read_text())
            return {str(k).upper(): v for k, v in data.items()}
        except Exception as exc:
            logger.warning("Failed to load human review overrides: %s", exc)
            return {}
