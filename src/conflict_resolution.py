from __future__ import annotations

import math
from dataclasses import dataclass

from evidence_schema import EvidenceCategory, EvidenceSource, TargetProfile


SOURCE_BASE_WEIGHTS = {
    EvidenceSource.DEPMAP: 0.9,
    EvidenceSource.OPENTARGETS: 0.85,
    EvidenceSource.PHAROS: 0.8,
}


@dataclass
class ConflictResolverConfig:
    disagreement_threshold: float = 0.35
    severe_disagreement_threshold: float = 0.55


class ConflictResolver:
    def __init__(self, config=None):
        self.config = config or ConflictResolverConfig()

    def resolve_profile(self, profile, category_weights, human_override=None):
        layer1 = self._layer1_source_weighted_concordance(profile)
        flags = self._layer2_biological_plausibility(profile, layer1)
        final_score = self._layer3_human_in_loop(profile, layer1["score"], flags, category_weights, human_override)
        profile.layer_scores = {
            "layer1_concordance": round(layer1["score"], 4),
            "layer2_penalty": round(layer1["penalty"], 4),
            "layer3_final": round(final_score, 4),
        }
        profile.conflict_flags = flags
        profile.aggregate_score = final_score
        return final_score

    def _source_weight(self, evidence):
        base = SOURCE_BASE_WEIGHTS.get(evidence.source, 0.7)
        quality_factor = float(evidence.confidence_score or 0.5)
        provenance = evidence.provenance or {}
        sample_size = float(provenance.get("sample_size") or 1.0)
        replicated = 1.0 if provenance.get("replicated", False) else 0.8
        freshness_days = float(provenance.get("freshness_days") or 365.0)
        freshness = max(0.5, min(1.0, 1.0 - (freshness_days / 730.0)))
        sample_factor = max(0.7, min(1.2, math.log10(sample_size + 10) / 1.5))
        return base * quality_factor * replicated * freshness * sample_factor

    def _layer1_source_weighted_concordance(self, profile):
        category_scores = {}
        disagreements = []
        weighted_total = 0.0
        total_weight = 0.0

        for evidence in profile.evidence:
            category_scores.setdefault(evidence.category, []).append(evidence)

        for category, evidences in category_scores.items():
            vals = [float(e.score) for e in evidences]
            mean_val = sum(vals) / len(vals)

            source_buckets = {}
            for e in evidences:
                source_buckets.setdefault(e.source.value, []).append(float(e.score))
            source_means = [sum(v) / len(v) for v in source_buckets.values()]

            if len(source_means) <= 1:
                disagreement = 0.0
            else:
                stdev = self._stdev(source_means)
                disagreement = min(stdev * 2.0, 1.0)

            disagreements.append(disagreement)

            score_num = 0.0
            score_den = 0.0
            for e in evidences:
                w = self._source_weight(e)
                score_num += w * float(e.score)
                score_den += w
            weighted_mean = score_num / score_den if score_den else mean_val

            category_scores[category] = {
                "weighted_mean": weighted_mean,
                "disagreement": disagreement,
                "n": len(evidences),
            }

        avg_disagreement = sum(disagreements) / len(disagreements) if disagreements else 0.0
        concordance_penalty = min(avg_disagreement * 0.35, 0.25)

        for category, data in category_scores.items():
            w = 1.0
            weighted_total += w * data["weighted_mean"]
            total_weight += w

        base_score = weighted_total / total_weight if total_weight else 0.0
        layer1_score = max(0.0, base_score - concordance_penalty)
        return {
            "score": layer1_score,
            "penalty": concordance_penalty,
            "categories": category_scores,
            "avg_disagreement": avg_disagreement,
        }

    def _layer2_biological_plausibility(self, profile, layer1):
        flags = []

        dep_common = self._find_value(profile, EvidenceSource.DEPMAP, "common_essential")
        dep_selectivity = self._find_score(profile, EvidenceSource.DEPMAP, "selectivity_score")
        pharos_tdl_score = self._find_score(profile, EvidenceSource.PHAROS, "target_development_level")
        ot_disease_score = self._find_score(profile, EvidenceSource.OPENTARGETS, "ot_max_disease_score")

        if dep_common is True and pharos_tdl_score is not None and pharos_tdl_score >= 0.75:
            flags.append({
                "id": "toxicity_vs_tractability",
                "severity": "major",
                "message": "Target is clinically tractable but appears common-essential in DepMap.",
            })

        if dep_selectivity is not None and dep_selectivity < 0.25 and ot_disease_score is not None and ot_disease_score > 0.75:
            flags.append({
                "id": "broad_essentiality_vs_disease_signal",
                "severity": "major",
                "message": "Strong disease association but low selectivity suggests toxicity risk.",
            })

        for category, item in layer1["categories"].items():
            disagreement = item["disagreement"]
            if disagreement >= self.config.severe_disagreement_threshold:
                flags.append({
                    "id": f"high_disagreement_{category.value}",
                    "severity": "major",
                    "message": f"High disagreement across sources in category {category.value}.",
                })
            elif disagreement >= self.config.disagreement_threshold:
                flags.append({
                    "id": f"medium_disagreement_{category.value}",
                    "severity": "minor",
                    "message": f"Moderate disagreement across sources in category {category.value}.",
                })

        return flags

    def _layer3_human_in_loop(self, profile, layer1_score, flags, category_weights, human_override):
        adjusted = layer1_score
        if flags:
            major = sum(1 for f in flags if f["severity"] == "major")
            minor = sum(1 for f in flags if f["severity"] == "minor")
            adjusted -= min(0.2, major * 0.07 + minor * 0.02)

        if human_override:
            note = human_override.get("note")
            if note:
                profile.add_human_note(note)

            resolved_ids = set(human_override.get("resolved_flags", []))
            if resolved_ids:
                flags[:] = [f for f in flags if f["id"] not in resolved_ids]

            multiplier = float(human_override.get("score_multiplier", 1.0))
            adjusted *= multiplier

            bonus = float(human_override.get("score_bonus", 0.0))
            adjusted += bonus

            category_multipliers = human_override.get("category_multipliers", {})
            if category_multipliers:
                temp_weights = dict(category_weights)
                for category_name, m in category_multipliers.items():
                    key = self._safe_category(category_name)
                    if key:
                        temp_weights[key] = temp_weights.get(key, 1.0) * float(m)
                baseline = profile.compute_aggregate_score(temp_weights)
                adjusted = (adjusted * 0.6) + (baseline * 0.4)

        return min(max(adjusted, 0.0), 1.0)

    def _find_score(self, profile, source, feature_name):
        for evidence in profile.evidence:
            if evidence.source == source and evidence.feature_name == feature_name:
                return float(evidence.score)
        return None

    def _find_value(self, profile, source, feature_name):
        for evidence in profile.evidence:
            if evidence.source == source and evidence.feature_name == feature_name:
                return evidence.value
        return None

    def _safe_category(self, category_name):
        for category in EvidenceCategory:
            if category.value == category_name or category.name == str(category_name).upper():
                return category
        return None

    def _stdev(self, values):
        if len(values) < 2:
            return 0.0
        mean_v = sum(values) / len(values)
        var = sum((v - mean_v) ** 2 for v in values) / len(values)
        return math.sqrt(var)
