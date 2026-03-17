#!/usr/bin/env python3

from __future__ import annotations

import argparse
import copy
import json
import random
import statistics
import time
import tracemalloc

from conflict_resolution import ConflictResolver
from evidence_schema import (
    ConfidenceLevel,
    EvidenceCategory,
    EvidenceItem,
    EvidenceSource,
    TargetProfile,
)
from pipeline import DEFAULT_WEIGHTS
from source_depmap import populate_profile_from_depmap
from source_opentargets import populate_profile_from_opentargets
from source_pharos import populate_profile_from_pharos


def make_item(
    gene,
    source,
    category,
    feature_name,
    score,
    value,
    confidence_score,
    disease_context,
    sample_size,
    freshness_days,
    replicated=True,
):
    return EvidenceItem(
        target_id=gene,
        source=source,
        data_type="synthetic",
        category=category,
        feature_name=feature_name,
        value=value,
        score=max(0.0, min(1.0, float(score))),
        disease_context=disease_context,
        confidence=ConfidenceLevel.HIGH if confidence_score >= 0.85 else ConfidenceLevel.MEDIUM,
        confidence_score=float(confidence_score),
        description=f"Synthetic evidence for {feature_name}",
        context_id=f"synthetic_ctx_{gene}",
        collection_id="synthetic_benchmark",
        metadata={"synthetic": True},
        provenance={
            "sample_size": int(sample_size),
            "replicated": bool(replicated),
            "freshness_days": float(freshness_days),
        },
    )


def build_synthetic_profile(gene, scenario, rng):
    profile = TargetProfile(gene_symbol=gene)

    if scenario == "concordant":
        dep_selectivity = rng.uniform(0.75, 0.95)
        dep_common = False
        tdl_score = rng.uniform(0.75, 1.0)
        ot_max = rng.uniform(0.7, 0.9)
        pharos_disease = max(0.0, ot_max - rng.uniform(0.0, 0.12))
    elif scenario == "moderate_conflict":
        dep_selectivity = rng.uniform(0.05, 0.22)
        dep_common = False
        tdl_score = rng.uniform(0.7, 0.95)
        ot_max = rng.uniform(0.82, 0.97)
        pharos_disease = rng.uniform(0.2, 0.45)
    else:
        dep_selectivity = rng.uniform(0.02, 0.15)
        dep_common = True
        tdl_score = rng.uniform(0.8, 1.0)
        ot_max = rng.uniform(0.88, 0.99)
        pharos_disease = rng.uniform(0.05, 0.25)

    dep_mean = rng.uniform(0.35, 0.85)

    profile.add_evidence(
        make_item(
            gene=gene,
            source=EvidenceSource.PHAROS,
            category=EvidenceCategory.DRUGGABILITY,
            feature_name="target_development_level",
            score=tdl_score,
            value="Tclin" if tdl_score >= 0.75 else "Tchem",
            confidence_score=0.9,
            disease_context="pan_disease",
            sample_size=1,
            freshness_days=20,
        )
    )

    profile.add_evidence(
        make_item(
            gene=gene,
            source=EvidenceSource.DEPMAP,
            category=EvidenceCategory.ESSENTIALITY,
            feature_name="selectivity_score",
            score=dep_selectivity,
            value=dep_selectivity,
            confidence_score=0.9,
            disease_context="cancer",
            sample_size=1200,
            freshness_days=120,
        )
    )

    profile.add_evidence(
        make_item(
            gene=gene,
            source=EvidenceSource.DEPMAP,
            category=EvidenceCategory.ESSENTIALITY,
            feature_name="crispr_mean_dependency",
            score=dep_mean,
            value=-dep_mean,
            confidence_score=0.88,
            disease_context="cancer",
            sample_size=1200,
            freshness_days=120,
        )
    )

    profile.add_evidence(
        make_item(
            gene=gene,
            source=EvidenceSource.DEPMAP,
            category=EvidenceCategory.ESSENTIALITY,
            feature_name="common_essential",
            score=0.0 if dep_common else 0.7,
            value=dep_common,
            confidence_score=0.9,
            disease_context="cancer",
            sample_size=1200,
            freshness_days=120,
        )
    )

    profile.add_evidence(
        make_item(
            gene=gene,
            source=EvidenceSource.OPENTARGETS,
            category=EvidenceCategory.DISEASE_ASSOCIATION,
            feature_name="ot_max_disease_score",
            score=ot_max,
            value=ot_max,
            confidence_score=0.86,
            disease_context="disease_agnostic",
            sample_size=300,
            freshness_days=35,
        )
    )

    profile.add_evidence(
        make_item(
            gene=gene,
            source=EvidenceSource.PHAROS,
            category=EvidenceCategory.DISEASE_ASSOCIATION,
            feature_name="pharos_disease_hint",
            score=pharos_disease,
            value=pharos_disease,
            confidence_score=0.78,
            disease_context="disease_agnostic",
            sample_size=30,
            freshness_days=60,
        )
    )

    return profile


def run_benchmark(n_targets, seed):
    rng = random.Random(seed)
    resolver = ConflictResolver()

    scenarios = ["concordant", "moderate_conflict", "severe_conflict"]
    expected_conflict = {"concordant": False, "moderate_conflict": True, "severe_conflict": True}

    records = []

    tracemalloc.start()
    t0 = time.perf_counter()

    for i in range(n_targets):
        scenario = scenarios[i % len(scenarios)]
        gene = f"SYN{i:04d}"

        profile = build_synthetic_profile(gene, scenario, rng)
        profile.compute_aggregate_score(DEFAULT_WEIGHTS)
        resolver.resolve_profile(profile, DEFAULT_WEIGHTS)

        has_conflict = len(profile.conflict_flags) > 0
        records.append(
            {
                "gene": gene,
                "scenario": scenario,
                "expected_conflict": expected_conflict[scenario],
                "detected_conflict": has_conflict,
                "final_score": float(profile.aggregate_score or 0.0),
                "flag_count": len(profile.conflict_flags),
                "layer1": float(profile.layer_scores.get("layer1_concordance", 0.0)),
                "penalty": float(profile.layer_scores.get("layer2_penalty", 0.0)),
                "layer3": float(profile.layer_scores.get("layer3_final", 0.0)),
            }
        )

    elapsed = time.perf_counter() - t0
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    total = len(records)
    detected = sum(1 for r in records if r["detected_conflict"])
    expected = sum(1 for r in records if r["expected_conflict"])
    correct = sum(1 for r in records if r["detected_conflict"] == r["expected_conflict"])

    by_scenario = {}
    for scenario in scenarios:
        subset = [r for r in records if r["scenario"] == scenario]
        by_scenario[scenario] = {
            "n": len(subset),
            "conflict_rate": (sum(1 for r in subset if r["detected_conflict"]) / len(subset)) if subset else 0.0,
            "mean_flag_count": statistics.mean([r["flag_count"] for r in subset]) if subset else 0.0,
            "mean_final_score": statistics.mean([r["final_score"] for r in subset]) if subset else 0.0,
            "mean_penalty": statistics.mean([r["penalty"] for r in subset]) if subset else 0.0,
        }

    return {
        "mode": "synthetic",
        "summary": {
            "n_targets": total,
            "elapsed_seconds": elapsed,
            "throughput_targets_per_second": total / elapsed if elapsed > 0 else 0.0,
            "peak_memory_mb": peak_mem / (1024 * 1024),
            "detected_conflicts": detected,
            "expected_conflicts": expected,
            "conflict_rate": detected / total if total else 0.0,
            "conflict_detection_accuracy": correct / total if total else 0.0,
        },
        "scenario_breakdown": by_scenario,
        "records": records,
    }


def build_real_profile(gene_symbol, use_ot_builtin=False):
    profile = TargetProfile(gene_symbol=gene_symbol.upper())
    try:
        populate_profile_from_pharos(profile)
    except Exception:
        pass
    try:
        populate_profile_from_depmap(profile, use_builtin=True)
    except Exception:
        pass
    try:
        populate_profile_from_opentargets(profile, use_builtin=use_ot_builtin)
    except Exception:
        pass
    return profile


def _clip(value):
    return max(0.0, min(1.0, float(value)))


def apply_synthetic_perturbation(profile, rng, noise_std=0.15, contradiction_rate=0.45):
    for evidence in profile.evidence:
        evidence.score = _clip(evidence.score + rng.gauss(0.0, noise_std))
        evidence.confidence_score = _clip((evidence.confidence_score or 0.5) * rng.uniform(0.85, 1.05))
        if evidence.provenance is None:
            evidence.provenance = {}
        evidence.provenance["freshness_days"] = float(evidence.provenance.get("freshness_days", 90.0)) * rng.uniform(1.0, 1.8)

    injected = rng.random() < contradiction_rate
    if injected:
        for evidence in profile.evidence:
            if evidence.source == EvidenceSource.DEPMAP and evidence.feature_name == "selectivity_score":
                evidence.value = rng.uniform(0.02, 0.12)
                evidence.score = float(evidence.value)
            if evidence.source == EvidenceSource.DEPMAP and evidence.feature_name == "common_essential":
                evidence.value = True
                evidence.score = 0.0
            if evidence.source == EvidenceSource.PHAROS and evidence.feature_name == "target_development_level":
                evidence.value = "Tclin"
                evidence.score = max(evidence.score, 0.85)
            if evidence.source == EvidenceSource.OPENTARGETS and evidence.feature_name == "ot_max_disease_score":
                evidence.value = rng.uniform(0.9, 0.99)
                evidence.score = float(evidence.value)

    return injected


def run_hybrid_benchmark(n_targets, seed, genes, use_ot_builtin=False, noise_std=0.15, contradiction_rate=0.45):
    rng = random.Random(seed)
    resolver = ConflictResolver()
    records = []

    tracemalloc.start()
    t0 = time.perf_counter()

    for i in range(n_targets):
        gene = genes[i % len(genes)].upper()
        baseline = build_real_profile(gene, use_ot_builtin=use_ot_builtin)
        if not baseline.evidence:
            continue

        baseline.compute_aggregate_score(DEFAULT_WEIGHTS)
        resolver.resolve_profile(baseline, DEFAULT_WEIGHTS)

        perturbed = copy.deepcopy(baseline)
        perturbed.gene_symbol = f"{gene}_NOISY_{i:04d}"
        injected = apply_synthetic_perturbation(
            perturbed,
            rng=rng,
            noise_std=noise_std,
            contradiction_rate=contradiction_rate,
        )
        perturbed.compute_aggregate_score(DEFAULT_WEIGHTS)
        resolver.resolve_profile(perturbed, DEFAULT_WEIGHTS)

        baseline_conflict = len(baseline.conflict_flags) > 0
        perturbed_conflict = len(perturbed.conflict_flags) > 0
        records.append(
            {
                "gene": gene,
                "baseline_score": float(baseline.aggregate_score or 0.0),
                "perturbed_score": float(perturbed.aggregate_score or 0.0),
                "score_drop": float((baseline.aggregate_score or 0.0) - (perturbed.aggregate_score or 0.0)),
                "baseline_conflict": baseline_conflict,
                "perturbed_conflict": perturbed_conflict,
                "baseline_flags": len(baseline.conflict_flags),
                "perturbed_flags": len(perturbed.conflict_flags),
                "injected_contradiction": injected,
                "detected_in_perturbed": perturbed_conflict,
                "perturbed_penalty": float(perturbed.layer_scores.get("layer2_penalty", 0.0)),
                "n_evidence": len(perturbed.evidence),
            }
        )

    elapsed = time.perf_counter() - t0
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    total = len(records)
    baseline_conflict_rate = (sum(1 for r in records if r["baseline_conflict"]) / total) if total else 0.0
    perturbed_conflict_rate = (sum(1 for r in records if r["perturbed_conflict"]) / total) if total else 0.0
    delta_conflict_rate = perturbed_conflict_rate - baseline_conflict_rate

    injected_total = sum(1 for r in records if r["injected_contradiction"])
    injected_detected = sum(1 for r in records if r["injected_contradiction"] and r["detected_in_perturbed"])
    non_injected_total = sum(1 for r in records if not r["injected_contradiction"])
    non_injected_clear = sum(1 for r in records if (not r["injected_contradiction"]) and (not r["detected_in_perturbed"]))

    expected_total = injected_total + non_injected_total
    correct_total = injected_detected + non_injected_clear

    by_scenario = {
        "baseline_real": {
            "n": total,
            "conflict_rate": baseline_conflict_rate,
            "mean_score": statistics.mean([r["baseline_score"] for r in records]) if records else 0.0,
        },
        "perturbed_real": {
            "n": total,
            "conflict_rate": perturbed_conflict_rate,
            "mean_score": statistics.mean([r["perturbed_score"] for r in records]) if records else 0.0,
            "mean_penalty": statistics.mean([r["perturbed_penalty"] for r in records]) if records else 0.0,
        },
        "injected_subset": {
            "n": injected_total,
            "conflict_rate": (
                sum(1 for r in records if r["injected_contradiction"] and r["detected_in_perturbed"]) / injected_total
            ) if injected_total else 0.0,
            "mean_score_drop": statistics.mean([r["score_drop"] for r in records if r["injected_contradiction"]]) if injected_total else 0.0,
        },
        "non_injected_subset": {
            "n": non_injected_total,
            "conflict_clear_rate": (
                sum(1 for r in records if (not r["injected_contradiction"]) and (not r["detected_in_perturbed"])) / non_injected_total
            ) if non_injected_total else 0.0,
            "mean_score_drop": statistics.mean([r["score_drop"] for r in records if not r["injected_contradiction"]]) if non_injected_total else 0.0,
        },
    }

    return {
        "mode": "hybrid",
        "summary": {
            "n_targets": total,
            "elapsed_seconds": elapsed,
            "throughput_targets_per_second": total / elapsed if elapsed > 0 else 0.0,
            "peak_memory_mb": peak_mem / (1024 * 1024),
            "baseline_conflict_rate": baseline_conflict_rate,
            "perturbed_conflict_rate": perturbed_conflict_rate,
            "delta_conflict_rate": delta_conflict_rate,
            "mean_score_drop": statistics.mean([r["score_drop"] for r in records]) if records else 0.0,
            "injected_cases": injected_total,
            "detected_in_injected": injected_detected,
            "hybrid_detection_accuracy": (correct_total / expected_total) if expected_total else 0.0,
            "noise_std": noise_std,
            "contradiction_rate": contradiction_rate,
            "opentargets_builtin": use_ot_builtin,
        },
        "scenario_breakdown": by_scenario,
        "records": records,
    }


def print_report(result):
    mode = result.get("mode", "synthetic")
    summary = result["summary"]
    print()
    print("=" * 92)
    print("Synthetic Conflict Stress Benchmark" if mode == "synthetic" else "Hybrid Robustness Benchmark")
    print("=" * 92)
    print(f"Mode:                       {mode}")
    print(f"Targets:                    {summary['n_targets']}")
    print(f"Elapsed time (s):           {summary['elapsed_seconds']:.4f}")
    print(f"Throughput (targets/s):     {summary['throughput_targets_per_second']:.2f}")
    print(f"Peak memory (MB):           {summary['peak_memory_mb']:.2f}")
    if mode == "synthetic":
        print(f"Detected conflicts:         {summary['detected_conflicts']}")
        print(f"Expected conflicts:         {summary['expected_conflicts']}")
        print(f"Conflict rate:              {summary['conflict_rate']:.2%}")
        print(f"Conflict detection accuracy:{summary['conflict_detection_accuracy']:.2%}")
        print()
        print(
            f"{'Scenario':<22}{'N':>8}{'ConflictRate':>16}"
            f"{'MeanFlags':>12}{'MeanPenalty':>14}{'MeanFinal':>12}"
        )
        print("-" * 92)
        for scenario, stats in result["scenario_breakdown"].items():
            print(
                f"{scenario:<22}{stats['n']:>8}{stats['conflict_rate']:>15.2%}"
                f"{stats['mean_flag_count']:>12.2f}{stats['mean_penalty']:>14.4f}{stats['mean_final_score']:>12.4f}"
            )
    else:
        print(f"Baseline conflict rate:     {summary['baseline_conflict_rate']:.2%}")
        print(f"Perturbed conflict rate:    {summary['perturbed_conflict_rate']:.2%}")
        print(f"Delta conflict rate:        {summary['delta_conflict_rate']:.2%}")
        print(f"Mean score drop:            {summary['mean_score_drop']:.4f}")
        print(f"Injected contradiction cases:{summary['injected_cases']}")
        print(f"Detected in injected cases: {summary['detected_in_injected']}")
        print(f"Hybrid detection accuracy:  {summary['hybrid_detection_accuracy']:.2%}")
        print(f"Noise std:                  {summary['noise_std']:.3f}")
        print(f"Contradiction rate:         {summary['contradiction_rate']:.2f}")
        print(f"Open Targets built-in mode: {summary['opentargets_builtin']}")
        print()
        print(f"{'Scenario':<22}{'N':>8}{'Rate/Value':>18}{'MeanScoreDrop':>18}")
        print("-" * 92)
        baseline = result["scenario_breakdown"]["baseline_real"]
        perturbed = result["scenario_breakdown"]["perturbed_real"]
        injected = result["scenario_breakdown"]["injected_subset"]
        clean = result["scenario_breakdown"]["non_injected_subset"]
        print(f"{'baseline_real':<22}{baseline['n']:>8}{baseline['conflict_rate']:>17.2%}{0.0:>18.4f}")
        print(f"{'perturbed_real':<22}{perturbed['n']:>8}{perturbed['conflict_rate']:>17.2%}{0.0:>18.4f}")
        print(f"{'injected_subset':<22}{injected['n']:>8}{injected['conflict_rate']:>17.2%}{injected['mean_score_drop']:>18.4f}")
        print(f"{'non_injected_subset':<22}{clean['n']:>8}{clean['conflict_clear_rate']:>17.2%}{clean['mean_score_drop']:>18.4f}")
    print("=" * 92)


def main():
    parser = argparse.ArgumentParser(description="Benchmark conflict handling and robustness under noisy evidence")
    parser.add_argument(
        "--mode",
        choices=["synthetic", "hybrid"],
        default="synthetic",
        help="synthetic: fully synthetic stress test; hybrid: real profiles + synthetic perturbations",
    )
    parser.add_argument("--n-targets", type=int, default=300, help="Number of synthetic targets")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--output", type=str, default="", help="Optional JSON output path")
    parser.add_argument(
        "--real-genes",
        nargs="+",
        default=["EGFR", "KRAS", "TP53", "BRCA1", "BRAF", "MYC", "CDK4", "MTOR"],
        help="Gene set used by hybrid mode",
    )
    parser.add_argument(
        "--hybrid-opentargets-builtin",
        action="store_true",
        help="Use Open Targets built-in fallback in hybrid mode",
    )
    parser.add_argument(
        "--noise-std",
        type=float,
        default=0.15,
        help="Standard deviation of score noise in hybrid mode",
    )
    parser.add_argument(
        "--contradiction-rate",
        type=float,
        default=0.45,
        help="Probability of injecting targeted contradictions in hybrid mode",
    )
    args = parser.parse_args()

    if args.mode == "synthetic":
        result = run_benchmark(n_targets=args.n_targets, seed=args.seed)
    else:
        result = run_hybrid_benchmark(
            n_targets=args.n_targets,
            seed=args.seed,
            genes=args.real_genes,
            use_ot_builtin=args.hybrid_opentargets_builtin,
            noise_std=args.noise_std,
            contradiction_rate=args.contradiction_rate,
        )

    print_report(result)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Saved JSON report to: {args.output}")


if __name__ == "__main__":
    main()
