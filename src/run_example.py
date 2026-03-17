#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline import TargetPipeline


def main():
    parser = argparse.ArgumentParser(description="Target Identification Pipeline")
    parser.add_argument("--genes", nargs="+", default=["EGFR", "KRAS", "TP53", "BRCA1", "MYC", "BRAF", "CDK4", "MTOR"])
    parser.add_argument("--offline", action="store_true", help="Skip PHAROS and Open Targets API calls")
    parser.add_argument("--no-opentargets", action="store_true", help="Disable Open Targets evidence source")
    parser.add_argument("--opentargets-builtin", action="store_true", help="Use built-in Open Targets sample data")
    parser.add_argument("--output", type=str, default=None, help="Path to save JSON results")
    parser.add_argument("--human-review", type=str, default="human_review_overrides.json", help="Path to human-in-the-loop override JSON")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s │ %(levelname)-7s │ %(name)s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    print()
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   Target Identification Pipeline – Minimal Prototype        ║")
    print("║   Inspired by DrugnomeAI (Raies et al., Commun. Biol. 2022)║")
    print("║   DOI: 10.1038/s42003-022-04245-4                          ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()

    pipeline = TargetPipeline(
        use_pharos=not args.offline,
        use_depmap=True,
        use_opentargets=(not args.no_opentargets) and (not args.offline),
        depmap_builtin=True,
        opentargets_builtin=args.opentargets_builtin,
        human_override_path=args.human_review,
    )

    gene_list = [g.upper() for g in args.genes]
    print(f"Evaluating {len(gene_list)} gene targets: {', '.join(gene_list)}")
    print(f"PHAROS API: {'DISABLED (offline mode)' if args.offline else 'ENABLED'}")
    print(f"Open Targets: {'DISABLED' if (args.no_opentargets or args.offline) else ('ENABLED (built-in)' if args.opentargets_builtin else 'ENABLED (live)')}")
    print(f"DepMap:      ENABLED (built-in demo data)")
    print(f"Human review overrides: {args.human_review}")
    print()

    results = pipeline.run(gene_list)

    print("\n" + "═" * 75)
    print("  DRUGGABILITY RANKING")
    print("═" * 75)
    print(pipeline.ranking_table())
    print()

    print("\n" + "═" * 75)
    print("  DETAILED EVIDENCE BREAKDOWN (Top 3)")
    print("═" * 75)
    for profile in results[:3]:
        print(pipeline.evidence_breakdown(profile.gene_symbol))

    output_path = args.output or os.path.join(os.path.dirname(os.path.abspath(__file__)), "results.json")
    with open(output_path, "w") as f:
        f.write(pipeline.to_json())
    print(f"Full results saved to: {output_path}")
    print()

    total_evidence = sum(len(p.evidence) for p in results)
    sources_used = set()
    for p in results:
        for e in p.evidence:
            sources_used.add(e.source.value)

    print("─" * 40)
    print(f"  Genes evaluated:     {len(results)}")
    print(f"  Total evidence items:{total_evidence}")
    print(f"  Data sources used:   {', '.join(sorted(sources_used))}")
    print(f"  Total conflict flags:{sum(len(p.conflict_flags) for p in results)}")
    print(f"  Top-ranked target:   {results[0].gene_symbol} (score={results[0].aggregate_score:.4f})")
    print("─" * 40)
    print()


if __name__ == "__main__":
    main()
