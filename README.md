# Target Identification Pipeline

A simple prototype that predicts how **druggable** a gene/protein target is by collecting evidence from multiple biological databases and combining them into a single score.

**Based on:** [DrugnomeAI](https://doi.org/10.1038/s42003-022-04245-4) by Raies et al., 2022

---

## What Does This Do?

In drug discovery, before you can make a drug, you need to find the right **target** — usually a protein in the body that, when blocked or activated, can treat a disease. But how do you know if a target is actually "druggable"?

This pipeline answers that by:
1. Pulling data about genes from **real biological databases**
2. Scoring each piece of evidence on a 0→1 scale
3. Combining all scores with a **weighted average** to rank genes

---

## How It Works

Here's the overall flow of the pipeline:

```
    ┌──────────────────────┐
    │   User gives a list  │
    │   of gene symbols    │
    │  (e.g. EGFR, KRAS)   │
    └──────────┬───────────┘
               │
               ▼
    ┌──────────────────────┐
    │     PIPELINE         │
    │  (pipeline.py)       │
    │                      │
    │  For each gene:      │
    │  ┌────────────────┐  │
    │  │ 1. Query PHAROS│──│──► PHAROS API (internet)
    │  └────────────────┘  │     returns: TDL, novelty,
    │  ┌────────────────┐  │     target family, description
    │  │ 2. Query DepMap│──│──► Built-in dataset
    │  └────────────────┘  │     returns: CRISPR scores,
    │  ┌────────────────┐  │     essentiality, selectivity
    │  │ 3. Score &     │  │
    │  │    Aggregate   │  │
    │  └────────────────┘  │
    └──────────┬───────────┘
               │
               ▼
    ┌──────────────────────┐
    │  Ranked list of      │
    │  genes by drug-      │
    │  gability score      │
    │  + JSON export       │
    └──────────────────────┘
```

---

## Project Structure

```
target_id_prototype/
│
├── evidence_schema.py     ← Data models (EvidenceItem, TargetProfile, enums)
├── source_pharos.py       ← Fetches data from PHAROS GraphQL API
├── source_depmap.py       ← Fetches data from DepMap (built-in dataset)
├── pipeline.py            ← Main pipeline that ties everything together
├── run_example.py         ← CLI script to run the pipeline
├── requirements.txt       ← Python dependencies
└── README.md              ← You are here!
```

---

## Data Sources

### 1. PHAROS (Illuminating the Druggable Genome)

PHAROS classifies every human protein into one of four **Target Development Levels**:

```
    Most studied                              Least studied
    ◄──────────────────────────────────────────────────────►

    ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌─────────┐
    │  Tclin  │  │  Tchem  │  │  Tbio   │  │  Tdark  │
    │ Score:  │  │ Score:  │  │ Score:  │  │ Score:  │
    │  1.00   │  │  0.75   │  │  0.40   │  │  0.10   │
    └─────────┘  └─────────┘  └─────────┘  └─────────┘
    Has approved   Has potent   Known        Barely
    drugs          compounds    biology      studied
```

**What we extract from PHAROS:**
| Feature | Category | What it means |
|---------|----------|---------------|
| Target Development Level | Druggability | How far along is this target in drug development? |
| Novelty score | Druggability | How novel/understudied is this target? |
| Target family | Protein Structure | Is it a Kinase, GPCR, Enzyme, etc.? |
| Function description | Phenotype | Curated text about what the protein does |

### DepMap (Cancer Dependency Map)

DepMap knocks out genes in ~1000 cancer cell lines using CRISPR and measures which ones die. This tells us which genes are **essential** for cancer survival.

```
    CRISPR Dependency Score
    ─────────────────────────────────────────────

     0.0                    -0.5                -1.0
      │                       │                   │
      ▼                       ▼                   ▼
    Not essential        Somewhat             Highly essential
    (gene knockout       essential            (cells die without
     has no effect)                            this gene)
```

**What we extract from DepMap:**
| Feature | What it means |
|---------|---------------|
| Mean CRISPR score | Average essentiality across all cell lines |
| Dependent cell line fraction | What % of cell lines need this gene to survive? |
| Selectivity score | Is it essential in only SOME cell lines? (good for drugs) |
| Common essential flag | Is it essential everywhere? (bad — would be toxic) |

---

## Scoring System

Each piece of evidence gets a score between **0.0** (not druggable) and **1.0** (very druggable).

To combine scores across categories, we use a **weighted average**:

```
                        Σ (weight_i × mean_score_i)
    Final Score  =  ─────────────────────────────────
                           Σ (weight_i)
```

The weights reflect how important each type of evidence is:

| Category | Weight | Why? |
|----------|--------|------|
| Druggability | 2.0 | Direct indicator — most important |
| Essentiality | 1.5 | Cancer dependency data is very relevant |
| Disease Association | 1.2 | More disease links = more useful target |
| Protein Interaction | 1.0 | Network context matters |
| Protein Structure | 1.0 | Structural info helps drug design |
| Genetic Constraint | 1.0 | Evolutionary importance |
| Expression | 0.8 | Where the gene is active |
| Pathway | 0.8 | What biological pathways it's in |
| Phenotype | 0.5 | Descriptive — less directly useful |

---

## How to Run

### Install

```bash
pip install -r requirements.txt
```

### Run

```bash
# Default run (PHAROS API + DepMap built-in data)
python3 run_example.py

# Offline mode (no internet needed)
python3 run_example.py --offline

# Pick specific genes
python3 run_example.py --genes EGFR BRAF PIK3CA

# Save results to a specific file
python3 run_example.py --output my_results.json
```

### Example Output

```
  DRUGGABILITY RANKING
═══════════════════════════════════════════════════════════
Rank  Gene       TDL      Score   #Evidence Sources
───────────────────────────────────────────────────────────
1     EGFR       Tclin    0.5838          8 depmap, pharos
2     BRAF       Tclin    0.5741          8 depmap, pharos
3     KRAS       Tclin    0.5565          8 depmap, pharos
4     BRCA1      Tchem    0.4410          8 depmap, pharos
5     TP53       Tchem    0.4075          8 depmap, pharos
```

---

## How to Add a New Data Source

The pipeline is designed to be extended. To add a new source (e.g., Open Targets):

**Step 1:** Create `source_opentargets.py`

```python
from evidence_schema import EvidenceItem, EvidenceSource, EvidenceCategory, ConfidenceLevel, TargetProfile

def fetch_opentargets_evidence(gene_symbol):
    items = []
    # ... call the Open Targets API ...
    # ... create EvidenceItem objects ...
    return items

def populate_profile_from_opentargets(profile):
    items = fetch_opentargets_evidence(profile.gene_symbol)
    for item in items:
        profile.add_evidence(item)
    return profile
```

**Step 2:** Add it to `pipeline.py`

```python
from source_opentargets import populate_profile_from_opentargets

# Inside the run() method, add:
if self.use_opentargets:
    populate_profile_from_opentargets(profile)
```

That's it! The schema handles deduplication and scoring automatically.

---

## How This Relates to DrugnomeAI

This prototype simplifies DrugnomeAI's approach to make it easier to understand:

```
    DrugnomeAI (full paper)          This Prototype (simplified)
    ──────────────────────           ──────────────────────────
    324+ features                    ~8 features
    20+ databases                    2 databases (PHAROS + DepMap)
    12 ML classifiers                Weighted average
    Genome-wide (20K genes)          User-specified gene list
    Trained on known drug targets    Score-based (no training)
```

The key idea from the paper that we keep is: **combining evidence from multiple independent sources gives better predictions than any single source alone.**

---

## References

1. Raies, A. et al. "DrugnomeAI is an ensemble machine-learning framework for predicting druggability of candidate drug targets." *Communications Biology* 5, 1291 (2022). [Link](https://doi.org/10.1038/s42003-022-04245-4)

2. Sheils, T.K. et al. "TCRD and Pharos 2021: mining the human proteome for disease biology." *Nucleic Acids Research* 49, D1334–D1346 (2021).

3. Tsherniak, A. et al. "Defining a cancer dependency map." *Cell* 170, 564–576 (2017).
