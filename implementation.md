# Phase 77 — Claude Reviewer for HTS Top-20 Hits

**Version:** 1.1 | **Tier:** Micro | **Date:** 2026-03-28

## Goal
Use Claude to review and triage the top-ranked HTS hits. Given ranked compounds with activity data, Claude provides structured assessments of drug-likeness, SAR context, and development priority.

CLI: `python main.py --input data/compounds.csv --n 3`

Outputs: hts_reviews.json, hts_review_report.txt

## Logic
1. Load compounds, sort by pIC50 descending (simulates HTS ranking)
2. Select top-N hits for review
3. For each hit, compute RDKit properties (MW, LogP, HBA, HBD, TPSA, rotatable bonds, RO5 violations)
4. Send compound profile to Claude for structured triage
5. Parse: priority (high/medium/low), liability_flags, rationale, next_step
6. Save structured reviews

## Key Concepts
- Integration pattern: RDKit computation + Claude review
- HTS hit triage workflow (rank -> compute -> assess -> prioritize)
- Structured priority assignment with explicit criteria
- Combining computed properties with AI domain assessment

## Verification Checklist
- [x] Top-N compounds ranked by pIC50
- [x] RDKit properties computed per compound (7 descriptors)
- [x] Claude review with priority, flags, rationale, next_step
- [x] JSON output parseable
- [x] --help works

## Results
| Rank | Compound | pIC50 | Priority | Flags |
|------|----------|-------|----------|-------|
| 1 | ind_006_CF3 | 8.55 | high | none |
| 2 | ind_007_CN | 8.35 | high | none |
| 3 | quin_006_CF3 | 8.25 | high | high_logp |

| Metric | Value |
|--------|-------|
| Priority distribution | high: 3 |
| Total tokens | in=958 out=386 |
| Est. cost | $0.0023 |

Key findings:
- All top-3 are pIC50 >= 8.0, so all correctly triaged as HIGH priority
- Claude identified high_logp flag for quin_006_CF3 (logP=3.25)
- Actionable next steps provided: ADME profiling, SAR around CF3, selectivity
- Integration pattern works: RDKit computes, Claude interprets

## Risks
- All top-3 were genuinely high priority — need lower-ranked compounds to test medium/low assignment
- rdkit-pypi not available on Python 3.13; use `rdkit` package instead
