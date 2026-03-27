# Phase 77 — Claude Reviewer for HTS Top-20 Hits

**Version:** 1.0 | **Tier:** Micro | **Date:** 2026-03-28

## Goal
Use Claude to review and triage the top-ranked HTS hits. Given ranked compounds with activity data, Claude provides structured assessments of drug-likeness, SAR context, and development priority.

CLI: `python main.py --input data/compounds.csv --n 3`

Outputs: hts_reviews.json, hts_review_report.txt

## Logic
1. Load compounds, sort by pIC50 descending (simulates HTS ranking)
2. Select top-N hits for review
3. For each hit, compute RDKit properties (MW, LogP, HBA, HBD, TPSA)
4. Send compound profile to Claude for structured triage
5. Parse: priority (high/medium/low), liability flags, rationale, next_step
6. Save structured reviews

## Key Concepts
- Integration pattern: RDKit computation + Claude review
- HTS hit triage workflow
- Structured priority assignment
- Combining computed properties with AI assessment

## Verification Checklist
- [ ] Top-N compounds ranked by pIC50
- [ ] RDKit properties computed per compound
- [ ] Claude review with priority, flags, rationale, next_step
- [ ] JSON output parseable
- [ ] --help works

## Risks
- Claude may assign all "high" priority — mitigate with explicit criteria
- RDKit import may fail if not in venv
