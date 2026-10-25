---
name: faers-disproportionality
description: Compute pharmacovigilance disproportionality statistics (reporting odds ratio with 95% CI, BCPNN information component IC/IC025) for drug–adverse-event pairs from FDA FAERS quarterly files, with FDA-recommended case deduplication and suspect-drug filtering. Use when the user asks whether an adverse event is reported disproportionately for a drug, wants a safety-signal table for a drug or drug class, or needs a reproducible FAERS analysis. Pair with chembl-database, uniprot-database and human-protein-atlas-database skills to explain signals mechanistically.
---

# FAERS disproportionality

## When to use
- "Is <event> a safety signal for <drug>?" / "Top signals for <drug>" / "Compare signals across <drug class>".
- Any analysis that needs counts of FAERS cases by drug and MedDRA preferred term with a full-database background.

## Workflow
1. Decide the quarters (default: the latest four available at https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html).
2. `uv run scripts/faers_prepare.py --quarters 2026q2 2026q1 2025q4 2025q3 --drugs <file>` downloads the ASCII zips, deduplicates
   cases (one version per `caseid`, most recent `fda_dt`), keeps primary/secondary suspect drugs, matches drugs on `prod_ai`.
3. `uv run scripts/faers_signals.py --drug <name>` prints ROR (95% CI), IC, IC025, case counts per PT, flagging signals with
   n>=3, ROR lower bound > 1 and IC025 > 0 (standard thresholds; state them in the report).
4. Always report: quarters used, number of deduplicated cases, that spontaneous reports cannot establish incidence or causality,
   and that signals are hypotheses for pharmacoepidemiological follow-up.

## Interpretation guardrails
- Disproportionality is relative to all other reports; indication, co-medication, stimulated reporting and notoriety bias inflate ROR.
- Compare suspect-only vs any-role counts if a signal looks implausible.
- For mechanism hypotheses, join to ChEMBL bioactivities (pChEMBL >= 6 ≈ 1 µM) and test target-level enrichment across drugs,
  then check tissue expression (Human Protein Atlas) and target function (UniProt). Many target-level enrichments are confounded by
  drug class; prefer associations supported by >= 3 pharmacological classes.

## References
- FDA FAERS quarterly data files and ASC_NTS readme; Norén et al. 2013 (IC shrinkage); Bate & Evans 2009 (disproportionality overview);
  Maciejewski et al. 2017 eLife (reverse translation of FAERS); Smit et al. 2021 (target–AE associations with plasma exposure).
