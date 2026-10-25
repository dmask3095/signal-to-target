#!/usr/bin/env python
# /// script
# dependencies = ["pandas", "scipy", "statsmodels"]
# ///
"""Print disproportionality signals for one drug from a processed FAERS pickle."""
import argparse, os, pickle, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from s2t.analyze import disproportionality
ap = argparse.ArgumentParser(); ap.add_argument("--drug", required=True); ap.add_argument("--faers", default="data/faers_processed.pkl"); ap.add_argument("--top", type=int, default=30)
a = ap.parse_args(); P = pickle.load(open(a.faers, "rb"))
sig = disproportionality(P, [a.drug], min_cases=1)
print(f"quarters={P.get('quarters')} dedup_cases={P['N_total']} drug_cases={len(P['drug_cases'][a.drug])}")
print(sig[sig.signal].sort_values("ror", ascending=False).head(a.top)[["pt","a","ror","ror_lo","ror_hi","ic","ic025"]].to_string(index=False))
