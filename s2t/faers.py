"""Download and process FAERS quarterly ASCII files (FDA), with FDA-recommended case deduplication.

Usage:
    python -m s2t.faers --quarters 2026q2 2026q1 2025q4 2025q3 --drugs drugs.txt --out data/faers_processed.pkl
"""
import argparse, io, pickle, re, subprocess, zipfile
from pathlib import Path
import pandas as pd

FDA = "https://fis.fda.gov/content/Exports/faers_ascii_{q}.zip"

def download(q: str, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True); p = dest / f"{q}.zip"
    if not p.exists():
        subprocess.run(["curl", "-sSL", "-o", str(p), FDA.format(q=q)], check=True)
    return p

def read_table(z: zipfile.ZipFile, prefix: str, cols):
    name = [x for x in z.namelist() if re.search(rf"(?i)ascii/{prefix}\d+q\d\.txt$", x)][0]
    return pd.read_csv(z.open(name), sep="$", usecols=cols, dtype=str, encoding="latin-1", on_bad_lines="skip")

def load_quarters(quarters, dest: Path):
    drug, reac, demo = [], [], []
    for q in quarters:
        z = zipfile.ZipFile(download(q, dest))
        drug.append(read_table(z, "DRUG", ["primaryid", "caseid", "role_cod", "drugname", "prod_ai"]))
        reac.append(read_table(z, "REAC", ["primaryid", "caseid", "pt"]))
        demo.append(read_table(z, "DEMO", ["primaryid", "caseid", "fda_dt"]))
    drug, reac, demo = (pd.concat(x) for x in (drug, reac, demo))
    # FDA guidance: one version per case, keep the most recent (highest fda_dt, then primaryid)
    demo["fda_dt"] = demo["fda_dt"].fillna("0")
    keep = set(demo.sort_values(["caseid", "fda_dt", "primaryid"]).drop_duplicates("caseid", keep="last")["primaryid"])
    drug = drug[drug.primaryid.isin(keep) & drug.role_cod.isin(["PS", "SS"])].copy()   # suspect drugs only
    reac = reac[reac.primaryid.isin(keep)].copy()
    drug["ai"] = drug["prod_ai"].fillna(drug["drugname"]).astype(str).str.upper()
    reac["pt"] = reac["pt"].str.strip().str.upper()
    reac = reac.drop_duplicates(["primaryid", "pt"])
    return drug, reac, len(keep)

def process(quarters, drugs, dest: Path):
    drug, reac, n_cases = load_quarters(quarters, dest)
    pt_counts = reac.groupby("pt").primaryid.nunique()
    drug_cases = {d: set(drug.loc[drug["ai"].str.contains(d.upper(), regex=False), "primaryid"]) for d in drugs}
    reac_by_pid = reac.groupby("primaryid").pt.apply(list)
    return {"N_total": n_cases, "pt_counts": pt_counts, "drug_cases": drug_cases, "reac_by_pid": reac_by_pid, "quarters": list(quarters)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--quarters", nargs="+", required=True); ap.add_argument("--drugs", required=True)
    ap.add_argument("--dest", default="data/faers"); ap.add_argument("--out", default="data/faers_processed.pkl")
    a = ap.parse_args(); drugs = [l.strip() for l in open(a.drugs) if l.strip()]
    out = process(a.quarters, drugs, Path(a.dest)); pickle.dump(out, open(a.out, "wb"))
    print(f"{out['N_total']} deduplicated cases, {len(out['pt_counts'])} PTs, {sum(1 for v in out['drug_cases'].values() if v)} drugs matched")
