"""Fetch drug bioactivities (pChEMBL) against human single-protein targets from ChEMBL.

Usage:
    python -m s2t.chembl --drugs drugs.txt --out data/drug_target_pchembl.csv
"""
import argparse, json, time
import pandas as pd, requests

CH = "https://www.ebi.ac.uk/chembl/api/data"
S = requests.Session(); S.headers["User-Agent"] = "signal-to-target/0.1 (research)"

def molecule(name):
    r = S.get(f"{CH}/molecule.json", params={"pref_name__iexact": name.upper(), "only": "molecule_chembl_id,pref_name,max_phase"}, timeout=40); r.raise_for_status()
    ms = r.json()["molecules"]
    if not ms:
        ms = S.get(f"{CH}/molecule/search.json", params={"q": name, "limit": 5}, timeout=40).json().get("molecules", [])
    def ph(m):
        try: return float(m.get("max_phase") or 0)   # ChEMBL returns max_phase as a string
        except ValueError: return 0.0
    ms = [m for m in ms if ph(m) >= 3] or ms
    return ms[0] if ms else None

def activities(chembl_id, pchembl_min=5.0, max_rows=4000):
    acts, url = [], f"{CH}/activity.json"
    params = {"molecule_chembl_id": chembl_id, "pchembl_value__gte": pchembl_min, "target_organism__iexact": "Homo sapiens",
              "limit": 1000, "only": "target_chembl_id,target_pref_name,pchembl_value,standard_type,assay_type"}
    while url and len(acts) < max_rows:
        j = S.get(url, params=params, timeout=60).json(); acts += j["activities"]
        nxt = j["page_meta"].get("next"); url = ("https://www.ebi.ac.uk" + nxt) if nxt else None; params = None
    return acts

def targets(tids):
    info = {}
    for i in range(0, len(tids), 40):
        chunk = tids[i:i + 40]
        r = S.get(f"{CH}/target.json", params={"target_chembl_id__in": ",".join(chunk), "limit": 40, "only": "target_chembl_id,pref_name,target_type,target_components"}, timeout=60)
        for t in r.json()["targets"]:
            comps = t.get("target_components") or []; acc = gene = None
            if t["target_type"] == "SINGLE PROTEIN" and comps:
                acc = comps[0].get("accession")
                gene = next((s["component_synonym"] for s in comps[0].get("target_component_synonyms", []) if s.get("syn_type") == "GENE_SYMBOL"), None)
            info[t["target_chembl_id"]] = {"name": t["pref_name"], "type": t["target_type"], "acc": acc, "gene": gene}
    return info

def build(drugs):
    rows, raw = [], {}
    for d in drugs:
        m = molecule(d)
        if not m: continue
        raw[d] = {"chembl": m, "activities": activities(m["molecule_chembl_id"])}
    tinfo = targets(sorted({a["target_chembl_id"] for v in raw.values() for a in v["activities"]}))
    for d, v in raw.items():
        for a in v["activities"]:
            t = tinfo.get(a["target_chembl_id"])
            if t and t["type"] == "SINGLE PROTEIN" and t["acc"] and a.get("pchembl_value") and a.get("assay_type") in ("B", "F"):
                rows.append((d, a["target_chembl_id"], t["gene"] or t["name"], t["acc"], float(a["pchembl_value"])))
    df = pd.DataFrame(rows, columns=["drug", "tid", "gene", "acc", "pchembl"]).groupby(["drug", "tid", "gene", "acc"]).pchembl.max().reset_index()
    return df, raw, tinfo

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--drugs", required=True); ap.add_argument("--out", default="data/drug_target_pchembl.csv")
    a = ap.parse_args(); drugs = [l.strip() for l in open(a.drugs) if l.strip()]
    df, raw, tinfo = build(drugs); df.to_csv(a.out, index=False); print(len(df), "drug-target pairs")
