"""Disproportionality (ROR, BCPNN IC) per drug-event, then target-event enrichment across drugs (Fisher, BH-FDR),
plus recovery of known safety-target associations as positive controls.

Usage:
    python -m s2t.analyze --faers data/faers_processed.pkl --targets data/drug_target_pchembl.csv --out results/
"""
import argparse, math, pickle
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

KNOWN = {  # curated from Bowes et al. 2012 (Nat Rev Drug Discov) safety panel + standard pharmacology; MedDRA PTs
 "KCNH2": ["ELECTROCARDIOGRAM QT PROLONGED","TORSADE DE POINTES","VENTRICULAR TACHYCARDIA","LONG QT SYNDROME","SUDDEN CARDIAC DEATH"],
 "HTR2B": ["CARDIAC VALVE DISEASE","MITRAL VALVE INCOMPETENCE","AORTIC VALVE INCOMPETENCE","TRICUSPID VALVE INCOMPETENCE","PULMONARY ARTERIAL HYPERTENSION"],
 "CHRM3": ["DRY MOUTH","URINARY RETENTION","CONSTIPATION","VISION BLURRED","MYDRIASIS"],
 "CHRM1": ["DRY MOUTH","URINARY RETENTION","CONSTIPATION","CONFUSIONAL STATE","MEMORY IMPAIRMENT","DELIRIUM"],
 "HRH1": ["SOMNOLENCE","SEDATION","WEIGHT INCREASED","FATIGUE"],
 "ADRA1A": ["ORTHOSTATIC HYPOTENSION","DIZZINESS POSTURAL","HYPOTENSION","SYNCOPE"],
 "DRD2": ["EXTRAPYRAMIDAL DISORDER","HYPERPROLACTINAEMIA","TARDIVE DYSKINESIA","DYSTONIA","AKATHISIA","GALACTORRHOEA","PARKINSONISM","NEUROLEPTIC MALIGNANT SYNDROME"],
 "SLC6A4": ["SEROTONIN SYNDROME","HYPONATRAEMIA","NAUSEA","SEXUAL DYSFUNCTION","LIBIDO DECREASED","INAPPROPRIATE ANTIDIURETIC HORMONE SECRETION","BRUXISM"],
 "SLC6A2": ["TACHYCARDIA","HYPERTENSION","HEART RATE INCREASED","PALPITATIONS","URINARY RETENTION"],
 "OPRM1": ["CONSTIPATION","RESPIRATORY DEPRESSION","MIOSIS","DRUG DEPENDENCE","PRURITUS","OPIOID USE DISORDER"],
 "PTGS1": ["GASTROINTESTINAL HAEMORRHAGE","GASTRIC ULCER","MELAENA","GASTROINTESTINAL ULCER","HAEMATEMESIS"],
 "PTGS2": ["GASTROINTESTINAL HAEMORRHAGE","GASTRIC ULCER","RENAL FAILURE","ACUTE KIDNEY INJURY","OEDEMA PERIPHERAL"],
 "ADRB1": ["BRADYCARDIA","HYPOTENSION","ATRIOVENTRICULAR BLOCK"],
 "ADRB2": ["TREMOR","BRONCHOSPASM","HYPOKALAEMIA"],
 "HTR2A": ["SEXUAL DYSFUNCTION","SOMNOLENCE","WEIGHT INCREASED"],
 "HTR2C": ["WEIGHT INCREASED","INCREASED APPETITE","HYPERGLYCAEMIA","DIABETES MELLITUS"],
 "GABRA1": ["SOMNOLENCE","SEDATION","AMNESIA","DRUG DEPENDENCE","RESPIRATORY DEPRESSION"],
 "SCN5A": ["ATRIOVENTRICULAR BLOCK","BRADYCARDIA","ELECTROCARDIOGRAM QRS COMPLEX PROLONGED"],
 "PDE5A": ["VISUAL IMPAIRMENT","CHROMATOPSIA","HEADACHE","FLUSHING","HEARING IMPAIRED","PRIAPISM"],
 "CACNA1C": ["OEDEMA PERIPHERAL","HYPOTENSION","GINGIVAL HYPERTROPHY"],
}

def disproportionality(P, drugs, min_cases=100, min_a=3):
    N, pt_counts, drug_cases, reac_by_pid = P["N_total"], P["pt_counts"], P["drug_cases"], P["reac_by_pid"]
    rows = []
    for d in drugs:
        pids = drug_cases.get(d, set()); n_drug = len(pids)
        if n_drug < min_cases: continue
        cnt = {}
        for p in pids:
            for pt in reac_by_pid.get(p, []): cnt[pt] = cnt.get(pt, 0) + 1
        for pt, a in cnt.items():
            if a < min_a: continue
            n_pt = int(pt_counts[pt]); b = n_drug - a; c = n_pt - a; dd = N - n_drug - c
            if b == 0 or c == 0: continue
            ror = (a * dd) / (b * c); se = math.sqrt(1/a + 1/b + 1/c + 1/dd)
            lo, hi = math.exp(math.log(ror) - 1.96*se), math.exp(math.log(ror) + 1.96*se)
            E = n_drug * n_pt / N; ic = math.log2((a + 0.5) / (E + 0.5))
            ic025 = ic - 3.3 * (a + 0.5) ** -0.5 - 2 * (a + 0.5) ** -1.5          # Norén et al. 2013 approximation
            rows.append((d, pt, a, n_drug, n_pt, ror, lo, hi, ic, ic025))
    sig = pd.DataFrame(rows, columns=["drug","pt","a","n_drug","n_pt","ror","ror_lo","ror_hi","ic","ic025"])
    sig["signal"] = (sig.a >= min_a) & (sig.ror_lo > 1) & (sig.ic025 > 0)
    return sig

def target_event(sig, act, study_drugs, pchembl_active=6.0, min_active=4, min_inactive=10, min_signal_drugs=3):
    A = act[act.drug.isin(study_drugs)]; SD = set(study_drugs)
    signals_by_pt = sig[sig.signal].groupby("pt").drug.apply(set)
    rows = []
    for g in A.gene.unique():
        act_drugs = set(A[(A.gene == g) & (A.pchembl >= pchembl_active)].drug); inactive = SD - act_drugs
        if len(act_drugs) < min_active or len(inactive) < min_inactive: continue
        for pt, sdr in signals_by_pt.items():
            if len(sdr) < min_signal_drugs: continue
            a = len(act_drugs & sdr); b = len(act_drugs) - a; c = len(inactive & sdr); dd = len(inactive) - c
            if a < 2: continue
            orr, p = stats.fisher_exact([[a, b], [c, dd]], alternative="greater")
            rows.append((g, pt, a, len(act_drugs), c, len(inactive), orr, p, sorted(act_drugs & sdr)))
    te = pd.DataFrame(rows, columns=["gene","pt","active_with_signal","n_active","inactive_with_signal","n_inactive","odds_ratio","p","drugs"])
    te["q"] = multipletests(te.p, method="fdr_bh")[1]
    return te.sort_values("p").reset_index(drop=True)

def recovery(te):
    out = []
    for g, pts in KNOWN.items():
        sub = te[te.gene == g]; hit = sub[sub.pt.isin(pts)].sort_values("p").head(1)
        if sub.empty: st = "not testable"
        elif hit.empty: st = "missed"
        else: st = "recovered (q<0.1)" if hit.q.iloc[0] < 0.1 else ("nominal (p<0.05)" if hit.p.iloc[0] < 0.05 else "missed")
        out.append((g, st, None if hit.empty else hit.pt.iloc[0], None if hit.empty else float(hit.p.iloc[0]), None if hit.empty else float(hit.q.iloc[0])))
    return pd.DataFrame(out, columns=["gene","status","best_pt","p","q"])

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--faers", default="data/faers_processed.pkl"); ap.add_argument("--targets", default="data/drug_target_pchembl.csv"); ap.add_argument("--out", default="results")
    a = ap.parse_args(); out = Path(a.out); out.mkdir(exist_ok=True)
    P = pickle.load(open(a.faers, "rb")); act = pd.read_csv(a.targets)
    drugs = sorted(P["drug_cases"]); sig = disproportionality(P, drugs); study = sorted(sig.drug.unique())
    te = target_event(sig, act, study); rec = recovery(te)
    sig.to_csv(out/"drug_event_disproportionality.csv", index=False); te.to_csv(out/"target_event_enrichment.csv", index=False); rec.to_csv(out/"positive_control_recovery.csv", index=False)
    print(f"{len(study)} drugs, {int(sig.signal.sum())} drug-event signals, {len(te)} target-event tests, {(te.q<0.05).sum()} at q<0.05"); print(rec.to_string(index=False))
