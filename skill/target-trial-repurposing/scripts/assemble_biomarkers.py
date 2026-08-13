#!/usr/bin/env python3
"""
assemble_biomarkers.py — organize biomarker outcomes (from extract_biomarkers.py) BY TARGET.

Writes two CSVs into --out-dir:
  1) <prefix>Biomarker_measures.csv   granular: one row per target x trial x measure (sorted by target)
  2) <prefix>Biomarker_by_target.csv  rollup:   per target, all trials' biomarkers aggregated

Usage:
  python3 assemble_biomarkers.py --in biomarker_rows.tsv --out-dir data --prefix "PD_toxin_target_"
"""
import argparse, csv, os
from collections import defaultdict, Counter

IN=["target","drug","nct","phase","status","outcome_type","biomarker_category","measure","time_frame","conditions"]
PDKEY=("parkinson","neurodeg","dementia","alzheimer","cognit","dopamin","lewy","tremor","synuclein")
CAT_ORDER=["PK","imaging","fluid_molecular","inflammation","genomic","other_biomarker"]

def load(path):
    recs=[]
    for ln in open(path).read().splitlines()[1:]:
        if not ln.strip(): continue
        c=ln.split("\t")
        if len(c)<10: continue
        recs.append(dict(zip(IN,c[:10])))
    seen=set(); out=[]
    for r in recs:
        k=(r["target"],r["drug"].upper(),r["nct"],r["outcome_type"],r["measure"])
        if k in seen: continue
        seen.add(k); out.append(r)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--in",dest="inp",required=True)
    ap.add_argument("--out-dir",required=True); ap.add_argument("--prefix",default="")
    a=ap.parse_args()
    recs=load(a.inp)
    recs.sort(key=lambda r:(r["target"],
        CAT_ORDER.index(r["biomarker_category"]) if r["biomarker_category"] in CAT_ORDER else 9,
        r["drug"].upper(), r["nct"]))

    g=os.path.join(a.out_dir,f"{a.prefix}Biomarker_measures.csv")
    with open(g,"w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["Target","Drug","NCT_ID","URL","Phase","Status","Outcome_type",
                    "Biomarker_category","Measure","Time_frame","Conditions","PD_neuro"])
        for r in recs:
            pdn="Yes" if any(k in r["conditions"].lower() for k in PDKEY) else "No"
            url=f"https://clinicaltrials.gov/study/{r['nct']}" if r["nct"] else ""
            w.writerow([r["target"],r["drug"],r["nct"],url,r["phase"],r["status"],r["outcome_type"],
                        r["biomarker_category"],r["measure"],r["time_frame"],r["conditions"],pdn])
    print("wrote",g,"rows",len(recs))

    bt=os.path.join(a.out_dir,f"{a.prefix}Biomarker_by_target.csv")
    by=defaultdict(list)
    for r in recs: by[r["target"]].append(r)
    with open(bt,"w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["Target","N_drugs","N_trials","N_biomarker_measures",
                    "PK","imaging","fluid_molecular","inflammation","genomic","other_biomarker",
                    "N_PDneuro_trials","Distinct_drugs","Distinct_biomarker_measures"])
        for tgt in sorted(by):
            rs=by[tgt]; cats=Counter(r["biomarker_category"] for r in rs)
            drugs=sorted(set(r["drug"] for r in rs))
            trials=set(r["nct"] for r in rs)
            pdtrials=set(r["nct"] for r in rs if any(k in r["conditions"].lower() for k in PDKEY))
            dmeas=sorted(set(r["measure"] for r in rs if r["measure"]))
            dm=" | ".join(dmeas)
            if len(dm)>1500: dm=dm[:1500]+f" ... (+{len(dmeas)} distinct total)"
            w.writerow([tgt,len(drugs),len(trials),len(rs),
                        cats.get("PK",0),cats.get("imaging",0),cats.get("fluid_molecular",0),
                        cats.get("inflammation",0),cats.get("genomic",0),cats.get("other_biomarker",0),
                        len(pdtrials),"; ".join(drugs),dm])
    print("wrote",bt,"targets",len(by))

if __name__=="__main__": main()
