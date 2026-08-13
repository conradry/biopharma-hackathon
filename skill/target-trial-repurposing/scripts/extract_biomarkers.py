#!/usr/bin/env python3
"""
extract_biomarkers.py — pull ONLY the biomarker-type outcome measures (full text,
categorized) from every drug's ClinicalTrials.gov trials, tagged by target.

Companion to extract_target_trials.py (which pulls ALL outcomes). This one filters
ctgov.design_outcomes to biomarker-matching measures and emits one row per
(target × trial × biomarker outcome measure), so results can be aggregated by target.

INPUT  (--drugs): the same drugs.tsv as extract_target_trials.py:
    search_term<TAB>display<TAB>targets      (targets = '; '-joined gene symbols)
OUTPUT (--out): TSV, columns = COLS below (exploded by target).

Usage:
    python3 extract_biomarkers.py --drugs drugs.tsv --out biomarker_rows.tsv
Then aggregate with assemble_biomarkers.py.
Requires: authenticated `paperclip` CLI.
"""
import argparse, subprocess, sys, time
PAD = "¤"
COLS = ["target","drug","nct","phase","status","outcome_type","biomarker_category",
        "measure","time_frame","conditions"]

# SQL-side filter: a measure is a candidate biomarker if it matches any of these.
BT = ("'%biomarker%','%synuclein%','%neurofilament%','%cerebrospinal%','% csf%','%glucocerebrosidase%',"
  "'%amyloid%','% tau %','%pharmacokinet%','% pk %','%plasma concentration%','%serum concentration%',"
  "'%cmax%','%auc%','%trough%','%receptor occupancy%','%pet imaging%','%pet scan%','%fdg-pet%','%fdg pet%',"
  "'%[123i]%','%[18f]%','%mibg%','%datscan%','%dat-spect%','%dopamine transporter%','%scintigraph%',"
  "'%cytokine%','%c-reactive%','% crp %','%interleukin%','%mutation%','%pd-l1%','%pd-1%','%gene expression%',"
  "'%circulating tumor%','%ctdna%','%egfr%'")

# Python-side categorization (priority order; specific analyte beats generic PK).
_CATS=[("imaging",("pet imaging","pet scan","pet/ct"," pet "," pet)","fdg-pet","fdg pet","spect imaging",
          "spect scan"," spect "," spect)","datscan","dat-spect","mibg","[18f]","[123i]","scintigraph",
          "regional glucose","glucose metabolism","dopamine transporter imaging")),
       ("fluid_molecular",("synuclein","neurofilament"," nfl ","cerebrospinal"," csf ","glucocerebrosidase",
          "amyloid","tau protein","p-tau","total tau","phospho-tau","bdnf")),
       ("inflammation",("cytokine","c-reactive"," crp ","interleukin","tnf-alpha","tnf-α","tnf level","il-6")),
       ("genomic",("mutation","pd-l1","pd-1 ","gene expression","circulating tumor","ctdna","egfr","genotype",
          "dna methylation","rna expression","microrna")),
       ("PK",("pharmacokinet"," pk ","plasma concentration","serum concentration","cmax","c max","auc",
          "area under the concentration","trough","steady state","receptor occupancy","half-life","clearance"))]
def categorize(m):
    s=" "+m.lower()+" "
    for cat,kws in _CATS:
        if any(k in s for k in kws): return cat
    return "other_biomarker"

def rp(col,start,n): return f"RPAD(COALESCE(SUBSTRING({col},{start},{n}),''),{n},'{PAD}')"

def build_sql(term, offset, limit=60):
    esc=term.replace("'","''")
    mc=", ".join(f"{rp('o.measure',1+55*k,55)} AS m{k+1}" for k in range(4))
    cc=", ".join(f"{rp('cond.c',1+55*k,55)} AS c{k+1}" for k in range(2))
    return f"""
WITH t AS (SELECT DISTINCT nct_id FROM ctgov.interventions WHERE name ILIKE '%{esc}%')
SELECT o.nct_id, s.phase, s.overall_status, o.outcome_type,
  {mc}, {rp('o.time_frame',1,40)} AS tf, {cc}
FROM t
JOIN ctgov.design_outcomes o ON o.nct_id=t.nct_id
JOIN ctgov.studies s ON s.nct_id=t.nct_id
LEFT JOIN LATERAL (SELECT translate(STRING_AGG(DISTINCT c.name,'; '),E'\\n\\r|\\t','    ') AS c
   FROM ctgov.conditions c WHERE c.nct_id=t.nct_id) cond ON true
WHERE o.measure ILIKE ANY(ARRAY[{BT}])
ORDER BY o.nct_id, o.outcome_type, o.id LIMIT {limit} OFFSET {offset}"""

def parse(text):
    return [ln.split(" | ") for ln in text.splitlines()
            if ln.startswith("NCT") and len(ln.split(" | "))==11]

def clean(x):
    x=x.replace(PAD,""); return "" if x=="NULL" else x

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--drugs",required=True); ap.add_argument("--out",required=True)
    ap.add_argument("--page",type=int,default=60)
    a=ap.parse_args()
    drugs=[]
    for ln in open(a.drugs).read().splitlines()[1:]:
        if not ln.strip(): continue
        c=ln.split("\t")
        drugs.append((c[0].strip(), c[1].strip() if len(c)>1 else c[0].strip(),
                      c[2].strip() if len(c)>2 else ""))
    open(a.out,"w").write("\t".join(COLS)+"\n")
    total=0
    for term,disp,targets in drugs:
        genes=[g.strip() for g in targets.split(";") if g.strip()] or [""]
        off=0; dc=0
        while True:
            sql=build_sql(term,off,a.page)
            try: r=subprocess.run(["paperclip","sql","-s","trials",sql],capture_output=True,text=True,timeout=90)
            except subprocess.TimeoutExpired:
                try: r=subprocess.run(["paperclip","sql","-s","trials",sql],capture_output=True,text=True,timeout=90)
                except subprocess.TimeoutExpired: print(f"[skip] {disp} off={off}",file=sys.stderr); break
            rows=parse(r.stdout+r.stderr)
            if not rows: break
            with open(a.out,"a") as f:
                for p in rows:
                    nct=p[0].strip(); phase=clean(p[1].strip()); status=clean(p[2].strip())
                    otype=p[3].strip()
                    measure=clean((p[4]+p[5]+p[6]+p[7]).rstrip(PAD))
                    tf=clean(p[8].rstrip(PAD).strip()); conds=clean((p[9]+p[10]).rstrip(PAD))
                    cat=categorize(measure)
                    for g in genes:
                        f.write("\t".join(x.replace("\t"," ") for x in
                            [g,disp,nct,phase,status,otype,cat,measure,tf,conds])+"\n")
                    total+=1
            dc+=len(rows); off+=a.page
            if len(rows)<a.page: break
        print(f"[{time.strftime('%H:%M:%S')}] {disp} [{','.join(genes)}]: {dc}  [total {total}]",flush=True)
    print(f"DONE total={total}")

if __name__=="__main__": main()
