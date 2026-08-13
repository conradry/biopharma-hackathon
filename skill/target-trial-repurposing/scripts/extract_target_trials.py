#!/usr/bin/env python3
"""
extract_target_trials.py — pull ALL ClinicalTrials.gov trials + outcome/biomarker
measures for a list of drugs, via the Paperclip trials SQL (AACT schema).

Handles the two Paperclip quirks that make naive extraction lossy:
  1. The SQL ASCII renderer hard-caps each cell at ~60 chars  -> long outcome
     text is chunked into 55-char SUBSTRING columns, sentinel-padded ('¤'), and
     reassembled here so no characters and no boundary spaces are lost.
  2. paperclip trims its own stdout at ~48KB  -> we page at 60 rows/query.
Also: free-text measures may contain newlines/pipes -> translate() strips them in SQL.

INPUT  (--drugs): TSV with a header and columns:  search_term<TAB>display<TAB>targets
   search_term = substring matched against ctgov.interventions.name (ILIKE '%term%')
   display     = pretty drug name for the output
   targets     = '; '-joined gene targets this drug hits (free text, carried through)
OUTPUT (--out): TSV, one row per (drug, trial), columns = COLS below.

Usage:
   python3 extract_target_trials.py --drugs drugs.tsv --out trials_rows.tsv
Requires: `paperclip` CLI on PATH, authenticated (`paperclip config`).
"""
import argparse, subprocess, sys, time, os

PAD = "¤"
COLS = ["targets","drug","nct","phase","status","start_date","enrollment","conditions",
        "n_primary","n_secondary","primary_outcomes","secondary_outcomes","biomarker"]

# Keyword net for the SQL-side biomarker flag (bool_or over ALL outcomes, untruncated).
# Broaden/narrow to taste; keep patterns specific to avoid false hits (e.g. NOT '% spect%').
BIOMARKER_TERMS = ("'%biomarker%','%synuclein%','%neurofilament%','%cerebrospinal%',"
    "'%glucocerebrosidase%','%amyloid%','%pharmacokinet%','% pk %','%plasma concentration%',"
    "'%serum concentration%','%pet imaging%','%pet scan%','%fdg-pet%','%[123i]%','%[18f]%',"
    "'%mibg%','%cytokine%','%egfr mutation%','%pd-l1%','%receptor occupancy%','%dat-spect%',"
    "'%datscan%','%dopamine transporter%','%gene expression%','%circulating tumor%'")

def rpad(expr):
    return f"RPAD(COALESCE({expr},''),55,'{PAD}')"

def build_sql(term, offset, limit=60):
    esc = term.replace("'", "''")
    cc = ", ".join(f"{rpad(f'SUBSTRING(cond.c,{1+55*k},55)')} AS c{k+1}" for k in range(2))
    pc = ", ".join(f"{rpad(f'SUBSTRING(od.pm,{1+55*k},55)')} AS p{k+1}" for k in range(4))
    sc = ", ".join(f"{rpad(f'SUBSTRING(od.sm,{1+55*k},55)')} AS s{k+1}" for k in range(4))
    return f"""
WITH t AS (SELECT DISTINCT s.nct_id, s.phase, s.overall_status, s.start_date, s.enrollment
  FROM ctgov.studies s JOIN ctgov.interventions i ON i.nct_id=s.nct_id AND i.name ILIKE '%{esc}%'
  ORDER BY s.nct_id LIMIT {limit} OFFSET {offset})
SELECT t.nct_id, t.phase, t.overall_status, t.start_date, t.enrollment,
  {cc}, od.np, od.ns, {pc}, {sc}, od.bmflag AS biomarker
FROM t
LEFT JOIN LATERAL (SELECT translate(STRING_AGG(DISTINCT c.name,'; '), E'\\n\\r|\\t','    ') AS c
  FROM ctgov.conditions c WHERE c.nct_id=t.nct_id) cond ON true
LEFT JOIN LATERAL (SELECT COUNT(*) FILTER (WHERE outcome_type='primary') AS np,
  COUNT(*) FILTER (WHERE outcome_type='secondary') AS ns,
  translate(STRING_AGG(DISTINCT measure,' ~ ') FILTER (WHERE outcome_type='primary'), E'\\n\\r|\\t','    ') AS pm,
  translate(STRING_AGG(DISTINCT measure,' ~ ') FILTER (WHERE outcome_type='secondary'), E'\\n\\r|\\t','    ') AS sm,
  bool_or(measure ILIKE ANY(ARRAY[{BIOMARKER_TERMS}])) AS bmflag
  FROM ctgov.design_outcomes o WHERE o.nct_id=t.nct_id) od ON true
ORDER BY t.nct_id"""

def parse_page(text):
    rows = []
    for ln in text.splitlines():
        if not ln.startswith("NCT"):
            continue
        parts = ln.split(" | ")
        if len(parts) != 18:
            continue
        rows.append(parts)
    return rows

def reassemble(p):
    def clean(x):
        x = x.replace(PAD, "")
        return "" if x == "NULL" else x
    conditions = (p[5]+p[6]).rstrip(PAD)
    primary    = (p[9]+p[10]+p[11]+p[12]).rstrip(PAD)
    secondary  = (p[13]+p[14]+p[15]+p[16]).rstrip(PAD)
    bm = p[17].strip()
    return dict(nct=p[0].strip(), phase=clean(p[1].strip()), status=clean(p[2].strip()),
        start_date=clean(p[3].strip()), enrollment=clean(p[4].strip()),
        conditions=clean(conditions), n_primary=clean(p[7].strip()), n_secondary=clean(p[8].strip()),
        primary_outcomes=clean(primary), secondary_outcomes=clean(secondary),
        biomarker=("Yes" if bm=="True" else ("No" if bm=="False" else "")))

def run_sql(sql, timeout=90):
    return subprocess.run(["paperclip","sql","-s","trials",sql],
        capture_output=True, text=True, timeout=timeout)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drugs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--page", type=int, default=60)
    a = ap.parse_args()

    drugs = []
    for ln in open(a.drugs).read().splitlines()[1:]:
        if not ln.strip(): continue
        c = ln.split("\t")
        drugs.append((c[0].strip(), c[1].strip() if len(c)>1 else c[0].strip(),
                      c[2].strip() if len(c)>2 else ""))

    with open(a.out, "w") as f: f.write("\t".join(COLS)+"\n")
    total = 0
    for term, disp, targets in drugs:
        off = 0; dc = 0
        while True:
            sql = build_sql(term, off, a.page)
            try: r = run_sql(sql)
            except subprocess.TimeoutExpired:
                try: r = run_sql(sql)
                except subprocess.TimeoutExpired:
                    print(f"[skip] {disp} offset={off} timeout", file=sys.stderr); break
            rows = parse_page(r.stdout + r.stderr)
            if not rows: break
            with open(a.out, "a") as f:
                for p in rows:
                    d = reassemble(p)
                    f.write("\t".join(x.replace("\t"," ") for x in
                        [targets, disp, d["nct"], d["phase"], d["status"], d["start_date"],
                         d["enrollment"], d["conditions"], d["n_primary"], d["n_secondary"],
                         d["primary_outcomes"], d["secondary_outcomes"], d["biomarker"]])+"\n")
            dc += len(rows); total += len(rows); off += a.page
            if len(rows) < a.page: break
        print(f"[{time.strftime('%H:%M:%S')}] {disp} ({targets}): {dc}  [total {total}]", flush=True)
    print(f"DONE total={total}")

if __name__ == "__main__":
    main()
