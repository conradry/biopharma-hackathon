#!/usr/bin/env python3
"""
derive_drugs.py — gene targets -> human UniProt -> approved drugs (ChEMBL) -> drugs.tsv

Produces the drugs.tsv that extract_target_trials.py consumes:
    search_term<TAB>display<TAB>targets

Two derivation modes:
  --mode direct     ChEMBL drugs_by_accession, max_phase=4 (approved, ANNOTATED DIRECT
                    mechanism of action). Best when targets are well-drugged.
  --mode bioactivity  ChEMBL bioactivities_by_accession, max_phase=4 approved compounds
                    with a measured pchembl_value against the target (repurposing net;
                    captures off-target/secondary activity). Slower, noisier.

Usage:
  python3 derive_drugs.py --genes SNCA,LRRK2,MAOB --mode direct --out drugs.tsv
  python3 derive_drugs.py --genes-file genes.txt   --mode bioactivity --out drugs.tsv

Notes:
  * salt forms are collapsed to the base ingredient (split_part(name,' ',1)) so the
    ILIKE trial search matches all salts of a drug.
  * paperclip chembl/proteins SQL previews large results; we page at LIMIT 25.
  * Review drugs.tsv before extraction: drop pathological generic search terms
    (e.g. 'SYNTHETIC', 'ESTROGENS') that would over-match ctgov.interventions.
"""
import argparse, subprocess, sys

def sql(source, q, timeout=120):
    r = subprocess.run(["paperclip","sql","-s",source,q], capture_output=True, text=True, timeout=timeout)
    return r.stdout + r.stderr

def parse_table(text, ncol):
    out=[]
    for ln in text.splitlines():
        if " | " not in ln: continue
        p=[x.strip() for x in ln.split(" | ")]
        if len(p)!=ncol: continue
        if p[0] in ("accession","gene_name","base","drug_name") or p[0].startswith("---"): continue
        out.append(p)
    return out

def map_genes(genes):
    inlist=",".join("'%s'"%g.replace("'","''") for g in genes)
    txt=sql("proteins", f"SELECT accession, gene_name FROM uniprot_v.proteins "
        f"WHERE organism='Homo sapiens' AND gene_name IN ({inlist}) ORDER BY gene_name")
    acc={}
    for p in parse_table(txt,2): acc[p[0]]=p[1]   # accession -> gene
    return acc

def derive(acc, mode):
    inlist=",".join("'%s'"%a for a in acc)
    base_genes={}
    off=0
    if mode=="direct":
        q=("SELECT split_part(d.drug_name,' ',1) AS base, p.gene_name "
           "FROM chembl_v.drugs_by_accession d JOIN uniprot_v.proteins p ON p.accession=d.accession "
           f"WHERE d.accession IN ({inlist}) AND d.max_phase=4 "
           "GROUP BY base, p.gene_name ORDER BY base, p.gene_name")
    else:
        q=("SELECT split_part(b.compound_name,' ',1) AS base, p.gene_name "
           "FROM chembl_v.bioactivities_by_accession b JOIN uniprot_v.proteins p ON p.accession=b.accession "
           f"WHERE b.accession IN ({inlist}) AND b.max_phase=4 AND b.compound_name IS NOT NULL "
           "AND b.pchembl_value IS NOT NULL GROUP BY base, p.gene_name ORDER BY base, p.gene_name")
    while True:
        txt=sql("chembl", q+f" LIMIT 25 OFFSET {off}")
        rows=parse_table(txt,2)
        if not rows: break
        for base,gene in rows:
            if base and not base.startswith("---"):
                base_genes.setdefault(base,set()).add(gene)
        off+=25
        if len(rows)<25: break
    return base_genes

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--genes"); ap.add_argument("--genes-file")
    ap.add_argument("--mode", choices=["direct","bioactivity"], default="direct")
    ap.add_argument("--out", required=True)
    a=ap.parse_args()
    if a.genes: genes=[g.strip() for g in a.genes.split(",") if g.strip()]
    elif a.genes_file: genes=[g.strip() for g in open(a.genes_file).read().replace(",","\n").split() if g.strip()]
    else: sys.exit("provide --genes or --genes-file")

    acc=map_genes(genes)
    missing=[g for g in genes if g not in acc.values()]
    if missing: print(f"[warn] no human UniProt match: {missing}", file=sys.stderr)
    print(f"mapped {len(acc)} targets", file=sys.stderr)

    base_genes=derive(acc, a.mode)
    with open(a.out,"w") as f:
        f.write("search_term\tdisplay\ttargets\n")
        for base in sorted(base_genes):
            disp=base.title() if base.isupper() else base
            f.write(f"{base}\t{disp}\t{'; '.join(sorted(base_genes[base]))}\n")
    print(f"wrote {len(base_genes)} drugs -> {a.out}")
    print("REVIEW drugs.tsv and delete any over-generic search_term before extraction.")

if __name__=="__main__": main()
