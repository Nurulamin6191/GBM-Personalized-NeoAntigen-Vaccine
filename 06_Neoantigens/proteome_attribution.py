#!/usr/bin/env python3
"""Proteome attribution: exact 9mer search in reviewed human proteome (UniProt UP000005640).
Regenerates 06_Neoantigens/proteome_attribution.tsv from refs-RUNTIME proteome FASTA.
Fetch: https://ftp.uniprot.org/.../UP000005640/UP000005640_9606.fasta.gz
"""
import csv
import gzip
import sys

FASTA = sys.argv[1] if len(sys.argv) > 1 else "refs-RUNTIME/proteome/UP000005640_9606.fasta.gz"
OUT = "06_Neoantigens/proteome_attribution.tsv"
PEPTIDES = ["ALAAAVAAV", "LLFGYPVYV", "SAKEMQRWV", "SAKEMQRWI", "EMQRWVMEK", "SYFPEITHI"]


def main():
    opener = gzip.open if FASTA.endswith(".gz") else open
    seqs, hid = {}, None
    with opener(FASTA, "rt") as f:
        for line in f:
            if line.startswith(">"):
                hid = line[1:].split()[0]
                seqs[hid] = []
            else:
                seqs[hid].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    print(f"proteins: {len(seqs)}", flush=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "human_proteome_hits", "example_accessions"])
        for pep in PEPTIDES:
            hits = [h for h, s in seqs.items() if pep in s]
            print(f"{pep}: {len(hits)} {hits[:3]}", flush=True)
            w.writerow([pep, len(hits), ";".join(hits[:10])])
    print(f"ATTRIBUTION_DONE -> {OUT}")


if __name__ == "__main__":
    main()
