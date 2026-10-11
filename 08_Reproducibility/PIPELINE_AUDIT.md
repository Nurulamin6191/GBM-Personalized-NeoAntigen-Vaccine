# Pipeline check (re-run 2026-10-11, after consolidation)

Method: static re-verification (raw FASTQs/BAMs retired with Last/; compute env
retired with rebuild/). Every check below runs on checked-in files only.

## Checks performed
- [x] All 8 shell stage scripts parse (`bash -n`).
- [x] All 9 Python scripts parse (`ast.parse`).
- [x] Table reconciliation: missense_dbsnp 4,265 = 2,126 + 1,815 + 324 (script-verified).
- [x] Key values exact: LLFGYPVYV 10.31053 / ALAAAVAAV 16.94121 (Galaxy vs fresh MHCflurry).
- [x] All 11 figure generators run clean; PNGs regenerated (fig01–fig11).
- [x] Doc cross-references audited; stale Galaxy/Last/rebuild paths fixed or manifest-noted.
- [x] `git status` clean; `core.fileMode false` set for the FUSE mount; remote in sync.
- [x] Table 1 rebuilds byte-identically with zero network (rare_genes.tsv path).

## Defects found and fixed in this pass
1. `.git` lost during folder consolidation → recovered via fresh clone + `git reset`
   (history intact; verified log + clean tree).
2. `.gitignore` dropped in the same move → restored with RUNTIME-dir entries.
3. `tesla_audit.py` IEDB column misalignment (13-field rows vs 12-field header) →
   positional parse; verified 10/10 presentation pass.
4. `fetch_rare_genes.py` format-string arity bug → fixed; 1,364 genes mapped.
5. `build_table1.py`/`gate_expression.py` depended on deleted VCF → refactored to
   `rare_genes.tsv`; Table 1 verified byte-identical without network.

## Residual risks (cannot re-execute without raw data/env)
- Stage scripts 01–05 (fetch/QC/align/call/type) are syntax-checked but not re-run
  end-to-end here; their outputs survive as curated logs/tables.
- MHCflurry/NetMHCpan runs survive as CSVs + IEDB API script (re-runnable anywhere).
- git identity is the repo owner (set locally after .git recovery).
