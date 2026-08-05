# ChEMBL mechanism-of-action linkset

`CHEMBL_MoA.py` pulls human drug mechanism-of-action records from ChEMBL and writes
`input.txt`; `chembl_human_MoA.config` turns it into a compound → gene linkset.

Smoke test without the full fetch:

```bash
python3 CHEMBL/CHEMBL_MoA.py --limit 40
```

## Gotchas

**Rows must stay full width.** linkset-creator v2.0 reads rows with Java's
`String.split("\t")`, which drops trailing empty fields and ignores quoting. A value
containing a newline, or an empty final column, gives a short row; the jar throws
`ArrayIndexOutOfBoundsException`, catches it internally and **exits 0 having
truncated the linkset there**. On 2026-08-05 that silently cut the output to 757 of
7561 rows. Hence the whitespace collapsing and the constant `row_end` column. The
only symptom in QC is a single edge with no `datasource`.

**The gene id is an HGNC accession, syscode `Hac`.** Column 23 holds `HGNC:1381`, so
`H` (HGNC symbol) and `L` (Entrez Gene) both map nothing. Measured over the full
dataset: `Hac` maps 99.8% and gives each gene node its symbol as label; column 22
with `H` maps 84.7% and leaves every label equal to its own id.

**One row per target component.** A ChEMBL target may be a protein complex: the
GABA-B receptor is GABBR1 plus GABBR2. Joining them into one `GABBR1;GABBR2` cell
creates a gene node that matches nothing, which is what capped mapping at ~85%.
Each component gets its own row instead, so a row is identified by `mechanism_id`
**plus** the gene, not `mechanism_id` alone. Components must be read one at a time —
zipping separate symbol and HGNC lists misaligns them, because a component can have
a symbol but no HGNC id. Keeping only the first component would cost ~6900 edges.
