# ChEMBL mechanism-of-action linkset

Compound → gene linkset for human, built from ChEMBL drug mechanism-of-action
records.

| | |
| --- | --- |
| Source | https://www.ebi.ac.uk/chembl/ (via `chembl_webresource_client`) |
| License | CC BY-SA 3.0 |
| DOI | not yet deposited |
| Workflow | `.github/workflows/Chembl_MoA.yml` |
| Species | human |

## How it works

`CHEMBL_MoA.py` pulls human mechanism records and writes `input.txt`,
`build/version.txt` and a generated `config/chembl_human_MoA.config`. Genes enter
as HGNC accessions (`target_syscode_in=Hac`), are labelled with their gene
symbol, and BridgeDb expands them to Entrez, Ensembl, HGNC and UniProt.

```bash
python3 CHEMBL/CHEMBL_MoA.py                 # full build
python3 CHEMBL/CHEMBL_MoA.py --limit 40      # smoke test, minutes not an hour
```

A full run takes about an hour in CI. The client caches responses, so a local
re-run is much faster than the first one.

## Targets

A ChEMBL target may be a protein complex — the GABA-B receptor is GABBR1 plus
GABBR2 — so each target component becomes its own compound → gene row. A row is
therefore identified by `mechanism_id` **plus** the gene, not `mechanism_id`
alone.

Components are read one at a time rather than by zipping separate symbol and
HGNC lists, which would misalign whenever a component has a symbol but no HGNC
id.

## Versioning

The release is read from ChEMBL's `status.json` at run time and flows into the
config, the graph version and the Zenodo metadata, so nothing needs editing when
ChEMBL releases. The config is generated from `chembl_human_MoA.config.template`
and is gitignored — do not hand-edit it.

## Release

Not yet deposited, and the deposit license is still undecided.
