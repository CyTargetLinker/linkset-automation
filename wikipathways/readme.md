# WikiPathways linkset

Pathway–gene linkset for human, mouse and rat, built from the WikiPathways GMT
release.

| | |
| --- | --- |
| Source | https://www.wikipathways.org (GMT files from https://data.wikipathways.org) |
| License | CC BY 4.0 as deposited; WikiPathways content itself is CC0 1.0 |
| DOI | [10.5281/zenodo.4500957](https://doi.org/10.5281/zenodo.4500957) (concept DOI, resolves to the latest version) |
| Workflow | `.github/workflows/create-linkset-wikipathways.yml` |

## How it works

`wp.py` downloads `wikipathways-<VERSION>-gmt-<Species>.gmt` for each species in
`SPECIES` and writes one `input_<code>.txt` per species. The GMT set name packs
four fields separated by `%`: `pathway name%WikiPathways_<date>%WP####%species`.

`wp_hsa.config`, `wp_mmu.config` and `wp_rno.config` drive the LinksetCreator.
Genes enter as Entrez (`target_syscode_in=L`) and are expanded through BridgeDb
to Ensembl, HGNC, Entrez and UniProt (`target_syscodes_out=En,H,L,S`).

This is the only resource with a full pipeline: build, QC, and deposit to
Zenodo.

## Updating to a new release

1. Bump `VERSION` in `wp.py`.
2. Bump `version=` in all three `wp_*.config` files to match.
3. Bump `version=` in the Zenodo metadata step of the workflow.

## Known limitation: gene labels are Entrez IDs, not symbols

The WikiPathways GMT contains only Entrez Gene IDs — there are no symbols in it.
`wp.py` therefore writes the same value into both the `gene_id` and
`gene_symbol` columns, and `target_label_column=4` points at that duplicate, so
gene nodes are labelled with numbers rather than symbols in Cytoscape.

For contrast, the GO linkset avoids this by joining symbols in from NCBI
`gene_info`, and an earlier SPARQL-based WikiPathways run produced real symbols
because the query returned a `GeneName` column. Fixing this needs either a
symbol lookup after the GMT parse or a return to the SPARQL endpoint for the
input file.
