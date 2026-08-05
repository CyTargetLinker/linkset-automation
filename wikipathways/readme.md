# WikiPathways linkset

Pathway–gene linksets for 13 species, built from the WikiPathways GMT release.

| | |
| --- | --- |
| Source | https://www.wikipathways.org (GMT files from https://data.wikipathways.org) |
| License | CC BY 4.0 as deposited; WikiPathways content itself is CC0 1.0 |
| DOI | [10.5281/zenodo.4500957](https://doi.org/10.5281/zenodo.4500957) (concept DOI, resolves to the latest version) |
| Workflow | `.github/workflows/create-linkset-wikipathways.yml` |
| Schedule | quarterly — 03:17 UTC on the 15th of Jan / Apr / Jul / Oct |
| Species | hsa mmu rno bta sce cel dre ptr cfa gga dme ath aga |

This is the only resource with a full pipeline: build, QC, and deposit to
Zenodo.

## How it works

`wp.py` resolves the release, downloads one GMT per species, and writes
`input_<code>.txt` plus a generated `config/wp_<code>.config` for each. Genes
enter as Entrez (`target_syscode_in=L` for every species) and BridgeDb expands
them to Ensembl, Entrez and UniProt, plus HGNC for human.

```bash
python3 wikipathways/wp.py                        # all species, newest release
python3 wikipathways/wp.py --species sce aga      # just these two
python3 wikipathways/wp.py --version 20260510     # a specific release
```

It also writes `build/manifest.tsv`, listing the species **actually built**.
The workflow iterates that rather than `species.tsv`, so a skipped species can
never leave a later step looking for a file that was never produced.

## Updating to a new release

**Nothing to do.** The release is resolved on every run and flows into the GMT
URLs, every generated config, and the Zenodo record version.

`data.wikipathways.org` keeps only a rolling 12-month window, so a pinned
release would eventually 404. `wp.py` walks the dated release directories
newest first and accepts one only once it holds a GMT for every `core` species.
`/current/` is not used: it is an alias that says nothing about whether the
release behind it is complete. An incomplete newest release is skipped for the
previous one, which is safe because publishing refuses to go backwards.

## Adding or removing a species

One line in [`species.tsv`](species.tsv); nothing else changes. The columns are
documented in its header.

Five species are listed but commented `off` — pig, horse, poplar, tomato and
maize have between one and five pathways each. `wp.py` also skips any `extra`
species below `MIN_PATHWAYS` (5).

## Generated configs

There are no committed per-species configs. Each `config/wp_<code>.config` is
generated from [`wp.config.template`](wp.config.template) on every run, and
`config/` is gitignored. Do not hand-edit them.

## Quality control

`scripts/qc_xgmml.py` checks structure and **BridgeDb mapping coverage**:
`--mapped-type gene --min-mapped-frac` fails a linkset unless that fraction of
gene nodes gained more than one identifier. A mapper that failed to load leaves
every node carrying only its input id — valid XML, useless for matching.

Coverage on release 20260710 was 100% for most species, 98.9% cow, 95.9% rat.
Floors are 0.90 (`core`) and 0.80 (`extra`).

Any QC failure fails the whole run and publishes nothing, rather than minting a
permanent DOI on an incomplete record.

## Publishing

Zenodo versions are immutable, so the run will not publish if the release is
already published, if it is *older* than the published one, or if a file from
the previous version is missing from this build. Needs the
`ZENODO_ACCESS_TOKEN` secret; `skip_zenodo` and `force_publish` make the gates
testable.

## Note on BridgeDb downloads

Use `scripts/fetch-bridge.sh`. Figshare redirects to a presigned S3 URL that
expires in ten seconds, and answers browser-like user agents with `HTTP 202
Accepted` and an empty body. 202 is a success code, so `wget` exits 0 and writes
a zero-byte file — which is how the previous download failed. The script sends
no spoofed user-agent and accepts only an explicit 200 that passes a zip test.

## Known limitation: gene labels are Entrez IDs, not symbols

The GMT contains only Entrez Gene IDs, so `wp.py` writes the same value into
both the `gene_id` and `gene_symbol` columns and gene nodes are labelled with
numbers in Cytoscape. Fixing it needs a symbol lookup after the GMT parse — NCBI
`gene_info` maps GeneID to Symbol — or a return to the SPARQL endpoint, whose
query returned a `GeneName` column.
