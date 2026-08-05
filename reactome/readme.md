# Reactome linkset

Pathway–gene linkset for human, built from Reactome's `NCBI2Reactome.txt`.

| | |
| --- | --- |
| Source | https://reactome.org (`NCBI2Reactome.txt`); gene symbols from NCBI Gene `gene_info` (public domain) |
| License | CC0 1.0 — Reactome dedicates its data to the public domain |
| DOI | not yet deposited — needs the one-time setup below |
| Workflow | `.github/workflows/create-linkset-reactome.yml` |
| Species | hsa |

## How it works

`reactome.py` resolves the release, downloads the pathway–gene table, and writes
`input_<code>.txt` plus a generated `config/reactome_<code>.config` for each
species. Genes enter as Entrez (`target_syscode_in=L`), are labelled with their
NCBI gene symbol, and BridgeDb expands them to Ensembl, Entrez and UniProt, plus
HGNC for human.

```bash
python3 reactome/reactome.py                    # all enabled species, current release
python3 reactome/reactome.py --species hsa      # just this one
python3 reactome/reactome.py --version 97       # a specific release
```

Release 097 human: 2310 pathways, 11337 genes, 48811 edges, 99.8% mapped.

## Input

`NCBI2Reactome.txt` rather than `ReactomePathways.gmt`. The GMT covers human
only, is served as a ZIP archive under a `.gmt` name, and is not available in the
versioned archive, so it cannot support a pinned build. The table carries Entrez
ids for every species and is a plain TSV in both current and archived releases.

Reactome lists a gene twice for a pathway when it is supported both by curation
and by orthology projection, so pairs are de-duplicated and their `TAS`/`IEA`
codes aggregated into one `evidence` column.

Genes NCBI has not named are dropped rather than kept under a bare identifier.
These are pathogen-side genes that Reactome files under the host species in its
infectious-disease pathways, and none of them map in BridgeDb. The share dropped
is gated by `max_reject_frac` in `species.tsv` (human is 3.4%).

## Species

Human only. The other eleven buildable species are listed and commented out in
`species.tsv`, because Reactome projects most non-human content from human by
orthology — all 1467 mouse pathways share their name and numeric stem with a
human pathway. Fly, chicken and yeast carry some natively curated pathways, so
they are worth judging separately rather than as one group.

Three further species have no BridgeDb bundle, and *M. tuberculosis* has no NCBI
per-species `gene_info`, the same blocker recorded for horse, tomato and poplar
in `wikipathways/species.tsv`.

## Versioning

The release is read from Reactome's ContentService at run time and zero-padded to
three digits (`097`), so nothing needs editing when Reactome releases. The
padding keeps the Zenodo version gate's string comparison correct past release
100. Reactome releases quarterly.

## Release

Not yet deposited. Dispatch the workflow once with `bootstrap_zenodo` to create
the first draft: it builds, QCs, uploads and sets the metadata, then stops
without publishing, because minting a DOI cannot be undone. Review it at
https://zenodo.org/me/uploads while logged in as the account that owns
`ZENODO_ACCESS_TOKEN`, publish by hand, then set `ZENODO_DEPOSITION_ID` and
`ZENODO_CONCEPT_ID` in the workflow and add a quarterly schedule.

Metadata lives in `zenodo-metadata.sh` so the bootstrap and later versions cannot
drift apart.
