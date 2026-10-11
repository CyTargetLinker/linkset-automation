# AOP-Wiki linksets

Two linksets from the AOP-Wiki, built from AOP-Wiki RDF.

| | |
| --- | --- |
| Source | https://aopwiki.org, via the AOP-Wiki RDF SPARQL endpoint (https://aopwiki.rdf.bigcat-bioinformatics.org/sparql) |
| License | CC BY-SA 4.0, the license on AOP-Wiki content |
| DOI | not yet deposited, needs the one-time setup below |
| Workflow | `.github/workflows/create-linkset-aopwiki.yml` |

| Linkset | Edges (2026-10-10) | Source nodes | Target nodes |
| --- | --- | --- | --- |
| `aopwiki_ke_gene.xgmml` | 3,708 | 455 key events | 1,687 genes, 99.9% mapped |
| `aopwiki_stressor_mie.xgmml` | 390 | 267 chemicals | 104 molecular initiating events |

Both use the same key event ids (`aop.events:<n>`), so a network extended with
both connects chemicals to genes through the key events they share.

```bash
python3 aopwiki/aopwiki.py
```

## Key event → gene

The AOP-Wiki has no gene table. These pairs are the HGNC genes that the AOP-Wiki
RDF build text-mines from each key event's title and description
(`edam:data_1025`), so they are **not curated**. The `detected_by` edge attribute
says whether named-entity recognition, the regular expression, or both found
the gene. Genes enter as HGNC accessions (`HGNC:108`, syscode `Hac`) and BridgeDb
expands them to Ensembl, Entrez, HGNC symbol and UniProt.

## Stressor → molecular initiating event

The AOP-Wiki attaches stressors to whole AOPs, not to key events. A prototypical
stressor is what starts the pathway, so each chemical stressor is linked to the
molecular initiating events of every AOP that lists it. The `aops` edge
attribute records those AOPs.

Chemicals are identified by InChIKey only, with CAS number and DSSTox id as node
attributes. They are not expanded through BridgeDb, because linkset-creator
cannot open the current BridgeDb metabolite database
([linksetCreator#5](https://github.com/CyTargetLinker/linksetCreator/issues/5)).
Stressors without a chemical entity (physical stressors, mixtures, unspecified
ones) are left out.

## Licensing

AOP-Wiki content is CC BY-SA 4.0, so the deposit is share-alike as well. Only
AOPs that carry that license in the RDF are used (529 of 600 in 2026-10-10). A key
event is included if it sits on at least one of them.

## Versioning

The version is the AOP-Wiki RDF release the endpoint serves (`pav:version`
`2026.10.10` becomes `2026-10-10`), which sorts correctly as a string for the
Zenodo version gate. The content is always what the endpoint serves now, so
`--version` only changes the label.

The endpoint truncates any result at 10,000 rows without saying so. Every
large query is paged with `ORDER BY`/`LIMIT`/`OFFSET` for that reason; do not add
one that isn't.

## Release

Not yet deposited. Dispatch the workflow once with `bootstrap_zenodo` to create
the first draft: it builds, QCs, uploads and sets the metadata, then stops
without publishing, because minting a DOI cannot be undone. Review it at
https://zenodo.org/me/uploads while logged in as the account that owns
`ZENODO_ACCESS_TOKEN`, publish by hand, then set `ZENODO_DEPOSITION_ID` and
`ZENODO_CONCEPT_ID` in the workflow and add a quarterly schedule.

Metadata lives in `zenodo-metadata.sh` so the bootstrap and later versions cannot
drift apart.
