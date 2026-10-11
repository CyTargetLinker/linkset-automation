# SIGNOR linkset

Protein → protein causal signalling linkset for human, built from SIGNOR.

| | |
| --- | --- |
| Source | https://signor.uniroma2.it (`getData.php?organism=9606`) |
| License | CC BY 4.0 |
| DOI | not yet deposited, needs the one-time setup below |
| Workflow | `.github/workflows/create-linkset-signor.yml` |
| Species | hsa |

## How it works

`signor.py` downloads SIGNOR's human network, keeps the UniProt → UniProt rows,
and writes `input_signor_protein.txt` plus a config generated into `config/`.
Both ends enter as UniProt (syscode `S`) and BridgeDb expands them to Ensembl,
Entrez and HGNC through the human bundle.

```bash
python3 signor/signor.py
python3 signor/signor.py --version 2026-10      # override the version label only
```

October 2026: 19,636 protein edges, 99.3% of nodes mapped.

## Chemicals

A ChEBI chemical → protein linkset (2,149 edges in October 2026) is written by
the same code but switched off in `LINKSETS`. It needs the BridgeDb metabolite
database, which linkset-creator v2.0 cannot open: its 2015 BridgeDb libraries
fail on connect with `No DataSource known for Cks` and write an empty linkset
([linksetCreator#5](https://github.com/CyTargetLinker/linksetCreator/issues/5)).
`signor_chemical.config.template` is ready for when a newer jar exists. A
separate linkset is needed because linkset-creator fixes the node type, input
syscode and BridgeDb file per side of a config, so proteins and chemicals cannot
share a side.

## Input

SIGNOR records one row per piece of evidence, so a pair often has several rows
with different effects, mechanisms or PMIDs. They are collapsed to one edge per
(source, target), joining the effects, mechanisms, PMIDs and SIGNOR ids. One
edge per pair is not a choice: linkset-creator keeps only the first row it sees
for a pair, and it also drops A → B when B → A came first
([linksetCreator#4](https://github.com/CyTargetLinker/linksetCreator/issues/4)).
The `sign` column is `+`, `-`, `+;-` when the rows disagree, or empty
when no effect is signed.

UniProt isoform and processed-chain suffixes (`P05771-2`, `P01019-PRO_0000032457`)
are stripped, since BridgeDb knows canonical accessions only.

Rows are not filtered on `TAX_ID`. That column is the species the experiment was
done in, not the species of the proteins, and the human file holds thousands of
rows on human proteins tagged `-1` or `10090`.

Left out: chemicals (see above), antibodies, and every complex, protein family,
phenotype, stimulus and microRNA, which carry SIGNOR-internal ids BridgeDb
cannot map.

## Versioning

The version is the current monthly release as `YYYY-MM`, read from the filename
SIGNOR's release endpoint offers (`Oct2026_release.txt` → `2026-10`). The format
sorts correctly as a string, which the Zenodo version gate relies on.

The content always comes from the live database, because the release file mixes
every species and cannot be sliced to human the way `getData` does. `--version`
therefore only changes the label and cannot rebuild an older release.

## Release

Not yet deposited. Dispatch the workflow once with `bootstrap_zenodo` to create
the first draft: it builds, QCs, uploads and sets the metadata, then stops
without publishing, because minting a DOI cannot be undone. Review it at
https://zenodo.org/me/uploads while logged in as the account that owns
`ZENODO_ACCESS_TOKEN`, publish by hand, then set `ZENODO_DEPOSITION_ID` and
`ZENODO_CONCEPT_ID` in the workflow and add a schedule.

Metadata lives in `zenodo-metadata.sh` so the bootstrap and later versions cannot
drift apart.
