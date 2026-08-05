# linkset-automation

[![WikiPathways linkset](https://img.shields.io/github/actions/workflow/status/CyTargetLinker/linkset-automation/create-linkset-wikipathways.yml?label=WikiPathways%20linkset)](https://github.com/CyTargetLinker/linkset-automation/actions/workflows/create-linkset-wikipathways.yml)
[![GO linkset](https://img.shields.io/github/actions/workflow/status/CyTargetLinker/linkset-automation/create-linkset-go.yml?label=GO%20linkset)](https://github.com/CyTargetLinker/linkset-automation/actions/workflows/create-linkset-go.yml)
[![TFLink linkset](https://img.shields.io/github/actions/workflow/status/CyTargetLinker/linkset-automation/create-linkset-tflink.yml?label=TFLink%20linkset)](https://github.com/CyTargetLinker/linkset-automation/actions/workflows/create-linkset-tflink.yml)
[![QC XGMML](https://img.shields.io/github/actions/workflow/status/CyTargetLinker/linkset-automation/qc-xgmml.yml?label=QC%20XGMML)](https://github.com/CyTargetLinker/linkset-automation/actions/workflows/qc-xgmml.yml)
[![BridgeDb cache](https://img.shields.io/github/actions/workflow/status/CyTargetLinker/linkset-automation/update-bridgedb.yml?label=BridgeDb%20cache)](https://github.com/CyTargetLinker/linkset-automation/actions/workflows/update-bridgedb.yml)

Under construction: This repo automatically generates CyTargetLinker linksets from different resources, starting with WikiPathways.

<img width="692" height="347" alt="image" src="https://github.com/user-attachments/assets/100e64c0-d6e5-47c9-8686-7be01c39c257" />

## Linkset resources

**This list is still being extended.** The resources below are the ones
covered so far; more are being added as they are automated, and the table is
updated once a resource is built and deposited.

| Resource | Details | Persistent DOI | License |
| --- | --- | --- | --- |
| [WikiPathways](https://www.wikipathways.org) (13 species) | [`wikipathways/`](wikipathways/readme.md) | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4500957.svg)](https://doi.org/10.5281/zenodo.4500957) | CC BY 4.0 |
| [Gene Ontology](https://geneontology.org) (BP / MF / CC, human) | [`go/`](go/readme.md) | not yet deposited | CC BY 4.0 |
| [TFLink](https://tflink.net) | [`tflink/`](tflink/readme.md) | not yet deposited | Free for non-commercial use |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) (mechanism of action) | [`CHEMBL/`](CHEMBL/readme.md) | not yet deposited | CC BY-SA 3.0 |

Each resource name links to its upstream source, and the Details column links to
that resource's own readme — how it is built, how to update it to a new release,
and what is still outstanding. Resource-specific documentation lives there
rather than in this file.

The DOI column holds the Zenodo *concept* DOI, which always resolves to the
latest published version; individual versions have their own DOIs. For resources
that are not yet deposited the license shown is that of the upstream source data,
and the license of the eventual deposit still needs to be decided.

## (Intended) Workflow
1. Run the bridgeDb update script to update bridgeDb versions for all config files.
2. Check the download link for the resource you want to update (check if an update is required) and copy it into the (data preprocessing) script
3. Run the preprocessing script to get the input.csv file for this version. The config file is available and should not change between versions.
4. Run the LinksetCreator
5. Do the quality control procedure with the obtained xgmml file (see **QC** below).
6. If positive, run the linkset upload-to-storage script
7. Update the CTL website to add the new linkset download link and description

The WikiPathways pipeline no longer needs steps 1-3: it resolves the release and
generates its configs on every run — see [`wikipathways/`](wikipathways/readme.md).

## QC
`scripts/qc_xgmml.py` validates that an XGMML linkset is well-formed and has the structure every CyTargetLinker linkset carries (a `<graph>`, nodes with an `identifiers` list + `type`, and edges whose `source`/`target` reference existing nodes plus `datasource`/`interaction`). Empty graphs fail by default.

```bash
python3 scripts/qc_xgmml.py path/to/linkset.xgmml          # one file
python3 scripts/qc_xgmml.py '*.xgmml'                       # a glob
python3 scripts/qc_xgmml.py --allow-empty --min-edges 1 a.xgmml b.xgmml

# Report, and optionally require, BridgeDb cross-reference coverage
python3 scripts/qc_xgmml.py --mapped-type gene --min-mapped-frac 0.9 linkset.xgmml
```
`--mapped-type` reports what fraction of nodes of that type gained more than one
identifier, and `--min-mapped-frac` fails the file below that fraction. This
catches a BridgeDb mapper that failed to load, which leaves every node carrying
only its input id — structurally valid, but useless for matching a network.

Exit code is non-zero if any file fails. The same check runs in CI via `.github/workflows/qc-xgmml.yml`, which can be triggered manually (point it at a `url` or repo `pattern`) or called from another workflow (`workflow_call`, e.g. on a built artifact).
