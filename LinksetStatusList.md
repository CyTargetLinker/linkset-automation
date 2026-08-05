# Linkset status

Which source databases this repository automates, and how far each one has got.

Status here describes the **automation pipeline**, not whether a CyTargetLinker
linkset exists at all. Several databases below already have linksets on the
[CyTargetLinker website](https://cytargetlinker.github.io/pages/linksets) that
were built by hand before this repository existed; they still appear as
candidates because the build is not yet automated.

Deposit details — persistent DOIs, licenses and per-resource build notes — live
in the resource table in [`README.md`](README.md#linkset-resources) rather than
being repeated here.

## In production

The workflow runs green and the output is deposited under a persistent DOI.

| Resource | Interactions | Workflow |
| --- | --- | --- |
| [WikiPathways](https://www.wikipathways.org) | pathway–gene, 13 species | [`create-linkset-wikipathways.yml`](.github/workflows/create-linkset-wikipathways.yml) |
| [TFLink](https://tflink.net) | transcription factor–target gene, 5 species | [`create-linkset-tflink.yml`](.github/workflows/create-linkset-tflink.yml) |

## Automated, not yet deposited

The workflow exists and the linkset builds, but nothing is published under a DOI
yet and the deposit license is still undecided.

| Resource | Interactions | Workflow | Outstanding |
| --- | --- | --- | --- |
| [Orphanet](https://www.orphadata.com/) | gene–rare disease | [`Create_OrphanetLinkset.yml`](.github/workflows/Create_OrphanetLinkset.yml) | Zenodo deposit; [`Orphanet/readme.md`](Orphanet/readme.md) is still a placeholder |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | compound–protein target (mechanism of action) | [`Chembl_MoA.yml`](.github/workflows/Chembl_MoA.yml) | Zenodo deposit |
| [Reactome](https://reactome.org/) | pathway–gene, human | [`create-linkset-reactome.yml`](.github/workflows/create-linkset-reactome.yml) | Zenodo bootstrap; whether orthology-inferred species should ship |

## Candidates

Databases we would like to automate. No workflow yet.

| Database | Interactions | Source |
| --- | --- | --- |
| DISEASES | disease–gene associations | https://diseases.jensenlab.org/Downloads |
| CTD | chemical–target, chemical–disease, chemical–pathway | https://ctdbase.org/downloads |
| TRRUST | transcription factor–target | https://www.grnpedia.org/trrust/ |
| CollecTRI | transcription factor–target | https://github.com/saezlab/collecTRI |

## License limitation

Redistributing the built linkset is not permitted by the upstream license. The
aim is to provide a script that lets users build the linkset themselves from
their own copy of the source data, without shipping the data files.

| Database | Interactions | Source |
| --- | --- | --- |
| DrugBank | drug–target | https://go.drugbank.com/ |
| DisGeNET | gene–disease associations | https://www.disgenet.org/ |

## Supporting workflows

Not linksets themselves, but part of the build:

- [`update-bridgedb.yml`](.github/workflows/update-bridgedb.yml) — refreshes the
  cached BridgeDb identifier-mapping bundle that every linkset build restores.
- [`qc-xgmml.yml`](.github/workflows/qc-xgmml.yml) — validates a built XGMML
  linkset; callable from another workflow or triggered manually.
