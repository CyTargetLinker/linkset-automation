# Linkset status

Which source databases this repository automates, and how far each one has got.

Status here describes the **automation pipeline**, not whether a CyTargetLinker
linkset exists at all. Several databases below already have linksets on the
[CyTargetLinker website](https://cytargetlinker.github.io/pages/linksets) that
were built by hand before this repository existed; they still appear as
candidates because the build is not yet automated.

Deposit details, meaning persistent DOIs, licenses and per-resource build
notes, live in the resource table in [`README.md`](README.md#linkset-resources)
rather than being repeated here.

## Complete automation

| Resource | Interactions | Workflow | Deposit | Lastest version | License
| --- | --- | --- | --- | --- |--- |
| [WikiPathways](https://www.wikipathways.org) | pathway–gene, 13 species | [`create-linkset-wikipathways.yml`](.github/workflows/create-linkset-wikipathways.yml) | [WikiPathways Linksets](https://zenodo.org/records/21806899) | Version 20260710 | [`CC0`](https://creativecommons.org/2008/04/16/cc0-betadiscussion-draft-2/) 
| [TFLink](https://tflink.net) | transcription factor–target gene, 5 species | [`create-linkset-tflink.yml`](.github/workflows/create-linkset-tflink.yml) | [TFLink Linksets](https://zenodo.org/records/21804829) | Version v1.0-SS | [`CC BY-NC 4.0`](https://creativecommons.org/licenses/by-nc/4.0/)|

## Automated, not yet deposited

The linkset generation is complete, but not yet deposited on Zenodo.

| Resource | Interactions | Workflow | Outstanding | Lastest version | License | 
| --- | --- | --- | --- | --- |--- |
| [Orphanet](https://www.orphadata.com/) | gene–rare disease | [`Create_OrphanetLinkset.yml`](.github/workflows/Create_OrphanetLinkset.yml) | [`Orphanet/readme.md`](Orphanet/readme.md) | July 2026 | [`CC BY 4.0`](https://creativecommons.org/licenses/by/4.0/deed.en) |
| [ChEMBL](https://www.ebi.ac.uk/chembl/) | compound–protein target (mechanism of action) | [`Chembl_MoA.yml`](.github/workflows/Chembl_MoA.yml) | [`CHEMBL/readme.md`](CHEMBL/readme.md) | August 2026 | [`CC0`](https://creativecommons.org/2008/04/16/cc0-betadiscussion-draft-2/) | 
| [Reactome](https://reactome.org/) | pathway–gene, human | [`create-linkset-reactome.yml`](.github/workflows/create-linkset-reactome.yml) | Zenodo bootstrap; whether orthology-inferred species should ship |  Na | [`CC BY 4.0`](https://creativecommons.org/licenses/by/4.0/deed.en) |
| [AOP-Wiki](https://aopwiki.org/) | key event–gene (text-mined), chemical stressor–molecular initiating event | [`create-linkset-aopwiki.yml`](.github/workflows/create-linkset-aopwiki.yml) | Zenodo bootstrap |  Na| [`CC BY-SA`](https://creativecommons.org/licenses/by-sa/2.0/)
| [SIGNOR](https://signor.uniroma2.it/) | protein–protein causal signalling, human | [`create-linkset-signor.yml`](.github/workflows/create-linkset-signor.yml) | Zenodo bootstrap; chemical–protein waits on [linksetCreator#5](https://github.com/CyTargetLinker/linksetCreator/issues/5) | Na | Na |

## Candidates

Databases we would develop automation. The consideration includes:
1. Regular updates from the database;
2. Open-licensed resources;
3. Data quality.

### Ready to build

Openly licensed, still curated upstream, and bounded as downloaded. Counts are
human, measured from the current release.

| Database | Interactions | update frequency | License | Source |
| --- | --- | --- | --- | --- |
| OpenTargets | gene_disease_drug | Every three month on 26th |  [`CC0 1.0`](https://creativecommons.org/publicdomain/zero/1.0/?ref=chooser-v1)  | https://www.guidetopharmacology.org/DATA/interactions.csv |


### Create but not update

Worth having, but not automation needed for regular updating is not clear.

| Database | Interactions | License| Source |
| --- | --- | --- | --- |
| CollecTRI | transcription factor–gene | freely available | https://diseases.jensenlab.org/Downloads |


## License limitation

Redistributing the built linkset is not permitted by the upstream license. The
aim is to provide a script that lets users build the linkset themselves from
their own copy of the source data, without shipping the data files.

| Database | Interactions | Source |
| --- | --- | --- |
| DrugBank | drug–target | https://go.drugbank.com/ |
| DisGeNET | gene–disease associations | https://www.disgenet.com/ |
| MalaCards | gene–disease associations | https://www.malacards.org/ |

## Supporting workflows

Not linksets themselves, but part of the build:

- [`update-bridgedb.yml`](.github/workflows/update-bridgedb.yml) refreshes the
  cached BridgeDb identifier-mapping bundle that every linkset build restores.
- [`qc-xgmml.yml`](.github/workflows/qc-xgmml.yml) validates a built XGMML
  linkset; callable from another workflow or triggered manually.
