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

Size is the binding constraint. WikiPathways human is 40k edges (32 MB) and
works; the GO linkset was removed at 226k edges (183 MB). What decides it is
edges per gene, not rows in the source.

### Ready to build

Openly licensed, still curated upstream, and bounded as downloaded. Counts are
human, measured from the current release.

| Database | Interactions | Size | License | Source |
| --- | --- | --- | --- | --- |
| SIGNOR | signed, directed causal signalling | 19.6k edges, 6.4k proteins, median degree 2 | CC BY 4.0 | https://signor.uniroma2.it/getData.php?organism=9606 |
| Guide to PHARMACOLOGY | ligand–target | 17.5k pairs, 1.9k targets, 9.6k ligands | ODbL, contents CC BY-SA 4.0 | https://www.guidetopharmacology.org/DATA/interactions.csv |
| Complex Portal | complex–component | 9.6k edges, 2.5k complexes | CC0 | https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/ |
| Ensembl Compara | gene–ortholog | 15.8k one-to-one pairs for human–rat | open | https://ftp.ensembl.org/pub/current/tsv/ensembl-compara/homologies/ |

SIGNOR is one TSV request per organism and carries an effect (up- or
down-regulates) on every edge; it is the only directed, signed resource here.
Guide to PHARMACOLOGY is share-alike, so its deposit cannot be CC BY. The
Complex Portal ComplexTAB files hold the curated complexes only — the hu.MAP
predicted set is not in them, and there is one file per taxon. Ensembl
homologies are filtered to `homology_type == ortholog_one2one`, which holds the
linkset to one edge per gene per species pair and lets a network built for one
species extend into another; the filename carries the release number, so the
release is resolved at run time as WikiPathways already does.

### Needs a decision first

Worth having, but not usable as downloaded.

| Database | Interactions | Decision | Source |
| --- | --- | --- | --- |
| HPO | gene–disease, gene–phenotype | `genes_to_disease` (~3.5 per gene), not `genes_to_phenotype` (~78 per gene) | https://github.com/obophenotype/human-phenotype-ontology/releases |
| DISEASES | disease–gene | confidence cutoff for the text-mining channel; the knowledge and experiments channels are bounded as they are | https://diseases.jensenlab.org/Downloads |
| CTD | chemical–gene, chemical–disease, chemical–pathway | ~3.8M curated chemical–gene interactions, so needs filtering by organism and interaction action. The terms of use do not say whether a derived file may be redistributed | https://ctdbase.org/downloads |
| Open Targets | target–disease | CC0 and quarterly, but 17M associations — needs a score threshold that can be defended | https://platform.opentargets.org/downloads |
| miRTarBase | microRNA–target gene | the strong-evidence subset (reporter assay, Western blot), not the full 3.8M | https://awi.cuhk.edu.cn/miRTarBase/ |
| TarBase | microRNA–target gene | CC BY 4.0, better licensed than miRTarBase; filter to the low-throughput methods | https://dianalab.e-ce.uth.gr/tarbasev9 |
| AOP-Wiki | key event–gene | no gene table upstream; would be derived through the KE–WikiPathways mappings, and coverage is partial | https://aopwiki.org |

### Not worth a workflow

Upstream no longer releases, so a scheduled build gains nothing over a one-off.

| Database | Last release |
| --- | --- |
| TRRUST | 2018 |
| CollecTRI | 2023 |
| CORUM | 2024 |
| HuRI | 2020 |
| SIDER | 2015 |
| NCBI HomoloGene | retired 2024 |

### Too large

These would repeat the GO problem: STRING (usable only above a 0.9 cutoff),
IntAct and BioGRID in full, RegNetwork (11M interactions), hTFtarget (3.2M),
InterPro/Pfam domains, Rhea, CellMarker (~32 cell types per gene), and predicted
microRNA targets such as miRDB and TargetScan (500–600 targets per microRNA).

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

- [`update-bridgedb.yml`](.github/workflows/update-bridgedb.yml) — refreshes the
  cached BridgeDb identifier-mapping bundle that every linkset build restores.
- [`qc-xgmml.yml`](.github/workflows/qc-xgmml.yml) — validates a built XGMML
  linkset; callable from another workflow or triggered manually.
