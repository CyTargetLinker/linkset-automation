# Orphanet Rare Disease–Gene Linkset

This directory contains the preprocessing code and configuration used to generate the **Orphanet rare disease–gene linkset** for [CyTargetLinker](https://cytargetlinker.github.io/).

## Contents

| File | Description |
| --- | --- |
| `orphadata_diseases_genes.py` | Downloads and combines Orphadata XML products into a tab-separated intermediate file named `input.txt`. |
| `Orphanet.config` | Configuration used by the CyTargetLinker Linkset Creator to convert `input.txt` into an XGMML linkset. |
| `readme.md` | Documentation for the Orphanet linkset workflow. |

## Data sources

The preprocessing script downloads the following Orphadata XML products:

## Generated data

The preprocessing script creates `input.txt`, a tab-separated file with one row per rare disease–gene association.

| Column | Description |
| --- | --- |
| `OrphaID` | Orphanet identifier, formatted as `ORPHA:<code>` |
| `DiseaseName` | Preferred rare disease name |
| `GeneSymbol` | Associated human gene symbol |
| `HGNC_ID` | HGNC identifier supplied by Orphadata, when available |
| `Synonyms` | Semicolon-separated disease synonyms |
| `Source` | Supporting PubMed identifiers and links extracted from the association evidence |
| `Prevalence` | Worldwide prevalence class when available; otherwise the available geographic prevalence class or classes |
| `Inheritance` | Semicolon-separated modes of inheritance |
| `AgeOfOnset` | Semicolon-separated age-of-onset categories |
| `ICD_10` | ICD-10 cross-references |
| `ICD_11` | ICD-11 cross-references with links to the WHO ICD browser |
| `OMIM` | OMIM cross-references with entry links |
| `UMLS` | UMLS cross-references |

Missing values are represented by `-`.

The [Orphanet workflow](https://github.com/CyTargetLinker/linkset-automation/blob/main/.github/workflows/Create_OrphanetLinkset.yml) converts this table into `Orphanet.xgmml`. 

## Requirements

For preprocessing:

- Python 3
- [`requests`](https://pypi.org/project/requests/)
