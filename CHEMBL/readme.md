# ChEMBL mechanism-of-action linkset

Compound → target linkset for human, built from ChEMBL compound's mechanism-of-action (MoA)
records.

| | |
| --- | --- |
| Source | https://www.ebi.ac.uk/chembl/ (via `chembl_webresource_client`) |
| License | CC BY-SA 3.0 |
| DOI | not yet deposited |
| Workflow | `.github/workflows/Chembl_MoA.yml` |
| Species | human |

## How it works

`CHEMBL_MoA.py` pulls human mechanism records and writes `input.txt`. A ChEMBL target can also be a protein complex; each complex was split into multiple targets.

```bash
python3 CHEMBL/CHEMBL_MoA.py                 # takes an hour
python3 CHEMBL/CHEMBL_MoA.py --limit 40      # test run, takes minutes
```
The file contains the following columns.

## Compound attributes

| Column | Description |
| --- | --- |
| `molecule_chembl_id` | ChEMBL molecule identifier |
| `compound_name` | Preferred compound name |
| `inchi_key` | Standard InChIKey |
| `pubchem_sid` | First available PubChem substance identifier |
| `max_phase` | Highest clinical development phase |
| `molecule_type` | ChEMBL molecule type |
| `parent_drug` | Parent ChEMBL molecule identifier |
| `first_approval` | First approval year |
| `indication_class` | Therapeutic or indication class |

## Edge attributes

| Column | Description |
| --- | --- |
| `mechanism_id` | ChEMBL mechanism identifier |
| `mechanism_of_action` | Description of the mechanism of action |
| `action_type` | Type of action, such as inhibitor or agonist |
| `direct_interaction` | Whether a direct interaction is reported |
| `molecular_mechanism` | Whether the record describes a molecular mechanism |
| `disease_efficacy` | Whether disease efficacy is reported |
| `selectivity_comment` | Selectivity information |
| `binding_site_comment` | Binding-site information |
| `ref_id` | Identifier of the first mechanism reference |
| `ref_type` | Reference type |
| `ref_url` | Reference URL |
| `row_end` | Placeholder used to prevent trailing empty columns from being removed by the Linkset Creator |

## Target attributes

| Column | Description |
| --- | --- |
| `target_chembl_id` | ChEMBL target identifier |
| `target_name` | Preferred target name |
| `target_type` | ChEMBL target type |
| `target_organism` | Target organism |
| `uniprot_accessions` | UniProt accession for the target component |
| `gene_symbol` | Human gene symbol |
| `hgnc_id` | HGNC identifier |

## Versioning

The release is read from ChEMBL's `status.json` at run time and flows into the
config, the graph version and the Zenodo metadata, so nothing needs editing when
ChEMBL releases. The config is generated from `chembl_human_MoA.config.template`
and is gitignored — do not hand-edit it.

## Release

Not yet deposited.
