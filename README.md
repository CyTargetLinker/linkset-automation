# linkset-automation
Under construction: This repo automatically generates CyTargetLinker linksets from different resources, starting with WikiPathways.

<img width="692" height="347" alt="image" src="https://github.com/user-attachments/assets/100e64c0-d6e5-47c9-8686-7be01c39c257" />

## (Intended) Workflow
1. Run the bridgeDb update script to update bridgeDb versions for all config files.
2. Check the download link for the resource you want to update (check if an update is required) and copy it into the (data preprocessing) script
3. Run the preprocessing script to get the input.csv file for this version. The config file is available and should not change between versions.
4. Run the LinksetCreator
5. Do the quality control procedure with the obtained xgmml file (see **QC** below).
6. If positive, run the linkset upload-to-storage script
7. Update the CTL website to add the new linkset download link and description

## QC
`scripts/qc_xgmml.py` validates that an XGMML linkset is well-formed and has the structure every CyTargetLinker linkset carries (a `<graph>`, nodes with an `identifiers` list + `type`, and edges whose `source`/`target` reference existing nodes plus `datasource`/`interaction`). Empty graphs fail by default.

```bash
python3 scripts/qc_xgmml.py path/to/linkset.xgmml          # one file
python3 scripts/qc_xgmml.py '*.xgmml'                       # a glob
python3 scripts/qc_xgmml.py --allow-empty --min-edges 1 a.xgmml b.xgmml
```
Exit code is non-zero if any file fails. The same check runs in CI via `.github/workflows/qc-xgmml.yml`, which can be triggered manually (point it at a `url` or repo `pattern`) or called from another workflow (`workflow_call`, e.g. on a built artifact).
