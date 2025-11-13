# linkset-automation
Under construction: This repo automatically generates CyTargetLinker linksets from different resources, starting with WikiPathways.

<img width="692" height="347" alt="image" src="https://github.com/user-attachments/assets/100e64c0-d6e5-47c9-8686-7be01c39c257" />

## (Intended) Workflow
1. Run the bridgeDb update script to update bridgeDb versions for all config files.
2. Check the download link for the resource you want to update (check if an update is required) and copy it into the (data preprocessing) script
3. Run the preprocessing script to get the input.csv file for this version. The config file is available and should not change between versions.
4. Run the LinksetCreator
5. Do the quality control procedure with the obtained xgmml file.
6. If positive, run the linkset upload-to-storage script
