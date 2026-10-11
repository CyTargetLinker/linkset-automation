#!/usr/bin/env bash
# Emit the Zenodo metadata JSON for the SIGNOR linkset deposition.
#
#   bash signor/zenodo-metadata.sh <version>
#
# Kept in one place because the workflow needs it twice, for the bootstrap draft
# and for every later version, and the two drifting apart would show up only as
# a published record with the wrong metadata.
set -euo pipefail

version="${1:?usage: zenodo-metadata.sh <version>}"

description="<p>This dataset contains linkset networks of causal interactions from the SIGNOR ${version} release, in XGMML format, for human.</p>\
<p>signor_hsa.xgmml holds protein to protein regulation and signor_hsa_chemical.xgmml holds ChEBI chemical to protein regulation. Each edge is one source and target pair, carrying the SIGNOR effects, a sign (+, - or +;- where the evidence disagrees), the mechanisms, whether the interaction is direct, and the supporting PMIDs and SIGNOR ids.</p>\
<p>Proteins carry UniProt accessions expanded to their BridgeDb cross-references (Ensembl, Entrez Gene, HGNC); chemicals carry ChEBI ids expanded through the BridgeDb metabolite mapping. Complexes, protein families, phenotypes, stimuli and microRNAs are not included, since they cannot be matched to a network through BridgeDb.</p>\
<p>SIGNOR data is distributed under CC BY 4.0.</p>"

jq -n \
  --arg title "SIGNOR Linksets" \
  --arg version "$version" \
  --arg date "$(date +%Y-%m-%d)" \
  --arg description "$description" \
  --argjson communities '[{"identifier": "cytargetlinker"}]' \
  --argjson creators '[
    {"name": "Kutmon, Martina", "affiliation": "Maastricht University", "orcid": "0000-0002-7699-8191"},
    {"name": "Martens, Marvin", "affiliation": "Maastricht University", "orcid": "0000-0003-2230-0840"},
    {"name": "Ehrhart, Friederike", "affiliation": "Maastricht University", "orcid": "0000-0002-7770-620X"}
  ]' \
  --argjson related_identifiers '[
    {"identifier": "10.1093/nar/gkac883", "relation": "isDerivedFrom", "resource_type": "publication-article", "scheme": "doi"},
    {"identifier": "10.1093/nar/gkad960", "relation": "references", "resource_type": "publication-article", "scheme": "doi"},
    {"identifier": "10.12688/f1000research.14613.2", "relation": "references", "resource_type": "publication-article", "scheme": "doi"}
  ]' \
  '{ metadata: {
       title: $title, version: $version, publication_date: $date,
       description: $description, access_right: "open",
       upload_type: "dataset", language: "eng",
       license: "cc-by-4.0", creators: $creators,
       related_identifiers: $related_identifiers,
       references: ["Lo Surdo P et al (2023) SIGNOR 3.0, the SIGnaling network open resource 3.0: 2022 update. Nucleic Acids Research 51:D631-D637. https://doi.org/10.1093/nar/gkac883"],
       communities: $communities } }'
