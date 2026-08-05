#!/usr/bin/env bash
# Emit the Zenodo metadata JSON for the Reactome linkset deposition.
#
#   bash reactome/zenodo-metadata.sh <version> <n_species>
#
# Kept in one place because the workflow needs it twice, for the bootstrap draft
# and for every later version, and the two drifting apart would show up only as
# a published record with the wrong metadata.
set -euo pipefail

version="${1:?usage: zenodo-metadata.sh <version> <n_species>}"
n_species="${2:?usage: zenodo-metadata.sh <version> <n_species>}"

description="<p>This dataset contains linkset networks of pathway-gene interactions from Reactome release ${version}, in XGMML format, for ${n_species} species.</p>\
<p>Built from Reactome's NCBI2Reactome table (lowest-level pathways), so gene nodes carry NCBI Entrez Gene IDs expanded to their BridgeDb cross-references, and are labelled with the official gene symbol from NCBI Gene.</p>\
<p>Genes NCBI has not named are omitted rather than kept under a bare identifier. These are almost entirely pathogen-side genes that Reactome files under the host species in its infectious-disease pathways, and none of them map in BridgeDb; the cost is that a small number of purely pathogen pathways are absent.</p>\
<p>Reactome data is dedicated to the public domain under CC0.</p>"

jq -n \
  --arg title "Reactome Linksets" \
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
    {"identifier": "10.1093/nar/gkad1025", "relation": "isDerivedFrom", "resource_type": "publication-article", "scheme": "doi"},
    {"identifier": "10.1093/nar/gkad960", "relation": "references", "resource_type": "publication-article", "scheme": "doi"},
    {"identifier": "10.12688/f1000research.14613.2", "relation": "references", "resource_type": "publication-article", "scheme": "doi"}
  ]' \
  '{ metadata: {
       title: $title, version: $version, publication_date: $date,
       description: $description, access_right: "open",
       upload_type: "dataset", language: "eng",
       license: "cc-zero", creators: $creators,
       related_identifiers: $related_identifiers,
       references: ["Milacic M et al (2024) The Reactome Pathway Knowledgebase 2024. Nucleic Acids Research 52:D672-D678. https://doi.org/10.1093/nar/gkad1025"],
       communities: $communities } }'
