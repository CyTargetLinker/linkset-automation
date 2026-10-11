#!/usr/bin/env bash
# Emit the Zenodo metadata JSON for the AOP-Wiki linkset deposition.
#
#   bash aopwiki/zenodo-metadata.sh <version>
#
# Kept in one place because the workflow needs it twice, for the bootstrap draft
# and for every later version, and the two drifting apart would show up only as
# a published record with the wrong metadata.
set -euo pipefail

version="${1:?usage: zenodo-metadata.sh <version>}"

description="<p>This dataset contains linkset networks derived from the AOP-Wiki, via AOP-Wiki RDF release ${version}, in XGMML format.</p>\
<p>aopwiki_ke_gene.xgmml links key events to the human genes text-mined from their titles and descriptions by the AOP-Wiki RDF build. These links are not curated; each edge records whether it was found by named-entity recognition, by regular expression, or both. Genes carry HGNC ids expanded to their BridgeDb cross-references (Ensembl, Entrez Gene, HGNC symbol, UniProt).</p>\
<p>aopwiki_stressor_mie.xgmml links chemical stressors to molecular initiating events. The AOP-Wiki attaches stressors to whole AOPs, so each chemical is linked to the molecular initiating events of the AOPs that list it, and each edge records those AOPs. Chemicals are identified by InChIKey, with CAS number and DSSTox id as attributes.</p>\
<p>Only AOPs that carry a license in the AOP-Wiki are used. AOP-Wiki content is distributed under CC BY-SA 4.0, so these linksets are too.</p>"

jq -n \
  --arg title "AOP-Wiki Linksets" \
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
    {"identifier": "10.1089/aivt.2021.0010", "relation": "isDerivedFrom", "resource_type": "publication-article", "scheme": "doi"},
    {"identifier": "10.1093/nar/gkad960", "relation": "references", "resource_type": "publication-article", "scheme": "doi"},
    {"identifier": "10.12688/f1000research.14613.2", "relation": "references", "resource_type": "publication-article", "scheme": "doi"}
  ]' \
  '{ metadata: {
       title: $title, version: $version, publication_date: $date,
       description: $description, access_right: "open",
       upload_type: "dataset", language: "eng",
       license: "cc-by-sa-4.0", creators: $creators,
       related_identifiers: $related_identifiers,
       references: ["Martens M, Evelo CT, Willighagen EL (2022) Providing Adverse Outcome Pathways from the AOP-Wiki in a Semantic Web Format to Increase Usability and Accessibility of the Content. Applied In Vitro Toxicology 8(1):2-13. https://doi.org/10.1089/aivt.2021.0010"],
       communities: $communities } }'
