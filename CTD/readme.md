The CTD.py script downloads chemical–gene interaction data from the Comparative Toxicogenomics Database (CTD) and creates a separate output file for each supported species. CTD updates its data regularly, usually early each month.

The workflow includes species with interaction data and available BridgeDb mappings for their gene identifiers. Currently supported species are Caenorhabditis elegans, Danio rerio, Gallus gallus, Homo sapiens, and Mus musculus.

Records are grouped by chemical ID and gene ID. For each pair, the script combines the interaction types, interaction actions, and PubMed IDs into one edge, and adds a count of the unique PubMed IDs.

reference:

Allan Peter Davis, Thomas C Wiegers, Daniela Sciaky, Fern Barkalow, Melissa Strong, Brent Wyatt, Jolene Wiegers, Roy McMorran, Sakib Abrar, Carolyn J Mattingly, Comparative Toxicogenomics Database’s 20th anniversary: update 2025, Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D1328–D1334, https://doi.org/10.1093/nar/gkae883
