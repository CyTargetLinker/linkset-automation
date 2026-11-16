The CTD.py script downloads chemical–gene interaction data from the Comparative Toxicogenomics Database (CTD) and creates a separate output file for each supported species. CTD updates its data regularly, usually early each month.

The workflow includes species with interaction data and available BridgeDb mappings for their gene identifiers. Currently supported species are Caenorhabditis elegans, Danio rerio, Gallus gallus, Homo sapiens, and Mus musculus.

Records are grouped by chemical ID and gene ID. For each pair, the script combines the interaction types, interaction actions, and PubMed IDs into one edge, and adds a count of the unique PubMed IDs.

reference:

Davis AP, Wiegers TC, Sciaky D, Barkalow F, Wyatt B, Wiegers J, McMorran R, Abrar S, Mattingly CJ Comparative Toxicogenomics Database: update 2027. Nucleic Acids Res 2026 Oct 23
