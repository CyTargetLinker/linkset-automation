import csv

import gseapy as gp
import mygene
import requests

# gmt version
version = "20251110"
gmt = "wikipathways-" + version + "-gmt-Homo_sapiens.gmt"

# Download from URL
url = "https://data.wikipathways.org/" + version + "/gmt/" + gmt
response = requests.get(url)

# Save and read
with open(gmt, "w") as f:
    f.write(response.text)

# Read GMT file
gene_sets = gp.read_gmt(gmt)

with open("input.txt", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")  # Tab-separated

    # Write header
    writer.writerow(
        [
            "pathway_name",
            "pathway_id",
            "gene_count",
            "gene_id",
            "gene_symbol",
            "version",
            "species",
        ]
    )

    # Access data
    for pathway_name, genes in gene_sets.items():
        parts = pathway_name.split("%")
        p_name = parts[0]  # "Clock controlled autophagy in bone metabolism"
        version = parts[1]  # "WikiPathways_20251110"
        pathway_id = parts[2]  # "WP5205"
        species = parts[3]  # "Homo sapiens"
        
        for gene in genes:
            writer.writerow(
                [p_name, pathway_id, len(genes), gene, gene, version, species]
            )
