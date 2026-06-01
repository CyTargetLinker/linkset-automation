import csv

import gseapy as gp
import requests

# WikiPathways GMT release
VERSION = "20260510"
BASE_URL = "https://data.wikipathways.org"

# (WikiPathways species token in the GMT filename, CyTargetLinker short code)
SPECIES = [
    ("Homo_sapiens", "hsa"),
    ("Mus_musculus", "mmu"),
    ("Rattus_norvegicus", "rno"),
]


def process(species_token, code):
    gmt = f"wikipathways-{VERSION}-gmt-{species_token}.gmt"
    url = f"{BASE_URL}/{VERSION}/gmt/{gmt}"
    print(f"{species_token} ({code}): downloading {url}")
    response = requests.get(url)
    response.raise_for_status()
    with open(gmt, "w") as f:
        f.write(response.text)

    gene_sets = gp.read_gmt(gmt)

    out = f"input_{code}.txt"
    rows = 0
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")  # Tab-separated
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
        # GMT set name: "Pathway name%WikiPathways_<date>%WP####%Species"
        for pathway_name, genes in gene_sets.items():
            parts = pathway_name.split("%")
            p_name = parts[0]
            gmt_version = parts[1]
            pathway_id = parts[2]
            species = parts[3]
            for gene in genes:
                writer.writerow(
                    [p_name, pathway_id, len(genes), gene, gene, gmt_version, species]
                )
                rows += 1
    print(f"  wrote {out} ({rows} interactions)")


def main():
    for species_token, code in SPECIES:
        process(species_token, code)


if __name__ == "__main__":
    main()
