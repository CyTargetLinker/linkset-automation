import csv
import gzip
import tempfile
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import requests


URL = "https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz"

OUTPUT_DIR = Path("CTD_outputs")
DATE = datetime.now().strftime("%Y-%m-%d")

SPECIES = {
    "caenorhabditis elegans": "caenorhabditis_elegans",
    "danio rerio": "danio_rerio",
    "gallus gallus": "gallus_gallus",
    "homo sapiens": "homo_sapiens",
    "mus musculus": "mus_musculus",
}

COLUMNS = [
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "GeneSymbol",
    "GeneID",
    "GeneForms",
    "Organism",
    "OrganismID",
    "Interaction",
    "InteractionActions",
    "PubMedIDs",
    "PubMedIDCount",
]


def split_values(value):
    """Split a CTD pipe-delimited field and preserve first-seen order."""
    return [part.strip() for part in value.split("|") if part.strip()]


def merge_values(existing, value):
    """Append unseen pipe-delimited values to an ordered collection."""
    for item in split_values(value):
        existing.setdefault(item, None)


def download_file(destination):
    print(f"Downloading {URL}")
    downloaded = 0

    with requests.get(URL, stream=True, timeout=(30, 300)) as response:
        response.raise_for_status()
        with destination.open("wb") as output:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    output.write(chunk)
                    downloaded += len(chunk)
                    print(
                        f"\rDownloaded: {downloaded / 1024 / 1024:.1f} MB",
                        end="",
                        flush=True,
                    )
    print("\nDownload complete.")


def split_by_species(source_path):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    field_names = None
    indices = None
    organism_index = None
    grouped = {name: OrderedDict() for name in SPECIES}
    total_rows = 0
    kept_rows = 0

    with gzip.open(source_path, mode="rt", encoding="utf-8-sig", newline="") as source:
        def data_lines():
            nonlocal field_names
            for line in source:
                if line.startswith("#"):
                    parsed = next(csv.reader([line[1:].strip()]))
                    if parsed and parsed[0].strip() == "ChemicalName":
                        field_names = [name.strip() for name in parsed]
                    continue
                if line.strip():
                    yield line

        for row in csv.reader(data_lines()):
            if indices is None:
                if field_names is None:
                    raise ValueError("Could not find the CTD column header.")
                required = [name for name in COLUMNS[:-1] if name not in field_names]
                if required:
                    raise ValueError("Missing CTD columns: " + ", ".join(required))
                indices = {name: field_names.index(name) for name in COLUMNS[:-1]}
                organism_index = indices["Organism"]

            total_rows += 1
            if len(row) != len(field_names):
                raise ValueError(
                    f"Unexpected field count at CTD data row {total_rows}: "
                    f"expected {len(field_names)}, found {len(row)}."
                )

            species = " ".join(row[organism_index].split()).casefold()
            if species not in SPECIES:
                continue

            chemical_id = row[indices["ChemicalID"]].strip()
            gene_id = row[indices["GeneID"]].strip()
            if not chemical_id or not gene_id:
                continue

            groups = grouped[species]
            key = (chemical_id, gene_id)
            if key not in groups:
                fixed = {
                    name: row[indices[name]].replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()
                    for name in COLUMNS[:8]
                }
                groups[key] = {
                    "fixed": fixed,
                    "Interaction": OrderedDict(),
                    "InteractionActions": OrderedDict(),
                    "PubMedIDs": OrderedDict(),
                }

            group = groups[key]
            for name in ("Interaction", "InteractionActions", "PubMedIDs"):
                merge_values(group[name], row[indices[name]])
            kept_rows += 1

            if total_rows % 100_000 == 0:
                print(f"Processed {total_rows:,} rows; retained {kept_rows:,}")

    if total_rows == 0:
        raise ValueError("No interaction records found.")

    print(f"\nProcessed: {total_rows:,} source rows")
    print(f"Retained: {kept_rows:,} rows before pair aggregation")
    print(f"Output folder: {OUTPUT_DIR.resolve()}")

    for species, token in SPECIES.items():
        groups = grouped[species]
        if not groups:
            print(f"  {token}: no matching records")
            continue

        filename = OUTPUT_DIR / f"CTD_{DATE}_{token}.txt"
        with filename.open("w", encoding="utf-8", newline="") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(COLUMNS)
            for group in groups.values():
                values = [group["fixed"][name] for name in COLUMNS[:8]]
                for name in ("Interaction", "InteractionActions", "PubMedIDs"):
                    values.append("|".join(group[name].keys()))
                values.append(str(len(group["PubMedIDs"])))
                writer.writerow(values)

        print(
            f"  {filename.name}: {len(groups):,} unique ChemicalID-GeneID pairs"
        )


def main():
    with tempfile.TemporaryDirectory(prefix="ctd_download_") as temp_dir:
        archive_path = Path(temp_dir) / "CTD_chem_gene_ixns.csv.gz"
        download_file(archive_path)
        split_by_species(archive_path)


if __name__ == "__main__":
    main()
