import csv
import gzip
import tempfile
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
]


def download_file(destination):
    print(f"Downloading {URL}")
    downloaded = 0

    with requests.get(
        URL,
        stream=True,
        timeout=(30, 300),
    ) as response:
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

    handles = {}
    writers = {}
    counts = {name: 0 for name in SPECIES}

    field_names = None
    column_indices = None
    organism_index = None
    total_rows = 0
    kept_rows = 0

    try:
        with gzip.open(
            source_path,
            mode="rt",
            encoding="utf-8-sig",
            newline="",
        ) as source:

            def data_lines():
                nonlocal field_names

                for line in source:
                    if line.startswith("#"):
                        header = line[1:].strip()
                        parsed = next(csv.reader([header]))

                        if parsed and parsed[0].strip() == "ChemicalName":
                            field_names = [
                                name.strip() for name in parsed
                            ]
                        continue

                    if line.strip():
                        yield line

            for row in csv.reader(data_lines()):
                if column_indices is None:
                    if field_names is None:
                        raise ValueError(
                            "Could not find the CTD column header."
                        )

                    missing = [
                        name for name in COLUMNS
                        if name not in field_names
                    ]
                    if missing:
                        raise ValueError(
                            "Missing columns: " + ", ".join(missing)
                        )

                    column_indices = [
                        field_names.index(name) for name in COLUMNS
                    ]
                    organism_index = field_names.index("Organism")

                total_rows += 1

                if len(row) != len(field_names):
                    raise ValueError(
                        f"Unexpected field count in data row {total_rows}: "
                        f"expected {len(field_names)}, found {len(row)}."
                    )

                # Exact species match, ignoring case and extra whitespace.
                species = " ".join(
                    row[organism_index].split()
                ).casefold()

                if species in SPECIES:
                    if species not in writers:
                        filename = (
                            f"CTD_{DATE}_{SPECIES[species]}.txt"
                        )
                        output_path = OUTPUT_DIR / filename

                        handle = output_path.open(
                            mode="w",
                            encoding="utf-8",
                            newline="",
                        )
                        handles[species] = handle

                        writer = csv.writer(
                            handle,
                            delimiter="\t",
                            lineterminator="\n",
                        )
                        writer.writerow(COLUMNS)
                        writers[species] = writer

                    values = [
                        row[index]
                        .replace("\t", " ")
                        .replace("\r", " ")
                        .replace("\n", " ")
                        for index in column_indices
                    ]

                    writers[species].writerow(values)
                    counts[species] += 1
                    kept_rows += 1

                if total_rows % 100_000 == 0:
                    print(
                        f"Processed {total_rows:,} rows; "
                        f"kept {kept_rows:,}"
                    )

        if total_rows == 0:
            raise ValueError("No interaction records found.")

    finally:
        for handle in handles.values():
            handle.close()

    print(f"\nProcessed: {total_rows:,} records")
    print(f"Kept: {kept_rows:,} records")
    print(f"Output folder: {OUTPUT_DIR.resolve()}")

    for species, filename_species in SPECIES.items():
        if counts[species]:
            print(
                f"  CTD_{DATE}_{filename_species}.txt: "
                f"{counts[species]:,} rows"
            )
        else:
            print(f"  {filename_species}: no matching records")


def main():
    with tempfile.TemporaryDirectory(prefix="ctd_download_") as temp_dir:
        archive_path = Path(temp_dir) / "CTD_chem_gene_ixns.csv.gz"
        download_file(archive_path)
        split_by_species(archive_path)


if __name__ == "__main__":
    main()
