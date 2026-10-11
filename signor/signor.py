#!/usr/bin/env python3
"""Preprocess SIGNOR's human causal network into linkset-creator input files.

Writes two linkset inputs, because linkset-creator fixes the node type, input
syscode and BridgeDb file per side of a config, so proteins and chemicals cannot
share a side:

    input_signor_protein.txt                    protein -> protein
    input_signor_chemical.txt                   ChEBI chemical -> protein
    signor/config/signor_{protein,chemical}.config

plus build/version.txt.

Content comes from getData.php?organism=9606, which is the live database sliced
to human proteins. It is not filtered on TAX_ID: that column is the species the
experiment was done in, and the human file carries thousands of rows on human
proteins with TAX_ID -1 (unspecified) or 10090 (mouse cells). The monthly release
file is used only for its name, which dates the build (Oct2026_release.txt ->
2026-10); its content mixes every species.

    ENTITYA TYPEA IDA DATABASEA ENTITYB TYPEB IDB DATABASEB EFFECT ...

Rows are collapsed to one edge per (source, target), joining the effects,
mechanisms, PMIDs and SIGNOR ids of the rows behind it. One edge per pair is not
a choice: linkset-creator keeps only the first row for a pair and drops the rest.
A pair whose rows disagree gets the sign "+;-".

    python3 signor/signor.py
    python3 signor/signor.py --version 2026-10      # override the version label only
"""

import argparse
import csv
import os
import pathlib
import re
import sys
import time

import requests

DATA_URL = "https://signor.uniroma2.it/getData.php?organism=9606"
RELEASE_URL = "https://signor.uniroma2.it/releases/getLatestRelease.php"
DATA_FILE = "signor_9606.tsv"

CONFIG_DIR = "signor/config"
BUILD_DIR = "build"

# Column order of getData.php, which sends no header. It matches the header of
# the monthly release file; getData appends a SCORE column and a trailing tab.
COLUMNS = [
    "ENTITYA", "TYPEA", "IDA", "DATABASEA", "ENTITYB", "TYPEB", "IDB", "DATABASEB",
    "EFFECT", "MECHANISM", "RESIDUE", "SEQUENCE", "TAX_ID", "CELL_DATA",
    "TISSUE_DATA", "MODULATOR_COMPLEX", "TARGET_COMPLEX", "MODIFICATIONA",
    "MODASEQ", "MODIFICATIONB", "MODBSEQ", "PMID", "DIRECT", "NOTES", "ANNOTATOR",
    "SENTENCE", "SIGNOR_ID",
]

# (name, source database, target database, minimum edges). A build below the
# floor means the download or the column layout changed, not that SIGNOR shrank.
LINKSETS = [
    ("protein", "UNIPROT", "UNIPROT", 10000),
    ("chemical", "ChEBI", "UNIPROT", 1000),
]

INPUT_COLUMNS = [
    "source_id", "source_label", "target_id", "target_label",
    "effect", "sign", "mechanism", "direct", "pmid", "signor_id",
]

# Isoform (P05771-2) and processed-chain (P01019-PRO_0000032457) suffixes.
# BridgeDb knows canonical accessions only, so these would never map.
UNIPROT_SUFFIX_RE = re.compile(r"-(\d+|PRO_\d+)$")

MONTHS = {name: f"{i:02d}" for i, name in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], start=1)}
RELEASE_NAME_RE = re.compile(r"(" + "|".join(MONTHS) + r")(\d{4})_release")

# linkset-creator splits rows with Java String.split("\t") rather than a CSV
# parser, so an embedded tab shifts every later column and a newline truncates
# the linkset.
UNSAFE_RE = re.compile(r"[\t\r\n\"]+")

DOWNLOAD_ATTEMPTS = 3


def resolve_version():
    """Return the current SIGNOR release as YYYY-MM.

    Read from the filename the release endpoint offers, without downloading the
    body. YYYY-MM sorts correctly as a string, which the Zenodo gate relies on.
    """
    try:
        with requests.get(RELEASE_URL, stream=True, timeout=60) as response:
            response.raise_for_status()
            disposition = response.headers.get("Content-Disposition", "")
    except requests.RequestException as error:
        sys.exit(f"ERROR: cannot reach {RELEASE_URL}: {error}")
    match = RELEASE_NAME_RE.search(disposition)
    if not match:
        sys.exit(f"ERROR: no release name in {RELEASE_URL} "
                 f"(Content-Disposition: {disposition!r})")
    month, year = match.groups()
    return f"{year}-{MONTHS[month]}"


def download(url, local, attempts=DOWNLOAD_ATTEMPTS):
    """Stream url to local, reusing an existing non-empty copy."""
    if os.path.exists(local) and os.path.getsize(local) > 0:
        print(f"  using existing {local} (delete it to refresh)")
        return local

    tmp = f"{local}.partial"
    delay = 2
    try:
        for attempt in range(1, attempts + 1):
            try:
                with requests.get(url, stream=True, timeout=600) as response:
                    response.raise_for_status()
                    with open(tmp, "wb") as handle:
                        for chunk in response.iter_content(chunk_size=1 << 20):
                            handle.write(chunk)
                os.replace(tmp, local)
                return local
            except requests.RequestException as error:
                print(f"  attempt {attempt}/{attempts} failed: {error}")
                if attempt == attempts:
                    sys.exit(f"ERROR: could not download {url}")
                time.sleep(delay)
                delay *= 2
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def clean(value):
    """Strip surrounding quotes and anything that would break the input file."""
    return UNSAFE_RE.sub(" ", value.strip().strip('"')).strip()


def normalise_id(identifier, database):
    if database == "UNIPROT":
        return UNIPROT_SUFFIX_RE.sub("", identifier)
    return identifier


def sign_of(effects):
    """'+', '-', '+;-' when the rows disagree, or '' when none is signed."""
    signs = set()
    for effect in effects:
        if effect.startswith("up-regulates"):
            signs.add("+")
        elif effect.startswith("down-regulates"):
            signs.add("-")
    return ";".join(sorted(signs))


def parse(path):
    """Yield one dict per row.

    Split on tabs rather than with the csv module: getData leaves fields
    unquoted, and a sentence that happens to start with a quote would make a
    CSV reader swallow the rest of the file.
    """
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(COLUMNS):
                sys.exit(f"ERROR: {path}:{number}: expected at least "
                         f"{len(COLUMNS)} columns, got {len(fields)}")
            yield {name: clean(value) for name, value in zip(COLUMNS, fields)}


def collapse(rows, source_db, target_db):
    """Group rows into edges keyed by (source, target)."""
    edges = {}
    for row in rows:
        if row["DATABASEA"] != source_db or row["DATABASEB"] != target_db:
            continue
        source = normalise_id(row["IDA"], source_db)
        target = normalise_id(row["IDB"], target_db)
        if not source or not target:
            continue
        edge = edges.setdefault((source, target), {
            "source_label": row["ENTITYA"] or source,
            "target_label": row["ENTITYB"] or target,
            "effect": set(), "mechanism": set(), "direct": False,
            "pmid": set(), "signor_id": set(),
        })
        if row["EFFECT"]:
            edge["effect"].add(row["EFFECT"])
        if row["MECHANISM"]:
            edge["mechanism"].add(row["MECHANISM"])
        edge["direct"] |= row["DIRECT"] == "t"
        edge["pmid"].update(p for p in row["PMID"].split("|") if p)
        if row["SIGNOR_ID"]:
            edge["signor_id"].add(row["SIGNOR_ID"])
    return edges


def write_input(name, edges):
    out = f"input_signor_{name}.txt"
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", quoting=csv.QUOTE_NONE,
                            escapechar=None, lineterminator="\n")
        writer.writerow(INPUT_COLUMNS)
        for (source, target), edge in sorted(edges.items()):
            writer.writerow([
                source, edge["source_label"], target, edge["target_label"],
                ";".join(sorted(edge["effect"])), sign_of(edge["effect"]),
                ";".join(sorted(edge["mechanism"])),
                "t" if edge["direct"] else "f",
                ";".join(sorted(edge["pmid"])), ";".join(sorted(edge["signor_id"])),
            ])
    return out


def write_config(name, version):
    template = pathlib.Path(f"signor/signor_{name}.config.template").read_text(encoding="utf-8")
    text = template.format(version=version)
    # ConfigFileReader keeps only lines splitting into exactly two parts on '=',
    # silently discarding the rest as an invalid attribute.
    for line in text.splitlines():
        if line.strip() and line.count("=") != 1:
            sys.exit(f"ERROR: generated config line is not a single key=value pair: {line!r}")
    out = pathlib.Path(CONFIG_DIR) / f"signor_{name}.config"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Build SIGNOR linkset inputs and configs.")
    parser.add_argument(
        "--version", metavar="YYYY-MM",
        help="version label (default: the current SIGNOR release). The content is "
             "always the live database; this does not fetch an older release.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    version = args.version or resolve_version()
    if not re.fullmatch(r"\d{4}-\d{2}", version):
        sys.exit(f"ERROR: version must be YYYY-MM, got {version!r}")
    print(f"SIGNOR {version}")

    download(DATA_URL, DATA_FILE)
    rows = list(parse(DATA_FILE))
    print(f"  {len(rows)} rows")

    for name, source_db, target_db, floor in LINKSETS:
        edges = collapse(rows, source_db, target_db)
        signs = {s: sum(1 for e in edges.values() if sign_of(e["effect"]) == s)
                 for s in ("+", "-", "+;-", "")}
        print(f"{name}: {len(edges)} edges ({signs['+']} up, {signs['-']} down, "
              f"{signs['+;-']} both, {signs['']} unsigned), "
              f"{len({k[0] for k in edges})} sources, "
              f"{len({k[1] for k in edges})} targets")
        if len(edges) < floor:
            sys.exit(f"ERROR: {name}: only {len(edges)} edges, below the floor of {floor}")
        out = write_input(name, edges)
        config = write_config(name, version)
        print(f"  -> {out}, {config}")

    build = pathlib.Path(BUILD_DIR)
    build.mkdir(exist_ok=True)
    (build / "version.txt").write_text(version + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
