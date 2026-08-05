#!/usr/bin/env python3
"""Preprocess WikiPathways GMT releases into linkset-creator input files.

For every species listed in wikipathways/species.tsv this writes:

    input_<code>.txt                 tab-delimited pathway-gene pairs
    wikipathways/wp_<code>.config    generated from wp.config.template

plus two files describing the build:

    build/version.txt                the WikiPathways release that was used
    build/manifest.tsv               one row per species ACTUALLY built

The manifest, not species.tsv, is what the GitHub workflow iterates over, so a
species skipped here can never leave the workflow looking for a config, a
BridgeDb download or an XGMML file that was never produced.

Release resolution
------------------
data.wikipathways.org keeps only a rolling 12-month window of monthly releases,
so the release is read from the GMT directory index rather than pinned in
source. Files are then fetched from the resolved /<version>/ path, never from
/current/, so a release flipping mid-run cannot yield a linkset built from two
releases.

Environment overrides (both wired to the workflow's manual-dispatch inputs):
    WP_VERSION   build a specific release instead of the current one
    WP_SPECIES   space-separated codes to build instead of all enabled ones
"""

import csv
import os
import pathlib
import re
import sys

import requests

BASE_URL = "https://data.wikipathways.org"
SPECIES_TSV = "wikipathways/species.tsv"
CONFIG_TEMPLATE = "wikipathways/wp.config.template"
BUILD_DIR = "build"

# wikipathways-<YYYYMMDD>-gmt-<Genus_species>.gmt
GMT_RE = re.compile(r"wikipathways-(\d{8})-gmt-([A-Za-z_]+)\.gmt")

# An 'extra' species with fewer pathways than this is not worth a linkset and is
# skipped; a 'core' species below it means something is wrong upstream.
MIN_PATHWAYS = 5

MANIFEST_COLUMNS = [
    "code", "wp_token", "organism", "bridgedb_species",
    "syscodes_out", "min_mapped_frac", "n_pathways", "n_rows",
]

INPUT_COLUMNS = [
    "pathway_name", "pathway_id", "gene_count",
    "gene_id", "gene_symbol", "version", "species",
]


def resolve_release():
    """Return (version, {wp_token: gmt_filename}) for the release to build."""
    pinned = os.environ.get("WP_VERSION", "").strip()
    index = f"{BASE_URL}/{pinned}/gmt/" if pinned else f"{BASE_URL}/current/gmt/"
    print(f"resolving release from {index}")

    response = requests.get(index, timeout=120)
    response.raise_for_status()

    found = {}
    for match in GMT_RE.finditer(response.text):
        found[match.group(2)] = (match.group(1), match.group(0))
    if not found:
        sys.exit(f"ERROR: no GMT filenames found at {index} (index format changed?)")

    versions = {version for version, _ in found.values()}
    if len(versions) != 1:
        # Building one linkset from two releases is far worse than not building.
        sys.exit(f"ERROR: {index} mixes releases {sorted(versions)}; refusing to build")
    version = versions.pop()
    if pinned and version != pinned:
        sys.exit(f"ERROR: asked for release {pinned} but the index served {version}")

    print(f"  release {version}, {len(found)} species available")
    return version, {token: name for token, (_, name) in found.items()}


def load_species(path=SPECIES_TSV):
    """Read species.tsv, honouring WP_SPECIES and dropping 'off' rows."""
    only = set(os.environ.get("WP_SPECIES", "").split())
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                sys.exit(f"ERROR: {path}: expected 6 columns, got {len(fields)}: {line!r}")
            code, token, bridgedb_species, syscodes_out, tier, min_frac = fields[:6]
            if tier == "off" or (only and code not in only):
                continue
            rows.append({
                "code": code,
                "wp_token": token,
                "organism": token.replace("_", " "),
                "bridgedb_species": bridgedb_species,
                "syscodes_out": syscodes_out,
                "tier": tier,
                "min_mapped_frac": min_frac,
            })

    if not rows:
        sys.exit(f"ERROR: no species selected from {path}"
                 + (f" (WP_SPECIES={' '.join(sorted(only))})" if only else ""))
    unknown = only - {row["code"] for row in rows}
    if unknown:
        sys.exit(f"ERROR: WP_SPECIES names species not enabled in {path}: "
                 f"{' '.join(sorted(unknown))}")
    return rows


def parse_gmt(text):
    """Parse GMT text into [(pathway_name, pathway_id, version, species, genes)].

    A GMT line is  <set name>\\t<url>\\t<gene>\\t<gene>...
    and the set name is  <name>%WikiPathways_<date>%<WP id>%<species>
    """
    pathways = []
    duplicate_genes = 0
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        # Some set names carry stray surrounding whitespace (e.g. mouse WP3673
        # is written as " Hfe effect on hepcidin production%..."), which would
        # otherwise end up in the pathway node's name and label.
        parts = [part.strip() for part in fields[0].split("%")]
        if len(parts) < 4:
            print(f"  WARNING: unparseable GMT set name, skipped: {fields[0]!r}")
            continue

        genes, seen = [], set()
        for gene in fields[2:]:
            gene = gene.strip()
            if not gene:
                continue
            if gene in seen:
                # A repeated gene would become a duplicate edge; linkset-creator
                # only de-duplicates when source and target ids are equal.
                duplicate_genes += 1
                continue
            seen.add(gene)
            genes.append(gene)
        if genes:
            pathways.append((parts[0], parts[2], parts[1], parts[3], genes))

    if duplicate_genes:
        print(f"  dropped {duplicate_genes} gene(s) repeated within a pathway")
    return pathways


def write_input(code, pathways):
    """Write input_<code>.txt. Returns the number of pathway-gene rows."""
    out = f"input_{code}.txt"
    rows = 0
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(INPUT_COLUMNS)
        for name, pathway_id, version, species, genes in pathways:
            for gene in genes:
                # The gene id doubles as the label: the GMT carries no symbol.
                writer.writerow([name, pathway_id, len(genes), gene, gene, version, species])
                rows += 1
    return rows


def write_config(species, version, template):
    """Generate wikipathways/wp_<code>.config from the template."""
    text = template.format(
        organism=species["organism"],
        version=version,
        code=species["code"],
        syscodes_out=species["syscodes_out"],
    )
    # ConfigFileReader splits on the first '=' but keeps only lines splitting
    # into exactly two parts, silently discarding anything else as an invalid
    # attribute. Catch that here rather than in a linkset missing a field.
    for line in text.splitlines():
        if line.strip() and line.count("=") != 1:
            sys.exit(f"ERROR: generated config line is not a single key=value pair: {line!r}")

    out = pathlib.Path(f"wikipathways/wp_{species['code']}.config")
    out.write_text(text, encoding="utf-8")
    return out


def main():
    version, available = resolve_release()
    species_list = load_species()
    template = pathlib.Path(CONFIG_TEMPLATE).read_text(encoding="utf-8")

    built = []
    for species in species_list:
        code, token, tier = species["code"], species["wp_token"], species["tier"]
        label = f"{code} ({species['organism']})"

        if token not in available:
            message = f"{label}: no GMT for {token} in release {version}"
            if tier == "core":
                sys.exit(f"ERROR: {message}")
            print(f"SKIP {message}")
            continue

        url = f"{BASE_URL}/{version}/gmt/{available[token]}"
        print(f"{label}: {url}")
        response = requests.get(url, timeout=300)
        response.raise_for_status()

        pathways = parse_gmt(response.text)
        if len(pathways) < MIN_PATHWAYS:
            message = (f"{label}: only {len(pathways)} pathway(s), "
                       f"below the minimum of {MIN_PATHWAYS}")
            if tier == "core":
                sys.exit(f"ERROR: {message}")
            print(f"SKIP {message}")
            continue

        n_rows = write_input(code, pathways)
        config = write_config(species, version, template)
        print(f"  {len(pathways)} pathways, {n_rows} pairs -> input_{code}.txt, {config}")

        species["n_pathways"] = len(pathways)
        species["n_rows"] = n_rows
        built.append(species)

    if not built:
        sys.exit("ERROR: no linkset inputs were produced")

    build = pathlib.Path(BUILD_DIR)
    build.mkdir(exist_ok=True)
    (build / "version.txt").write_text(version + "\n", encoding="utf-8")
    with open(build / "manifest.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(MANIFEST_COLUMNS)
        for species in built:
            writer.writerow([species[column] for column in MANIFEST_COLUMNS])

    print(f"\nWikiPathways {version}: built inputs for {len(built)} species "
          f"({' '.join(species['code'] for species in built)})")


if __name__ == "__main__":
    main()
