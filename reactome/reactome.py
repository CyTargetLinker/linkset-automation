#!/usr/bin/env python3
"""Preprocess a Reactome release into linkset-creator input files.

For every species listed in reactome/species.tsv this writes:

    input_<code>.txt                            tab-delimited pathway-gene pairs
    reactome/config/reactome_<code>.config      generated from the template

plus build/version.txt and build/manifest.tsv describing the build, exactly as
wikipathways/wp.py does; the manifest, not species.tsv, is what the workflow
iterates over.

Why NCBI2Reactome.txt and not ReactomePathways.gmt
--------------------------------------------------
The GMT looks like the obvious analogue of the WikiPathways input and is not.
It covers human only (every pathway id is R-HSA), the bytes served at the .gmt
URL are actually a ZIP archive, and the plain .gmt path 403s in the versioned
archive, so it cannot support a pinned, reproducible build. NCBI2Reactome.txt
carries Entrez ids for all species, is a plain TSV in both current and archived
releases, and lands directly on the pipeline's existing syscode L.

    gene_id  pathway_id  url  pathway_name  evidence_code  species

Three properties of that file drive the parsing below, all measured on release
097: a (gene, pathway) pair may appear twice, once TAS and once IEA (4076 human
pairs, an 8% edge inflation if written straight out); a few gene ids are not
Entrez ids at all but GenBank/RefSeq accessions (82 human rows); and some
pathway names carry stray surrounding whitespace (23 human names).

    python3 reactome/reactome.py                    # all enabled species, current release
    python3 reactome/reactome.py --species hsa      # just this one
    python3 reactome/reactome.py --version 97       # a specific release
"""

import argparse
import csv
import gzip
import os
import pathlib
import re
import sys
import time

import requests

VERSION_URL = "https://reactome.org/ContentService/data/database/version"
CURRENT_URL = "https://reactome.org/download/current"
ARCHIVE_URL = "https://download.reactome.org"
DATA_FILE = "NCBI2Reactome.txt"

SPECIES_TSV = "reactome/species.tsv"
CONFIG_TEMPLATE = "reactome/reactome.config.template"
CONFIG_DIR = "reactome/config"
BUILD_DIR = "build"

# NCBI2Reactome.txt column order.
COL_GENE, COL_PATHWAY, COL_NAME, COL_EVIDENCE, COL_SPECIES = 0, 1, 3, 4, 5
N_COLUMNS = 6

# An Entrez GeneID is all digits. Anything else in that column is an accession
# Reactome carries for a non-gene entity and BridgeDb cannot map.
ENTREZ_RE = re.compile(r"^\d+$")

NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO"

NCBI_CLADES = frozenset({
    "Archaea_Bacteria", "Fungi", "Invertebrates", "Mammalia",
    "Non-mammalian_vertebrates", "Plants", "Protozoa", "Viruses",
})

PLACEHOLDER_SYMBOLS = frozenset({"-", "NEWENTRY"})

# A value carrying any of these would corrupt the input file: linkset-creator
# splits rows with Java String.split("\t") rather than a CSV parser, so an
# embedded tab shifts every later column and a newline truncates the linkset.
UNSAFE_IN_FIELD = frozenset("\t\r\n\"")

# A pathway with fewer genes than this is not worth a node.
MIN_PATHWAYS = 5

DOWNLOAD_ATTEMPTS = 3

MANIFEST_COLUMNS = [
    "code", "ncbi_token", "organism", "bridgedb_species",
    "syscodes_out", "min_mapped_frac", "n_pathways", "n_rows",
    # Append here, never insert: the workflow reads this file positionally.
    "n_genes", "reject_frac",
]

INPUT_COLUMNS = [
    "pathway_name", "pathway_id", "gene_count",
    "gene_id", "gene_symbol", "version", "species", "evidence",
]


def probe(url):
    """True if url exists. Reactome answers an absent release with 403, not 404."""
    try:
        response = requests.get(url, headers={"Range": "bytes=0-0"}, timeout=60)
    except requests.RequestException as error:
        print(f"  cannot reach {url} ({error})")
        return False
    return response.status_code in (200, 206)


def resolve_release(pinned=None):
    """Return (version, base_url) for the release to build.

    The version is Reactome's release integer, zero-padded to three digits. The
    padding is not cosmetic: the Zenodo gate compares version strings with the
    shell's `\\<`, under which "100" sorts before "97" and publishing would
    silently stop at release 100.

    The archive is preferred over /current/ because it is pinnable, which is
    what makes a build reproducible; /current/ is only a fallback for a release
    that ContentService already announces but the archive has not yet published.
    """
    if pinned:
        try:
            release = int(pinned)
        except ValueError:
            sys.exit(f"ERROR: --version must be a release number, got {pinned!r}")
    else:
        try:
            response = requests.get(VERSION_URL, timeout=60)
            response.raise_for_status()
            release = int(response.text.strip())
        except (requests.RequestException, ValueError) as error:
            sys.exit(f"ERROR: cannot read the Reactome release from {VERSION_URL}: {error}")

    version = f"{release:03d}"
    archive = f"{ARCHIVE_URL}/{release}"
    print(f"resolving release {version}")

    if probe(f"{archive}/{DATA_FILE}"):
        print(f"  using the versioned archive {archive}")
        return version, archive
    if pinned:
        sys.exit(f"ERROR: release {release} is not in the archive at {archive}")

    print(f"  WARNING: release {release} is not in the archive yet, "
          f"falling back to {CURRENT_URL}")
    return version, CURRENT_URL


def load_species(selected=None, path=SPECIES_TSV):
    """Read species.tsv, dropping 'off' rows."""
    only = set(selected or ())
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                sys.exit(f"ERROR: {path}: expected 9 columns, got {len(fields)}: {line!r}")
            (code, reactome_species, ncbi_token, clade, bridgedb_species,
             syscodes_out, tier, min_frac, max_reject) = fields[:9]
            if tier == "off" or (only and code not in only):
                continue
            if clade not in NCBI_CLADES:
                sys.exit(f"ERROR: {path}: {code} has unknown ncbi_clade {clade!r}; "
                         f"expected one of {', '.join(sorted(NCBI_CLADES))}")
            rows.append({
                "code": code,
                "reactome_species": reactome_species,
                "ncbi_token": ncbi_token,
                "ncbi_clade": clade,
                "organism": reactome_species,
                "bridgedb_species": bridgedb_species,
                "syscodes_out": syscodes_out,
                "tier": tier,
                "min_mapped_frac": min_frac,
                "max_reject_frac": float(max_reject),
            })

    if not rows:
        sys.exit(f"ERROR: no species selected from {path}"
                 + (f" (--species {' '.join(sorted(only))})" if only else ""))
    unknown = only - {row["code"] for row in rows}
    if unknown:
        sys.exit(f"ERROR: --species names species not enabled in {path}: "
                 f"{' '.join(sorted(unknown))}")
    return rows


def download(url, local, attempts=DOWNLOAD_ATTEMPTS):
    """Stream url to local, reusing an existing non-empty copy."""
    if os.path.exists(local) and os.path.getsize(local) > 0:
        print(f"  using existing {local} (delete it to refresh)")
        return local

    tmp = f"{local[:-3]}.partial.gz" if local.endswith(".gz") else f"{local}.partial"
    delay = 2
    try:
        for attempt in range(1, attempts + 1):
            try:
                with requests.get(url, stream=True, timeout=600) as response:
                    if response.status_code in (403, 404):
                        # Neither is transient. Reactome returns 403 for a
                        # release that does not exist, and NCBI 404s a species
                        # it publishes no per-species gene_info for.
                        sys.exit(f"ERROR: {url} is not available "
                                 f"(HTTP {response.status_code})")
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


def usable_symbol(gene_id, symbol):
    """True if gene_info's Symbol is a name rather than an NCBI placeholder."""
    if not symbol or symbol in PLACEHOLDER_SYMBOLS:
        return False
    if symbol == f"LOC{gene_id}":
        return False
    return not (UNSAFE_IN_FIELD & set(symbol))


def load_symbols(species):
    """Return {GeneID: Symbol} from this species' NCBI gene_info file.

    Deliberately does NOT filter on tax_id: a per-species file also carries
    strain and subspecies rows, and GeneIDs are globally unique anyway.
    """
    url = (f"{NCBI_GENE_INFO_URL}/{species['ncbi_clade']}/"
           f"{species['ncbi_token']}.gene_info.gz")
    local = f"gene_info_{species['code']}.gz"
    print(f"  symbols: {url}")
    download(url, local)

    symbols = {}
    try:
        with gzip.open(local, "rt", encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if not row or row[0].startswith("#") or len(row) < 3:
                    continue
                gene_id, symbol = row[1], row[2]
                if usable_symbol(gene_id, symbol):
                    symbols[gene_id] = symbol
    except (OSError, EOFError, UnicodeDecodeError) as error:
        sys.exit(f"ERROR: cannot read {local}: {error} (delete it and re-run)")

    if not symbols:
        sys.exit(f"ERROR: {local} contains no usable gene symbols")
    return symbols


def evidence_label(codes):
    """Collapse one pair's evidence codes into a single field.

    Reactome writes a pair supported both by curation and by projection as two
    rows, TAS and IEA. Written straight out those become two identical edges, so
    the codes are aggregated here and the pair emitted once. Curated first, so
    the strongest evidence leads.
    """
    ordered = [code for code in ("TAS", "IEA") if code in codes]
    ordered += sorted(codes - {"TAS", "IEA"})
    # Never empty: an empty final column makes linkset-creator's split("\t")
    # return a short row, which truncates the linkset silently.
    return "+".join(ordered) or "NA"


def parse_ncbi2reactome(path, species_name):
    """Return ({pathway_id: {name, genes: {gene_id: {codes}}}}, stats)."""
    pathways = {}
    stats = {"rows": 0, "rejected_id": 0, "rejected_unsafe": 0, "duplicate_pairs": 0}

    with open(path, encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < N_COLUMNS or row[COL_SPECIES] != species_name:
                continue
            stats["rows"] += 1

            gene_id = row[COL_GENE].strip()
            if not ENTREZ_RE.match(gene_id):
                # A GenBank/RefSeq accession rather than a gene, e.g. NC_012920
                # or MN908947.3. BridgeDb maps none of them.
                stats["rejected_id"] += 1
                continue

            pathway_id = row[COL_PATHWAY].strip()
            name = row[COL_NAME].strip()
            if UNSAFE_IN_FIELD & set(name) or UNSAFE_IN_FIELD & set(pathway_id):
                stats["rejected_unsafe"] += 1
                continue

            entry = pathways.setdefault(pathway_id,
                                        {"name": name or pathway_id, "genes": {}})
            codes = entry["genes"].setdefault(gene_id, set())
            if codes:
                stats["duplicate_pairs"] += 1
            code = row[COL_EVIDENCE].strip()
            if code:
                codes.add(code)

    stats["pathways"] = len(pathways)
    return pathways, stats


def apply_symbols(pathways, symbols):
    """Drop genes NCBI has not named, and pathways left empty. Returns stats.

    For WikiPathways an unnamed gene keeps its Entrez id as a label. Here it is
    rejected outright, because the ids gene_info does not carry are pathogen-side
    genes filed under the host species in Reactome's infectious-disease pathways
    (Escherichia ruysiae skp, secY, secE) plus non-Entrez accessions. None of
    them map in BridgeDb, so keeping them would only add identifier-less nodes
    unreachable from any gene query. It is counted and gated rather than silent,
    because it does remove whole pathways: on release 097, 33 human pathways
    disappear entirely and 69 lose at least half their genes.
    """
    genes_before = {gene for entry in pathways.values() for gene in entry["genes"]}
    emptied = []

    for pathway_id, entry in list(pathways.items()):
        entry["genes"] = {gene: codes for gene, codes in entry["genes"].items()
                          if gene in symbols}
        if not entry["genes"]:
            emptied.append(entry["name"])
            del pathways[pathway_id]

    genes_after = {gene for entry in pathways.values() for gene in entry["genes"]}
    rejected = len(genes_before) - len(genes_after)
    return {
        "genes_before": len(genes_before),
        "genes_after": len(genes_after),
        "rejected_genes": rejected,
        "reject_frac": rejected / len(genes_before) if genes_before else 0.0,
        "emptied_pathways": emptied,
    }


def write_input(code, pathways, symbols, version, species_name):
    """Write input_<code>.txt. Returns the number of pathway-gene rows."""
    out = f"input_{code}.txt"
    rows = 0
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(INPUT_COLUMNS)
        for pathway_id, entry in pathways.items():
            for gene_id, codes in entry["genes"].items():
                writer.writerow([
                    entry["name"], pathway_id, len(entry["genes"]),
                    gene_id, symbols[gene_id], version, species_name,
                    evidence_label(codes),
                ])
                rows += 1
    return rows


def write_config(species, version, template):
    """Generate reactome/config/reactome_<code>.config from the template."""
    text = template.format(
        organism=species["organism"],
        version=version,
        code=species["code"],
        syscodes_out=species["syscodes_out"],
    )
    # ConfigFileReader keeps only lines splitting into exactly two parts on '=',
    # silently discarding the rest as an invalid attribute.
    for line in text.splitlines():
        if line.strip() and line.count("=") != 1:
            sys.exit(f"ERROR: generated config line is not a single key=value pair: {line!r}")

    out = pathlib.Path(CONFIG_DIR) / f"reactome_{species['code']}.config"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build Reactome linkset inputs and configs.")
    parser.add_argument(
        "--version", metavar="N",
        help="Reactome release number to build (default: the current release)")
    parser.add_argument(
        "--species", nargs="+", metavar="CODE",
        help="species codes to build, e.g. hsa (default: all enabled)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    species_list = load_species(args.species)
    version, base_url = resolve_release(args.version)
    template = pathlib.Path(CONFIG_TEMPLATE).read_text(encoding="utf-8")

    url = f"{base_url}/{DATA_FILE}"
    print(f"pathway-gene table: {url}")
    download(url, DATA_FILE)

    built = []
    for species in species_list:
        code, tier = species["code"], species["tier"]
        label = f"{code} ({species['organism']})"
        print(f"{label}:")

        pathways, stats = parse_ncbi2reactome(DATA_FILE, species["reactome_species"])
        print(f"  {stats['rows']} rows -> {stats['pathways']} pathways; "
              f"collapsed {stats['duplicate_pairs']} duplicate pair(s), "
              f"rejected {stats['rejected_id']} non-Entrez id(s)"
              + (f", {stats['rejected_unsafe']} unsafe name(s)"
                 if stats["rejected_unsafe"] else ""))

        if len(pathways) < MIN_PATHWAYS:
            message = (f"{label}: only {len(pathways)} pathway(s), "
                       f"below the minimum of {MIN_PATHWAYS}")
            if tier == "core":
                sys.exit(f"ERROR: {message}")
            print(f"SKIP {message}")
            continue

        symbols = load_symbols(species)
        joined = apply_symbols(pathways, symbols)
        print(f"  {len(symbols)} symbols in gene_info; rejected "
              f"{joined['rejected_genes']}/{joined['genes_before']} "
              f"({joined['reject_frac']:.2%}) gene(s) NCBI has not named, "
              f"emptying {len(joined['emptied_pathways'])} pathway(s)")
        if joined["reject_frac"] > species["max_reject_frac"]:
            sys.exit(f"ERROR: {label}: rejected {joined['reject_frac']:.2%} of genes, "
                     f"above the ceiling of {species['max_reject_frac']:.0%}")

        n_rows = write_input(code, pathways, symbols, version,
                             species["reactome_species"])
        config = write_config(species, version, template)
        print(f"  {len(pathways)} pathways, {n_rows} pairs -> input_{code}.txt, {config}")

        species["n_pathways"] = len(pathways)
        species["n_rows"] = n_rows
        species["n_genes"] = joined["genes_after"]
        species["reject_frac"] = f"{joined['reject_frac']:.4f}"
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

    print(f"\nReactome {version}: built inputs for {len(built)} species "
          f"({' '.join(species['code'] for species in built)})")


if __name__ == "__main__":
    main()
