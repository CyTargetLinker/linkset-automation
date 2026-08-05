#!/usr/bin/env python3
"""Preprocess WikiPathways GMT releases into linkset-creator input files.

For every species listed in wikipathways/species.tsv this writes:

    input_<code>.txt                        tab-delimited pathway-gene pairs
    wikipathways/config/wp_<code>.config    generated from wp.config.template

Gene nodes are labelled with their NCBI gene symbol, looked up per species in
NCBI gene_info; the GMT itself carries only Entrez ids. A gene NCBI has not
named keeps its Entrez id as the label.

plus two files describing the build:

    build/version.txt                the WikiPathways release that was used
    build/manifest.tsv               one row per species ACTUALLY built

The manifest, not species.tsv, is what the GitHub workflow iterates over, so a
species skipped here can never leave the workflow looking for a config, a
BridgeDb download or an XGMML file that was never produced.

Release resolution
------------------
data.wikipathways.org keeps only a rolling 12-month window of monthly releases,
so the release is discovered rather than pinned in source.

Resolution deliberately walks the dated release directories listed on the site
root (20260710/, 20260610/, ...) newest first, rather than trusting /current/.
/current/ is a moving alias and says nothing about whether that release is
complete; a dated directory that contains a GMT for every core species is
positive evidence that the release actually shipped. An incomplete newest
directory (a release still being published) is skipped with a warning and the
previous one is used, which is safe because the workflow refuses to publish a
release older than the one already on Zenodo.

Files are then fetched from the resolved /<version>/ path, so a release
appearing mid-run cannot yield a linkset built from two releases.

Both arguments are optional and are wired to the workflow's manual-dispatch
inputs:

    python3 wikipathways/wp.py                              # all species, current release
    python3 wikipathways/wp.py --species sce aga            # just these two
    python3 wikipathways/wp.py --version 20260510           # a specific release
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

BASE_URL = "https://data.wikipathways.org"
SPECIES_TSV = "wikipathways/species.tsv"
CONFIG_TEMPLATE = "wikipathways/wp.config.template"
# Generated per-species configs live together in their own directory, so the
# whole directory can be gitignored and cleared without touching sources.
CONFIG_DIR = "wikipathways/config"
BUILD_DIR = "build"

# wikipathways-<YYYYMMDD>-gmt-<Genus_species>.gmt
GMT_RE = re.compile(r"wikipathways-(\d{8})-gmt-([A-Za-z_]+)\.gmt")

# Dated release directories linked from the site root, e.g. href="20260710"
RELEASE_RE = re.compile(r'href="(\d{8})/?"')

# How many dated directories to look at before giving up. More than a couple of
# incomplete releases in a row is a problem to report, not to work around.
MAX_RELEASE_CANDIDATES = 3

# An 'extra' species with fewer pathways than this is not worth a linkset and is
# skipped; a 'core' species below it means something is wrong upstream.
MIN_PATHWAYS = 5

# NCBI per-species gene_info gives GeneID -> Symbol. The GMT carries no symbol,
# so without this every gene node is labelled with a bare Entrez number. The
# filename is the WikiPathways token for all 13 species -- including
# Canis_familiaris, which NCBI does NOT file under Canis_lupus_familiaris -- so
# only the clade directory differs, and it is not derivable from the token.
NCBI_GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO"

# Validated when species.tsv is read, so a typo fails before any download.
NCBI_CLADES = frozenset({
    "Archaea_Bacteria", "Fungi", "Invertebrates", "Mammalia",
    "Non-mammalian_vertebrates", "Plants", "Protozoa", "Viruses",
})

# What NCBI writes for a gene it has not named. LOC<GeneID> is handled
# separately because it embeds the gene's own id.
PLACEHOLDER_SYMBOLS = frozenset({"-", "NEWENTRY"})

# A symbol carrying any of these would corrupt the input file: csv.writer would
# quote the field, but linkset-creator splits rows with Java String.split("\t")
# rather than a CSV parser, so an embedded tab shifts every later column.
UNSAFE_IN_SYMBOL = frozenset("\t\r\n\"")

# Fail a species below this fraction of its genes carrying a symbol. Measured on
# release 20260710: 100% for mmu/sce/cel/dme/ath, 98.9-99.9% elsewhere, and
# 81.2% for aga (NCBI has named only a third of Anopheles genes). Well below the
# worst real value, but far above the 0% of a gene_info file that downloaded but
# was empty, truncated, or for the wrong organism.
MIN_SYMBOL_FRAC = 0.5

DOWNLOAD_ATTEMPTS = 3

MANIFEST_COLUMNS = [
    "code", "wp_token", "organism", "bridgedb_species",
    "syscodes_out", "min_mapped_frac", "n_pathways", "n_rows",
    # Append here, never insert: the workflow reads this file positionally
    # (`cut -f3`, and `read -r code token organism bridgedb_species rest`).
    "n_genes", "symbol_frac",
]

INPUT_COLUMNS = [
    "pathway_name", "pathway_id", "gene_count",
    "gene_id", "gene_symbol", "version", "species",
]


def list_releases():
    """Return the dated release directories on the site root, newest first."""
    response = requests.get(f"{BASE_URL}/", timeout=120)
    response.raise_for_status()
    releases = sorted({match.group(1) for match in RELEASE_RE.finditer(response.text)},
                      reverse=True)
    if not releases:
        sys.exit(f"ERROR: no dated release directories found at {BASE_URL}/ "
                 "(index format changed?)")
    return releases


def read_release(version):
    """Return {wp_token: gmt_filename} for one dated release, or {} if unusable."""
    index = f"{BASE_URL}/{version}/gmt/"
    try:
        response = requests.get(index, timeout=120)
        response.raise_for_status()
    except requests.RequestException as error:
        print(f"  {version}: cannot read {index} ({error})")
        return {}

    found = {}
    for match in GMT_RE.finditer(response.text):
        found[match.group(2)] = (match.group(1), match.group(0))
    if not found:
        print(f"  {version}: no GMT files listed")
        return {}

    stamped = {stamp for stamp, _ in found.values()}
    if stamped != {version}:
        # A dated directory holding differently-stamped files is not a release
        # we can reason about; building from it could mix two releases.
        print(f"  {version}: contains files stamped {sorted(stamped)}")
        return {}

    return {token: name for token, (_, name) in found.items()}


def resolve_release(required_tokens, pinned=None):
    """Return (version, {wp_token: gmt_filename}) for the newest usable release.

    A release counts as published only once it carries a GMT for every core
    species; /current/ is not consulted because it is an alias that says
    nothing about whether the release behind it is complete.
    """
    candidates = [pinned] if pinned else list_releases()[:MAX_RELEASE_CANDIDATES]
    print(f"resolving release from {BASE_URL}/ "
          f"(candidates: {', '.join(candidates)})")

    for version in candidates:
        available = read_release(version)
        if not available:
            continue
        missing = sorted(required_tokens - available.keys())
        if missing:
            print(f"  {version}: incomplete, missing {', '.join(missing)}")
            continue
        if not version.endswith("10"):
            # WikiPathways releases on the 10th; a different day is worth
            # noticing but is not a reason to refuse to build.
            print(f"  WARNING: release {version} is not dated the 10th")
        print(f"  using release {version} ({len(available)} species available)")
        return version, available

    if pinned:
        sys.exit(f"ERROR: release {pinned} is missing or incomplete")
    sys.exit(f"ERROR: none of the {len(candidates)} newest releases "
             f"({', '.join(candidates)}) is complete")


def load_species(selected=None, path=SPECIES_TSV):
    """Read species.tsv, dropping 'off' rows.

    `selected` narrows the result to those codes; None returns the full enabled
    set, which the caller uses to ask what a release must contain regardless of
    which subset this run was asked to build.
    """
    only = set(selected or ())
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                sys.exit(f"ERROR: {path}: expected 7 columns, got {len(fields)}: {line!r}")
            (code, token, clade, bridgedb_species,
             syscodes_out, tier, min_frac) = fields[:7]
            if tier == "off" or (only and code not in only):
                continue
            if clade not in NCBI_CLADES:
                sys.exit(f"ERROR: {path}: {code} has unknown ncbi_clade {clade!r}; "
                         f"expected one of {', '.join(sorted(NCBI_CLADES))}")
            rows.append({
                "code": code,
                "wp_token": token,
                "ncbi_clade": clade,
                "organism": token.replace("_", " "),
                "bridgedb_species": bridgedb_species,
                "syscodes_out": syscodes_out,
                "tier": tier,
                "min_mapped_frac": min_frac,
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
    """Stream url to local, reusing an existing non-empty copy.

    Written to a sibling temp file and renamed only once the transfer finished,
    so an interrupted run cannot leave a truncated file that the next run
    happily reuses. The temp name keeps the .gz suffix so .gitignore covers it.
    """
    if os.path.exists(local) and os.path.getsize(local) > 0:
        print(f"  using existing {local} (delete it to refresh)")
        return local

    tmp = f"{local[:-3]}.partial.gz" if local.endswith(".gz") else f"{local}.partial"
    delay = 2
    try:
        for attempt in range(1, attempts + 1):
            try:
                with requests.get(url, stream=True, timeout=600) as response:
                    if response.status_code == 404:
                        # Not transient, so do not spend retries on it: either
                        # the ncbi_clade is wrong, or NCBI publishes no
                        # per-species file for this organism (it does so for a
                        # shortlist only -- horse, tomato and poplar have none).
                        sys.exit(f"ERROR: {url} does not exist (404). Check the "
                                 f"ncbi_clade column against the listing at "
                                 f"{NCBI_GENE_INFO_URL}/")
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
    """True if gene_info's Symbol is a name rather than an NCBI placeholder.

    NCBI never leaves Symbol empty: an unnamed gene is given LOC<GeneID>. The
    test is equality with this row's own GeneID rather than a "LOC" prefix, so a
    real gene whose name merely starts with LOC is never discarded.

    Systematic-looking names that are NOT placeholders are deliberately kept:
    Arabidopsis AGI codes (AT2G01050) and yeast systematic names (YCL068C) are
    the identifiers those communities actually use.
    """
    if not symbol or symbol in PLACEHOLDER_SYMBOLS:
        return False
    if symbol == f"LOC{gene_id}":
        return False
    return not (UNSAFE_IN_SYMBOL & set(symbol))


def load_symbols(species):
    """Return {GeneID: Symbol} from this species' NCBI gene_info file.

    Deliberately does NOT filter on tax_id. A per-species file also carries
    strain and subspecies rows: the yeast file holds 6478 rows of 559292 against
    36 of 4932, and the mouse file spans five taxa, so filtering on the nominal
    id would discard nearly everything. Entrez GeneIDs are globally unique, so a
    row for another taxon can never be matched by the wrong species.

    gene_info columns: 0 tax_id  1 GeneID  2 Symbol  3 LocusTag  4 Synonyms ...
    """
    url = f"{NCBI_GENE_INFO_URL}/{species['ncbi_clade']}/{species['wp_token']}.gene_info.gz"
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
        # gzip.BadGzipFile is an OSError; a truncated file raises EOFError.
        sys.exit(f"ERROR: cannot read {local}: {error} (delete it and re-run)")

    if not symbols:
        sys.exit(f"ERROR: {local} contains no usable gene symbols")
    return symbols


def symbol_coverage(pathways, symbols):
    """Return (distinct genes, how many of them have a symbol)."""
    genes = {gene for pathway in pathways for gene in pathway[4]}
    return len(genes), sum(1 for gene in genes if gene in symbols)


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


def write_input(code, pathways, symbols):
    """Write input_<code>.txt. Returns the number of pathway-gene rows."""
    out = f"input_{code}.txt"
    rows = 0
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(INPUT_COLUMNS)
        for name, pathway_id, version, species, genes in pathways:
            for gene in genes:
                # The GMT has no symbol; it comes from NCBI gene_info, and a
                # gene NCBI has not named keeps its Entrez id as the label.
                writer.writerow([name, pathway_id, len(genes), gene,
                                 symbols.get(gene, gene), version, species])
                rows += 1
    return rows


def write_config(species, version, template):
    """Generate wikipathways/config/wp_<code>.config from the template."""
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

    out = pathlib.Path(CONFIG_DIR) / f"wp_{species['code']}.config"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build WikiPathways linkset inputs and configs.")
    parser.add_argument(
        "--version", metavar="YYYYMMDD",
        help="WikiPathways release to build (default: the newest complete one)")
    parser.add_argument(
        "--species", nargs="+", metavar="CODE",
        help="species codes to build, e.g. hsa sce (default: all enabled)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    species_list = load_species(args.species)
    # A release counts as published once every core species is in it. Read that
    # from the whole table rather than the --species subset: whether a release
    # shipped does not depend on which species this run was asked to build.
    core_tokens = {species["wp_token"]
                   for species in load_species()
                   if species["tier"] == "core"}
    version, available = resolve_release(core_tokens, pinned=args.version)
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

        symbols = load_symbols(species)
        n_genes, n_symbols = symbol_coverage(pathways, symbols)
        frac = n_symbols / n_genes
        print(f"  {len(symbols)} symbols in gene_info; "
              f"{n_symbols}/{n_genes} ({frac:.1%}) of this species' genes named")
        if frac < MIN_SYMBOL_FRAC:
            sys.exit(f"ERROR: {label}: only {n_symbols}/{n_genes} ({frac:.1%}) of genes "
                     f"have a symbol, below the minimum of {MIN_SYMBOL_FRAC:.0%}")

        n_rows = write_input(code, pathways, symbols)
        config = write_config(species, version, template)
        print(f"  {len(pathways)} pathways, {n_rows} pairs -> input_{code}.txt, {config}")

        species["n_pathways"] = len(pathways)
        species["n_rows"] = n_rows
        species["n_genes"] = n_genes
        species["symbol_frac"] = f"{frac:.4f}"
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
