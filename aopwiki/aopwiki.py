#!/usr/bin/env python3
"""Preprocess AOP-Wiki RDF into linkset-creator input files.

Writes two linkset inputs, each with a config generated into aopwiki/config/:

    input_aopwiki_ke_gene.txt           key event -> gene
    input_aopwiki_stressor_mie.txt      chemical stressor -> molecular initiating event

plus build/version.txt. Both come from the AOP-Wiki RDF SPARQL endpoint.

Key event -> gene: the AOP-Wiki itself has no gene table. These pairs are the
HGNC genes the AOP-Wiki RDF build text-mines from each key event's title and
description (edam:data_1025, the union of the NER and regex detections), so they
are uncurated. Each edge records which method found it.

Stressor -> MIE: the AOP-Wiki attaches stressors to whole AOPs, not to key
events. A prototypical stressor is what starts the pathway, so each chemical
stressor is linked to the molecular initiating event(s) of every AOP that lists
it, and the edge records those AOPs. Chemicals are identified by InChIKey only:
linkset-creator cannot open the current BridgeDb metabolite database
(CyTargetLinker/linksetCreator#5), so they are not expanded to cross-references.
CAS and DTXSID travel as node attributes.

Only AOPs that carry a license in the RDF (CC BY-SA 4.0) are used, so nothing
unlicensed is redistributed: key events must sit on at least one licensed AOP.

    python3 aopwiki/aopwiki.py
"""

import argparse
import csv
import pathlib
import re
import sys
import time

import requests

ENDPOINT = "https://aopwiki.rdf.bigcat-bioinformatics.org/sparql"
CONFIG_DIR = "aopwiki/config"
BUILD_DIR = "build"

PREFIXES = """
PREFIX aopo: <http://aopkb.org/aop_ontology#>
PREFIX cheminf: <http://semanticscience.org/resource/>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX edam: <http://edamontology.org/>
PREFIX ncit: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#>
PREFIX pav: <http://purl.org/pav/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX void: <http://rdfs.org/ns/void#>
PREFIX aopwiki: <https://aopwiki.rdf.bigcat-bioinformatics.org/>
"""

Q_VERSION = PREFIXES + """
SELECT ?version WHERE { ?dataset a void:Dataset ; pav:version ?version }
"""

# Key events that sit on at least one licensed AOP, with their text-mined genes.
Q_KE_GENE = PREFIXES + """
SELECT DISTINCT ?ke ?title ?gene ?symbol WHERE {
  ?ke a aopo:KeyEvent ; dc:title ?title ; edam:data_1025 ?gene .
  OPTIONAL { ?gene rdfs:label ?symbol }
  FILTER EXISTS { ?aop aopo:has_key_event ?ke ; dcterms:license ?license }
}
ORDER BY ?ke ?gene ?title ?symbol
"""

Q_KE_GENE_METHOD = PREFIXES + """
SELECT DISTINCT ?ke ?gene ?method WHERE {
  ?ke a aopo:KeyEvent .
  { ?ke aopwiki:geneDetectedByNER ?gene . BIND("NER" AS ?method) }
  UNION
  { ?ke aopwiki:geneDetectedByRegex ?gene . BIND("regex" AS ?method) }
}
ORDER BY ?ke ?gene ?method
"""

Q_STRESSOR_MIE = PREFIXES + """
SELECT DISTINCT ?inchikey ?chemical ?cas ?dtxsid ?stressor ?aop ?mie ?mie_title WHERE {
  ?aop a aopo:AdverseOutcomePathway ; dcterms:license ?license ;
       ncit:C54571 ?s ; aopo:has_molecular_initiating_event ?mie .
  ?mie dc:title ?mie_title .
  ?s aopo:has_chemical_entity ?c ; dc:title ?stressor .
  ?c cheminf:CHEMINF_000059 ?inchikey .
  OPTIONAL { ?c rdfs:label ?chemical }
  OPTIONAL { ?c cheminf:CHEMINF_000446 ?cas }
  OPTIONAL { ?c cheminf:CHEMINF_000568 ?dtxsid }
}
ORDER BY ?inchikey ?mie ?aop ?stressor ?chemical ?cas ?dtxsid ?mie_title
"""

KE_GENE_COLUMNS = ["ke_id", "ke_title", "hgnc_id", "gene_symbol", "detected_by"]
STRESSOR_COLUMNS = [
    "inchikey", "chemical_name", "casrn", "dtxsid",
    "ke_id", "ke_title", "aops", "stressor_names",
]

# (name, minimum edges). A build below the floor means the endpoint or the RDF
# model changed, not that the AOP-Wiki shrank.
FLOORS = {"ke_gene": 2000, "stressor_mie": 200}

INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

# linkset-creator splits rows with Java String.split("\t") rather than a CSV
# parser, so an embedded tab shifts every later column and a newline truncates
# the linkset.
UNSAFE_RE = re.compile(r"[\t\r\n\"]+")

QUERY_ATTEMPTS = 3

# The endpoint silently truncates any result at 10,000 rows, so large queries
# are paged. Paging needs a stable order, hence the ORDER BY on each of them.
PAGE_SIZE = 5000


def sparql_all(query):
    """Run an ORDER BY query page by page and return every row."""
    rows = []
    while True:
        page = sparql(f"{query}\nLIMIT {PAGE_SIZE} OFFSET {len(rows)}")
        rows.extend(page)
        if len(page) < PAGE_SIZE:
            return rows


def sparql(query, attempts=QUERY_ATTEMPTS):
    """Run a SELECT query and return its bindings as a list of {var: value}."""
    delay = 5
    for attempt in range(1, attempts + 1):
        try:
            response = requests.post(
                ENDPOINT, data={"query": query},
                headers={"Accept": "application/sparql-results+json"}, timeout=300)
            response.raise_for_status()
            bindings = response.json()["results"]["bindings"]
            return [{k: v["value"] for k, v in row.items()} for row in bindings]
        except (requests.RequestException, ValueError, KeyError) as error:
            print(f"  attempt {attempt}/{attempts} failed: {error}")
            if attempt == attempts:
                sys.exit(f"ERROR: SPARQL query against {ENDPOINT} failed")
            time.sleep(delay)
            delay *= 2


def local(iri):
    """The last path segment of an identifiers.org IRI: .../hgnc/108 -> 108."""
    return iri.rstrip("/").rsplit("/", 1)[-1]


def ke_id(iri):
    return f"aop.events:{local(iri)}"


def aop_id(iri):
    return f"aop:{local(iri)}"


def clean(value):
    return UNSAFE_RE.sub(" ", value or "").strip()


def resolve_version():
    """Return the AOP-Wiki RDF release as YYYY-MM-DD (from pav:version 2026.10.10)."""
    rows = sparql(Q_VERSION)
    versions = sorted({row["version"] for row in rows if "version" in row})
    if len(versions) != 1:
        sys.exit(f"ERROR: expected one pav:version on the void:Dataset, got {versions}")
    match = re.fullmatch(r"(\d{4})\.(\d{1,2})\.(\d{1,2})", versions[0])
    if not match:
        sys.exit(f"ERROR: unexpected AOP-Wiki RDF version {versions[0]!r}")
    year, month, day = match.groups()
    return f"{year}-{int(month):02d}-{int(day):02d}"


def build_ke_gene():
    methods = {}
    for row in sparql_all(Q_KE_GENE_METHOD):
        methods.setdefault((row["ke"], row["gene"]), set()).add(row["method"])

    edges = {}
    for row in sparql_all(Q_KE_GENE):
        key = (row["ke"], row["gene"])
        if key in edges:
            continue
        hgnc = local(row["gene"])
        if not hgnc.isdigit():
            continue
        edges[key] = [
            ke_id(row["ke"]), clean(row["title"]),
            f"HGNC:{hgnc}", clean(row.get("symbol")) or f"HGNC:{hgnc}",
            ";".join(sorted(methods.get(key, ()))),
        ]
    return sorted(edges.values())


def build_stressor_mie():
    edges = {}
    rejected = set()
    for row in sparql_all(Q_STRESSOR_MIE):
        inchikey = local(row["inchikey"])
        if not INCHIKEY_RE.match(inchikey):
            rejected.add(inchikey)
            continue
        key = (inchikey, row["mie"])
        edge = edges.setdefault(key, {
            "chemical": clean(row.get("chemical")) or inchikey,
            "cas": set(), "dtxsid": set(),
            "mie_title": clean(row["mie_title"]),
            "aops": set(), "stressors": set(),
        })
        if row.get("cas"):
            edge["cas"].add(local(row["cas"]))
        if row.get("dtxsid"):
            edge["dtxsid"].add(local(row["dtxsid"]))
        edge["aops"].add(aop_id(row["aop"]))
        edge["stressors"].add(clean(row["stressor"]))
    if rejected:
        print(f"  skipped {len(rejected)} malformed InChIKey(s): "
              f"{' '.join(sorted(rejected)[:5])}")

    def by_number(values):
        return ";".join(sorted(values, key=lambda v: int(v.split(":")[1])))

    return [
        [inchikey, edge["chemical"], ";".join(sorted(edge["cas"])),
         ";".join(sorted(edge["dtxsid"])), ke_id(mie), edge["mie_title"],
         by_number(edge["aops"]), ";".join(sorted(edge["stressors"]))]
        for (inchikey, mie), edge in sorted(edges.items())
    ]


def write_input(name, columns, rows):
    out = f"input_aopwiki_{name}.txt"
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", quoting=csv.QUOTE_NONE,
                            lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    return out


def write_config(name, version):
    template = pathlib.Path(f"aopwiki/aopwiki_{name}.config.template").read_text(encoding="utf-8")
    text = template.format(version=version)
    # ConfigFileReader keeps only lines splitting into exactly two parts on '=',
    # silently discarding the rest as an invalid attribute.
    for line in text.splitlines():
        if line.strip() and line.count("=") != 1:
            sys.exit(f"ERROR: generated config line is not a single key=value pair: {line!r}")
    out = pathlib.Path(CONFIG_DIR) / f"aopwiki_{name}.config"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return out


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Build AOP-Wiki linkset inputs and configs.")
    parser.add_argument(
        "--version", metavar="YYYY-MM-DD",
        help="version label (default: the AOP-Wiki RDF release the endpoint serves). "
             "The content is always what the endpoint serves now.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    version = args.version or resolve_version()
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", version):
        sys.exit(f"ERROR: version must be YYYY-MM-DD, got {version!r}")
    print(f"AOP-Wiki RDF {version}")

    linksets = [
        ("ke_gene", KE_GENE_COLUMNS, build_ke_gene),
        ("stressor_mie", STRESSOR_COLUMNS, build_stressor_mie),
    ]
    for name, columns, build in linksets:
        rows = build()
        sources = {row[0] for row in rows}
        targets = {row[2] if name == "ke_gene" else row[4] for row in rows}
        print(f"{name}: {len(rows)} edges, {len(sources)} sources, {len(targets)} targets")
        if len(rows) < FLOORS[name]:
            sys.exit(f"ERROR: {name}: only {len(rows)} edges, below the floor of {FLOORS[name]}")
        out = write_input(name, columns, rows)
        config = write_config(name, version)
        print(f"  -> {out}, {config}")

    build = pathlib.Path(BUILD_DIR)
    build.mkdir(exist_ok=True)
    (build / "version.txt").write_text(version + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
