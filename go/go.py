"""GO -> gene preprocessing for CyTargetLinker linksets.

Downloads NCBI gene2go and writes one tab-delimited input_<code>.txt per species
in the column order expected by go_<code>.config:

    go_id  go_term  gene_id  gene_label  evidence  aspect  qualifier  pubmed

Only experimental evidence codes are kept (curated experimental annotations);
negated (NOT-qualified) annotations are dropped. gene2go carries Entrez Gene IDs
(BridgeDb syscode L), so the linkset matches the WikiPathways gene identifiers.

GO-hierarchy propagation (true-path rule)
-----------------------------------------
gene2go records only the most specific (direct) annotation for each gene. By the
GO true-path rule a gene annotated to a term is implicitly annotated to all of
that term's ancestors, so a direct-only linkset leaves high-level terms looking
almost empty (e.g. "response to decreased oxygen levels" had 1 direct human
gene but ~hundreds with descendants). We therefore download go-basic.obo and, for
every direct annotation, also emit an edge from the gene to each ancestor term
reachable via is_a / part_of. (regulates is intentionally excluded: a regulator
of process X is not a participant in X, so propagating gene annotations across it
would be biologically wrong.) Ancestor term names come from the ontology.

Each GO-term -> gene pair is written exactly once. A pair can be reached from
several direct annotations (directly and via different descendants), so the
evidence codes, qualifiers and PubMed IDs of every contributing annotation are
merged ('|'-joined) into one edge. (linksetCreator does not de-duplicate edges
itself because GO IDs and Entrez IDs never collide.)

gene2go has no gene symbol, so the readable label comes from NCBI gene_info
(GeneID -> Symbol). Without it the genes would appear in Cytoscape as bare
Entrez numbers. If a symbol is missing the Entrez ID is used as the label.

gene2go has no release version, so VERSION is just a date stamp for the linkset
name. The ontology does carry one: go-basic.obo's data-version header names the
GO release (e.g. "releases/2026-05-19"). GO is CC BY 4.0 and its citation policy
asks redistributed data to state the release it came from, so check_version()
warns when VERSION has drifted from the ontology actually downloaded.

To add species, append an entry to SPECIES (tax_id, code, gene_info_url) and add
a matching go_<code>.config.
"""

import csv
import gzip
import os

import requests

VERSION = "20260601"
GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
LOCAL = "gene2go.gz"

# go-basic is the filtered, cycle-free ontology meant for annotation propagation.
GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
GO_OBO_LOCAL = "go-basic.obo"

# Relationship types over which a gene annotation propagates to ancestor terms.
PROPAGATE_RELATIONS = {"is_a", "part_of"}

# Drop terms annotated to more than this many genes. After true-path propagation
# the ontology roots and other near-root terms link to most of the genome
# (e.g. "molecular_function" -> ~16k genes), which is biologically uninformative
# and bloats the linkset. This is the standard "too generic" cutoff used in GO
# enrichment. Set to None to keep every term.
MAX_GENES = 500

# (NCBI tax_id, CyTargetLinker short code, per-species gene_info.gz URL) — human only for now.
# gene_info supplies the readable gene Symbol used as the node label.
SPECIES = [
    (
        "9606",
        "hsa",
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    ),
]

# GO experimental evidence codes (core experimental + high-throughput experimental).
# See http://geneontology.org/docs/guide-go-evidence-codes/
EXPERIMENTAL = {
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",   # core experimental
    "HTP", "HDA", "HMP", "HGI", "HEP",          # high-throughput experimental
}

# gene2go columns (tab-separated; header line starts with '#')
# 0 tax_id  1 GeneID  2 GO_ID  3 Evidence  4 Qualifier  5 GO_term  6 PubMed  7 Category

# The three GO aspects (gene2go "Category") are written as separate linksets:
# propagation never crosses namespaces, so each is self-contained. Maps the
# gene2go category to the output file suffix and a matching go_<code>_<suffix>.config.
ASPECTS = {
    "Process": "bp",     # Biological Process
    "Function": "mf",    # Molecular Function
    "Component": "cc",   # Cellular Component
}


def download(url, local):
    """Download url to local, reusing an existing non-empty local copy."""
    if os.path.exists(local) and os.path.getsize(local) > 0:
        print(f"using existing {local}")
        return local
    print(f"downloading {url}")
    with requests.get(url, stream=True, timeout=600) as resp:
        resp.raise_for_status()
        with open(local, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    print(f"  saved {local}")
    return local


def load_symbols(url, tax):
    """Return {GeneID: Symbol} for one species from an NCBI gene_info.gz file.

    gene_info columns: 0 tax_id  1 GeneID  2 Symbol  ... A '-' symbol means none.
    """
    local = f"gene_info_{tax}.gz"
    download(url, local)
    symbols = {}
    with gzip.open(local, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 3:
                continue
            if row[0] != tax:
                continue
            gene_id, symbol = row[1], row[2]
            if symbol and symbol != "-":
                symbols[gene_id] = symbol
    print(f"  loaded {len(symbols)} symbols for tax {tax}")
    return symbols


def parse_obo(path):
    """Parse go-basic.obo into (parents, names, alt_ids, data_version).

    parents      : {GO_id -> set(direct parent GO_ids)} over PROPAGATE_RELATIONS
    names        : {GO_id -> term name}
    alt_ids      : {secondary GO_id -> primary GO_id}
    data_version : the ontology release from the OBO header (e.g.
                   "releases/2026-05-19"), or None if the header lacks one.
    Obsolete terms are skipped.
    """
    parents, names, alt_ids = {}, {}, {}
    data_version = None
    cur, par, alts, name, obsolete = None, set(), [], None, False

    def flush():
        if cur and not obsolete:
            names[cur] = name or cur
            parents[cur] = par
            for a in alts:
                alt_ids[a] = cur

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                flush()
                cur, par, alts, name, obsolete = None, set(), [], None, False
            elif line.startswith("[") and line.endswith("]"):
                # a non-Term stanza (e.g. [Typedef]) ends term parsing
                flush()
                cur, par, alts, name, obsolete = None, set(), [], None, True
            elif line.startswith("id: GO:"):
                cur = line[4:].strip()
            elif line.startswith("name: "):
                name = line[6:].strip()
            elif line.startswith("alt_id: GO:"):
                alts.append(line[8:].strip())
            elif line.startswith("is_a: GO:") and "is_a" in PROPAGATE_RELATIONS:
                par.add(line[6:].split("!", 1)[0].strip())
            elif line.startswith("relationship: "):
                rel = line[len("relationship: "):].split("!", 1)[0].split()
                if len(rel) >= 2 and rel[0] in PROPAGATE_RELATIONS and rel[1].startswith("GO:"):
                    par.add(rel[1])
            elif line == "is_obsolete: true":
                obsolete = True
            elif cur is None and line.startswith("data-version: "):
                # header stanza, before the first [Term]
                data_version = line[len("data-version: "):].strip()
    flush()
    print(f"  parsed {len(names)} GO terms ({len(alt_ids)} alt ids)")
    print(f"  GO release: {data_version or 'unknown'}")
    return parents, names, alt_ids, data_version


def build_ancestor_fn(parents):
    """Return anc(term) -> set of all ancestor terms (memoized, iterative)."""
    cache = {}

    def anc(term):
        if term in cache:
            return cache[term]
        seen = set()
        stack = list(parents.get(term, ()))
        while stack:
            p = stack.pop()
            if p in seen:
                continue
            seen.add(p)
            stack.extend(parents.get(p, ()))
        cache[term] = seen
        return seen

    return anc


def _merge(existing, value):
    """Append a '|'-delimited token (or tokens) to a set-like ordered list."""
    for tok in value.split("|"):
        if tok and tok not in existing:
            existing.append(tok)


def check_version(data_version):
    """Warn when VERSION does not name the GO release actually downloaded.

    GO is CC BY 4.0 and its citation policy asks that redistributed data state
    the release it came from. VERSION is what the go_<code>_<aspect>.config
    files record as the linkset version, so if it drifts from the ontology on
    disk the deposit cites a release it was not built from. The configs are
    static files, so correcting this is a manual edit — this only makes the
    drift visible instead of silent.
    """
    if not data_version:
        print("WARNING: go-basic.obo has no data-version header; "
              f"cannot verify VERSION={VERSION}")
        return
    # data-version looks like "releases/2026-05-19"; VERSION like "20260601"
    release = data_version.rsplit("/", 1)[-1].replace("-", "")
    if release != VERSION:
        print(f"WARNING: VERSION={VERSION} does not match the GO release "
              f"actually downloaded ({data_version}). Update VERSION and the "
              f"version= line in go/go_<code>_<aspect>.config to {release} "
              "so the linkset cites the release it was built from.")


def main():
    download(GENE2GO_URL, LOCAL)
    download(GO_OBO_URL, GO_OBO_LOCAL)
    parents, names, alt_ids, data_version = parse_obo(GO_OBO_LOCAL)
    check_version(data_version)
    ancestors = build_ancestor_fn(parents)

    taxmap = {tax: code for tax, code, _ in SPECIES}
    symbol_maps = {code: load_symbols(url, tax) for tax, code, url in SPECIES}

    # Per species: (go_id, gene_id) -> aggregated annotation.
    pairs = {code: {} for _, code, _ in SPECIES}

    def record(table, go_id, gene_id, term_name, aspect, evidence, qualifier, pubmed):
        key = (go_id, gene_id)
        entry = table.get(key)
        if entry is None:
            entry = {"go_term": term_name, "aspect": aspect,
                     "evidence": [], "qualifier": [], "pubmed": []}
            table[key] = entry
        _merge(entry["evidence"], evidence)
        if qualifier:
            _merge(entry["qualifier"], qualifier)
        if pubmed and pubmed != "-":
            _merge(entry["pubmed"], pubmed)

    # Stream the gzip line-by-line; the full file is ~10 GB uncompressed.
    with gzip.open(LOCAL, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 8:
                continue
            tax, gene_id, go_id, evidence, qualifier, go_term, pubmed, category = row[:8]
            code = taxmap.get(tax)
            if code is None:
                continue
            if evidence not in EXPERIMENTAL:
                continue
            if "NOT" in qualifier:  # drop negated annotations
                continue

            go_id = alt_ids.get(go_id, go_id)  # normalize secondary ids
            table = pairs[code]
            # the direct term plus every ancestor (true-path propagation)
            for term in (go_id, *ancestors(go_id)):
                record(table, term, gene_id, names.get(term, go_term),
                       category, evidence, qualifier, pubmed)

    for code, table in pairs.items():
        symbols = symbol_maps[code]

        # Drop terms that ended up annotated to more than MAX_GENES genes
        # (ontology roots and other near-root, uninformative terms).
        dropped_terms = set()
        if MAX_GENES is not None:
            gene_count = {}
            for go_id, _ in table:
                gene_count[go_id] = gene_count.get(go_id, 0) + 1
            dropped_terms = {g for g, n in gene_count.items() if n > MAX_GENES}
            if dropped_terms:
                dropped_edges = sum(gene_count[g] for g in dropped_terms)
                print(f"{code}: dropping {len(dropped_terms)} terms with >{MAX_GENES} "
                      f"genes ({dropped_edges} edges)")

        # One output file per GO aspect (BP / MF / CC).
        files, writers = {}, {}
        for suffix in ASPECTS.values():
            f = open(f"input_{code}_{suffix}.txt", "w", newline="", encoding="utf-8")
            files[suffix] = f
            w = csv.writer(f, delimiter="\t")
            w.writerow(
                ["go_id", "go_term", "gene_id", "gene_label",
                 "evidence", "aspect", "qualifier", "pubmed"]
            )
            writers[suffix] = w

        written = {s: 0 for s in ASPECTS.values()}
        genes = {s: set() for s in ASPECTS.values()}
        for (go_id, gene_id), e in table.items():
            if go_id in dropped_terms:
                continue
            suffix = ASPECTS.get(e["aspect"])
            if suffix is None:  # unknown category — skip
                continue
            label = symbols.get(gene_id, gene_id)
            # pubmed is the last column: Java's String.split("\t") drops
            # trailing empty fields, so never leave it empty (use gene2go's
            # "-" placeholder) or linksetCreator reads a too-short row.
            writers[suffix].writerow([
                go_id, e["go_term"], gene_id, label,
                "|".join(e["evidence"]), e["aspect"],
                "|".join(e["qualifier"]) or "-", "|".join(e["pubmed"]) or "-",
            ])
            written[suffix] += 1
            genes[suffix].add(gene_id)

        for suffix, f in files.items():
            f.close()
            print(f"{code}/{suffix}: wrote input_{code}_{suffix}.txt "
                  f"({written[suffix]} pairs, {len(genes[suffix])} genes)")


if __name__ == "__main__":
    main()
