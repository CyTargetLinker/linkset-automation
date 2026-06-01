"""GO -> gene preprocessing for CyTargetLinker linksets.

Downloads NCBI gene2go and writes one tab-delimited input_<code>.txt per species
in the column order expected by go_<code>.config:

    go_id  go_term  gene_id  gene_id(label)  evidence  aspect  qualifier  pubmed

Only experimental evidence codes are kept (curated experimental annotations);
negated (NOT-qualified) annotations are dropped. gene2go carries Entrez Gene IDs
(BridgeDb syscode L), so the linkset matches the WikiPathways gene identifiers.

gene2go has no release version, so VERSION is just a date stamp for the linkset
name. To add species, append (tax_id, code) tuples to SPECIES and add a config.
"""

import csv
import gzip
import os

import requests

VERSION = "20260601"
GENE2GO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
LOCAL = "gene2go.gz"

# (NCBI tax_id, CyTargetLinker short code) — human only for now
SPECIES = [
    ("9606", "hsa"),
]

# GO experimental evidence codes (core experimental + high-throughput experimental).
# See http://geneontology.org/docs/guide-go-evidence-codes/
EXPERIMENTAL = {
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",   # core experimental
    "HTP", "HDA", "HMP", "HGI", "HEP",          # high-throughput experimental
}

# gene2go columns (tab-separated; header line starts with '#')
# 0 tax_id  1 GeneID  2 GO_ID  3 Evidence  4 Qualifier  5 GO_term  6 PubMed  7 Category


def download():
    """Download gene2go.gz, reusing an existing local copy if present."""
    if os.path.exists(LOCAL) and os.path.getsize(LOCAL) > 0:
        print(f"using existing {LOCAL}")
        return LOCAL
    print(f"downloading {GENE2GO_URL}")
    with requests.get(GENE2GO_URL, stream=True, timeout=600) as resp:
        resp.raise_for_status()
        with open(LOCAL, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    print(f"  saved {LOCAL}")
    return LOCAL


def main():
    local = download()
    writers = {}
    files = {}
    counts = {code: 0 for _, code in SPECIES}
    taxmap = {tax: code for tax, code in SPECIES}

    for code in taxmap.values():
        f = open(f"input_{code}.txt", "w", newline="", encoding="utf-8")
        files[code] = f
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["go_id", "go_term", "gene_id", "gene_label",
             "evidence", "aspect", "qualifier", "pubmed"]
        )
        writers[code] = w

    # Stream the gzip line-by-line; the full file is ~10 GB uncompressed.
    with gzip.open(local, "rt", encoding="utf-8") as fh:
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
            writers[code].writerow(
                [go_id, go_term, gene_id, gene_id, evidence, category, qualifier, pubmed]
            )
            counts[code] += 1

    for code, f in files.items():
        f.close()
        print(f"{code}: wrote input_{code}.txt ({counts[code]} interactions)")


if __name__ == "__main__":
    main()
