"""TFLink preprocessing for CyTargetLinker linksets.

Downloads the TFLink "simpleFormat" interaction table for each species and writes
one tab-delimited input_<code>.txt per species in the column order expected by
tflink_<code>.config:

    NCBI.GeneID.TF  Name.TF  NCBI.GeneID.Target  Name.Target  Detection.method  PubmedID

The linkset-creator jar then turns each input_<code>.txt + tflink_<code>.config into
an XGMML linkset.

CONFIDENCE selects which TFLink file is used:
  "SS"  small-scale  - manually curated, high-confidence (default; Cytoscape-friendly)
  "LS"  large-scale  - high-throughput (e.g. ChIP-seq); millions of human interactions
  "All" both combined
Zebrafish is omitted because TFLink has no small-scale (SS) data for it.
"""

import os

import pandas as pd
import requests

VERSION = "v1.0"
BASE_URL = "https://cdn.netbiol.org/tflink/download_files"
CONFIDENCE = "SS"

# (TFLink organism name as used in the file name, CyTargetLinker short code)
# Rattus norvegicus is deliberately absent: TFLink v1.0 has 8 small-scale rat
# interactions, which builds a 12-node network rather than a usable linkset.
# Zebrafish has no small-scale data at all. Add rat back if a later TFLink
# release carries real small-scale coverage for it.
SPECIES = [
    ("Homo_sapiens", "hsa"),
    ("Mus_musculus", "mmu"),
    ("Drosophila_melanogaster", "dme"),
    ("Caenorhabditis_elegans", "cel"),
    ("Saccharomyces_cerevisiae", "sce"),
]

# Columns kept from the raw TFLink table, in linkset (source, target, edge) order.
OUTPUT_COLUMNS = [
    "NCBI.GeneID.TF",
    "Name.TF",
    "NCBI.GeneID.Target",
    "Name.Target",
    "Detection.method",
    "PubmedID",
]


def download(organism):
    """Download a species' simpleFormat table, returning the local file name.

    Most files are served as plain .tsv; the large (LS/All human and mouse) tables
    are served gzipped, so fall back to .tsv.gz when the plain file is absent.
    """
    stem = f"TFLink_{organism}_interactions_{CONFIDENCE}_simpleFormat_{VERSION}.tsv"
    # reuse a previously downloaded file if one is already present
    for suffix in ("", ".gz"):
        local = stem + suffix
        if os.path.exists(local) and os.path.getsize(local) > 0:
            print(f"  using existing {local}")
            return local
    for suffix in ("", ".gz"):
        url = f"{BASE_URL}/{stem}{suffix}"
        resp = requests.get(url, stream=True, timeout=120)
        if resp.status_code == 200:
            local = stem + suffix
            with open(local, "wb") as f:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
            print(f"  downloaded {url} -> {local}")
            return local
        resp.close()
    raise RuntimeError(f"could not download TFLink {CONFIDENCE} table for {organism}")


def main():
    for organism, code in SPECIES:
        print(f"{organism} ({code}):")
        local = download(organism)
        # compression is inferred from the .gz / .tsv suffix
        df = pd.read_csv(local, sep="\t", usecols=OUTPUT_COLUMNS, dtype=str)
        print(f"  {len(df)} interactions")
        out = f"input_{code}.txt"
        df[OUTPUT_COLUMNS].to_csv(out, sep="\t", index=False)
        print(f"  wrote {out}")


if __name__ == "__main__":
    main()
