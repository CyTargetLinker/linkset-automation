# TFLink linkset

Transcription factor → target gene linkset, built from TFLink's small-scale
(high-confidence) `simpleFormat` files.

| | |
| --- | --- |
| Source | https://tflink.net |
| License | "TFLink is freely available for non-commercial use" |
| DOI | not yet deposited |
| Workflow | `.github/workflows/create-linkset-tflink.yml` |

## Species

Six species have small-scale data and are built: human, mouse, rat, fruit fly,
*C. elegans* and yeast. Zebrafish is excluded because TFLink has no small-scale
data for it.

Unlike the other resources, **both ends of an edge are genes** — the TF and its
target — so both are mapped through BridgeDb.

## Confidence level

`CONFIDENCE = "SS"` (small-scale) is the default and is what the high-confidence
linkset is built from. Set it to `"LS"` or `"All"` to build large-scale variants
instead, but note the size: human large-scale is roughly 6.7 million
interactions.

Most files are served as plain `.tsv`; the large tables (LS/All for human and
mouse) are gzipped, so `tflink.py` falls back to `.tsv.gz` when the plain file is
absent.

## Versioning

TFLink is frozen at v1.0 (2022). Bump `VERSION` in `tflink.py` when a new release
appears, and `version=` in the six `tflink_*.config` files to match.

## Not done yet

There is no Zenodo deposit step — the workflow builds and QCs all six species
and uploads them as artifacts, but nothing is archived or gets a DOI.

Before anything is deposited, the license needs a decision. "Free for
non-commercial use" is materially more restrictive than the CC licenses on the
other resources, and it is not obvious that it is compatible with an open Zenodo
deposit.
