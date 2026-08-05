# TFLink linkset

Transcription factor → target gene linkset, built from TFLink's small-scale
(high-confidence) `simpleFormat` files.

| | |
| --- | --- |
| Source | https://tflink.net |
| License | CC BY-NC 4.0 (non-commercial, attribution required) |
| DOI | not yet deposited — needs the one-time setup below |
| Workflow | `.github/workflows/create-linkset-tflink.yml` |

## Licensing

TFLink's FAQ states only that "TFLink is freely available for non-commercial
use", with no license named on the site or the download page. The TFLink paper
(Liska *et al.* 2022, [10.1093/database/baac083](https://doi.org/10.1093/database/baac083))
is published under **CC BY-NC 4.0**, which permits non-commercial re-use,
distribution and reproduction provided the original work is cited.

The linkset is therefore deposited as **CC BY-NC 4.0** — the most permissive
license consistent with the source. Note this is more restrictive than the
CC BY 4.0 used for WikiPathways, and an NC license is not "open" in the
Open Definition sense. If TFLink's terms ever need to be relied on more firmly
than a paper's license statement, the FAQ suggests contacting the authors.

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

## Release

The workflow runs the full chain: build → QC → deposit to Zenodo, gated the same
way as WikiPathways. It will not publish when `skip_zenodo` is set, when Zenodo
already carries this version (unless `force_publish`), or when the build has
lost a species relative to what is already published — Zenodo versions are
immutable, so a shrinking record is refused.

The deposit version string is `<VERSION>-<CONFIDENCE>`, e.g. `v1.0-SS`, because
the confidence level changes what the linkset contains.

### One-time setup, still outstanding

TFLink has no Zenodo record yet, so `ZENODO_DEPOSITION_ID` and
`ZENODO_CONCEPT_ID` are empty in the workflow `env:` block. Until they are set,
every run builds and QCs normally and then skips the publish steps with a notice
in the job summary — it does not fail.

To finish it:

1. Run the workflow with **`bootstrap_zenodo`** ticked. It builds and QCs as
   usual, then creates a Zenodo deposition, uploads the six linksets and sets
   the metadata — and stops. The draft is **not published**: minting a DOI is
   irreversible, so the first public record is left for a human to check and
   publish. The draft URL and both IDs are printed in the job summary.
2. Review the draft on Zenodo and publish it.
3. Put the two IDs from the job summary into `ZENODO_DEPOSITION_ID` and
   `ZENODO_CONCEPT_ID` in the workflow `env:` block.
4. Add the concept DOI to the resource table in the top-level `README.md`.

After that the workflow versions the record on its own, exactly as the
WikiPathways one does. Bootstrap refuses to run once `ZENODO_DEPOSITION_ID` is
set, so it cannot create a second record by accident.

## BridgeDb downloads

Bridge files come through `scripts/fetch-bridge.sh`, not a plain `wget`.
Figshare answers `HTTP 202 Accepted` with an empty body when it decides the
caller is a browser, and since 202 is a success code a naive download writes a
zero-byte file and carries on. The script accepts only an explicit 200 that
passes a zip test.

Each species' bridge file is deleted straight after its linkset is built —
six of them do not fit on a runner alongside the outputs.
