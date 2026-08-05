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

Five species are built: human, mouse, yeast, fruit fly and *C. elegans*.
Small-scale interaction counts in TFLink v1.0:

| Species | Code | SS interactions |
| --- | --- | --- |
| *Homo sapiens* | `hsa` | 16,634 |
| *Mus musculus* | `mmu` | 8,687 |
| *Saccharomyces cerevisiae* | `sce` | 5,349 |
| *Drosophila melanogaster* | `dme` | 699 |
| *Caenorhabditis elegans* | `cel` | 109 |
| *Rattus norvegicus* | — | 8 (excluded) |
| *Danio rerio* | — | none |

Rat is excluded because 8 interactions build a 12-node network rather than a
usable linkset — not for lack of BridgeDb coverage, which is fine for rat. Add
it back to `SPECIES` in `tflink.py`, with a matching config, if a later TFLink
release carries real small-scale rat data. Zebrafish has no small-scale data at
all.

Unlike the other resources, **both ends of an edge are genes** — the TF and its
target — so both are mapped through BridgeDb.

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

### Record

Bootstrapped and published 2026-08-05, in the `cytargetlinker` community.

| | |
| --- | --- |
| Concept DOI | [10.5281/zenodo.21804828](https://doi.org/10.5281/zenodo.21804828) — always the latest version |
| First version | `10.5281/zenodo.21804829` (`v1.0-SS`) |
| Deposition ID | `21804829` |
| Concept ID | `21804828` |

Both IDs are set in the workflow `env:` block, so later runs version the record
by themselves, exactly as the WikiPathways one does.

The account that owns the record is the one `ZENODO_ACCESS_TOKEN` belongs to —
the same account that owns the WikiPathways deposit. Drafts are only visible to
that account, so check https://zenodo.org/me/uploads while logged in as it
rather than looking for a draft URL.

The `bootstrap_zenodo` input that created the record refuses to run now that
`ZENODO_DEPOSITION_ID` is set, so it cannot mint a second record by accident.

## BridgeDb downloads

Bridge files come through `scripts/fetch-bridge.sh`, not a plain `wget`.
Figshare answers `HTTP 202 Accepted` with an empty body when it decides the
caller is a browser, and since 202 is a success code a naive download writes a
zero-byte file and carries on. The script accepts only an explicit 200 that
passes a zip test.

Each species' bridge file is deleted straight after its linkset is built —
six of them do not fit on a runner alongside the outputs.
