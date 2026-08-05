# WikiPathways linkset

Pathway–gene linksets for 13 species, built from the WikiPathways GMT release.

| | |
| --- | --- |
| Source | https://www.wikipathways.org (GMT files from https://data.wikipathways.org) |
| License | CC BY 4.0 as deposited; WikiPathways content itself is CC0 1.0 |
| DOI | [10.5281/zenodo.4500957](https://doi.org/10.5281/zenodo.4500957) (concept DOI, resolves to the latest version) |
| Workflow | `.github/workflows/create-linkset-wikipathways.yml` |
| Schedule | quarterly — 03:17 UTC on the 15th of Jan / Apr / Jul / Oct |
| Species | hsa mmu rno bta sce cel dre ptr cfa gga dme ath aga |

This is the only resource with a full pipeline: build, QC, and deposit to
Zenodo.

## How it works

`wp.py` resolves the WikiPathways release, downloads one GMT per species, and
writes `input_<code>.txt` plus a generated `config/wp_<code>.config` for each.
The GMT set name packs four fields separated by `%`:
`pathway name%WikiPathways_<date>%WP####%species`.

Genes enter as Entrez (`target_syscode_in=L` — every species' GMT uses Entrez
ids, including the plants and insects) and are expanded through BridgeDb to
Ensembl, Entrez and UniProt, plus HGNC for human.

It takes two optional arguments:

```bash
python3 wikipathways/wp.py                        # all species, newest release
python3 wikipathways/wp.py --species sce aga      # just these two
python3 wikipathways/wp.py --version 20260510     # a specific release
```

`wp.py` also writes `build/version.txt` and `build/manifest.tsv`. The manifest
lists the species that were **actually built**, and is what every later workflow
step iterates — so a species skipped during preprocessing can never leave a
later step looking for a file that was never produced.

## Updating to a new release

**Nothing to do.** The release is resolved on every run and flows into the GMT
URLs, the `version=` line of every generated config, and the Zenodo record
version.

`data.wikipathways.org` keeps only a **rolling 12-month window** of monthly
releases, so a hardcoded release would eventually 404. `wp.py` therefore walks
the **dated release directories** on <https://data.wikipathways.org/>
(`20260710/`, `20260610/`, …) newest first.

It deliberately does not use `/current/`: that is a moving alias and says
nothing about whether the release behind it is complete, whereas a dated
directory holding a GMT for **every `core` species** is positive evidence that
the release shipped. (The stray `20230810` directory is missing mouse and rat,
and is correctly rejected on exactly this test.) WikiPathways releases on the
10th of the month, and a release dated otherwise is flagged in the log.

If the newest directory is incomplete — a release still being published — it is
skipped with a warning and the previous one used. That is safe because the
workflow refuses to publish a release older than the one already on Zenodo.
Three directories are tried before giving up, since more than a couple of
incomplete releases in a row is a problem to report rather than work around.

Downloads then come from the resolved `/<version>/gmt/` path, so a release
appearing mid-run cannot produce a linkset assembled from two releases.

## Adding or removing a species

Everything is driven by [`species.tsv`](species.tsv) — adding a species is one
line there, and nothing else changes.

```
code	wp_token	bridgedb_species	syscodes_out	tier	min_mapped_frac
hsa	Homo_sapiens	Homo sapiens	En,H,L,S	core	0.90
```

- **`wp_token`** — the species as spelled in the GMT filename.
- **`bridgedb_species`** — the key in bridgedb/data `gene.json`. Stated
  explicitly rather than derived from `wp_token`, because the two sources can
  drift (Ensembl now calls the dog *Canis lupus familiaris*) and because the
  `.bridge` filenames are not derivable from the species name anyway
  (*Equus caballus* is `Qc_`, *Populus trichocarpa* is `Pi_`).
- **`syscodes_out`** — BridgeDb systems to expand gene ids into. `H` (HGNC) is
  human-only.
- **`tier`** — `core` means a missing GMT fails the run, `extra` means it is
  skipped with a warning, `off` means it is not built.
- **`min_mapped_frac`** — the QC floor for BridgeDb coverage (below).

Five species are listed but commented `off` because WikiPathways has too little
content for them to be useful linksets — pig, horse, poplar, tomato and maize
each have between one and five pathways. `wp.py` also skips any `extra` species
below `MIN_PATHWAYS` (5), so one that shrinks again drops out on its own.

## Configs are generated

There are no committed per-species configs. Each `config/wp_<code>.config` is
generated from [`wp.config.template`](wp.config.template) on every run, and the
whole `config/` directory is gitignored. This is what removes the need to
hand-maintain the release version across a dozen files.

The generator asserts that every line of a rendered config is exactly one
`key=value` pair. LinksetCreator's `ConfigFileReader` splits on the first `=`
but keeps only lines that split into exactly two parts, silently discarding
anything else — so a value containing `=` would vanish without an error.

## Quality control

Each linkset is validated by `scripts/qc_xgmml.py` before anything is published.
Beyond the structural checks it verifies **BridgeDb mapping coverage**:
`--mapped-type gene --min-mapped-frac` fails a linkset unless that fraction of
gene nodes gained more than one identifier. A mapper that fails to load leaves
every node carrying only its input id — structurally valid, but useless for
matching a user's network.

Measured on release 20260710: 100% for most species, 98.9% cow, 95.9% rat. The
floors are 0.90 for `core` and 0.80 for `extra`. Arabidopsis and Anopheles reach
100% despite their Ensembl Plants / Metazoa 56 mapping files, which had been
expected to map the Entrez-keyed GMTs poorly.

**Any QC failure fails the whole run and publishes nothing.** A partial publish
would mint a permanent DOI on an incomplete record; a failed run just needs
re-dispatching.

## Publishing

Three gates protect the Zenodo record, all because published versions are
immutable. The run will not publish if:

1. the concept record already publishes this WikiPathways release — otherwise a
   quarterly schedule would mint a new DOI whether or not anything changed;
2. the built release is **older** than the published one, which could otherwise
   happen after falling back from an incomplete release;
3. a file present in the last published version is missing from this build, so a
   transient upstream 404 cannot silently drop a species.

Publishing needs the `ZENODO_ACCESS_TOKEN` secret. The `skip_zenodo` and
`force_publish` dispatch inputs make the gates testable.

## Note on BridgeDb downloads

Use `scripts/fetch-bridge.sh`; do not hand-roll a download. Figshare redirects
to a presigned S3 URL that expires in ten seconds, and it answers requests
carrying a full browser user-agent with `HTTP 202 Accepted` and an empty body.
202 is a success code, so `wget` exits 0 and writes a zero-byte file, and
`--retry-on-http-error` cannot match a 2xx — which is exactly how the previous
download failed, looping until it gave up. The script sends no spoofed
user-agent and accepts a download only on an explicit 200 that also passes a zip
test.

## Known limitation: gene labels are Entrez IDs, not symbols

The WikiPathways GMT contains only Entrez Gene IDs — there are no symbols in it.
`wp.py` therefore writes the same value into both the `gene_id` and
`gene_symbol` columns, and `target_label_column=4` points at that duplicate, so
gene nodes are labelled with numbers rather than symbols in Cytoscape.

For contrast, the GO linkset avoids this by joining symbols in from NCBI
`gene_info`, and an earlier SPARQL-based WikiPathways run produced real symbols
because the query returned a `GeneName` column. Fixing this needs either a
symbol lookup after the GMT parse or a return to the SPARQL endpoint for the
input file.
