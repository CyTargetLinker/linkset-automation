# Reactome pathway-gene linkset

`reactome.py` builds one linkset input per species listed in `species.tsv`, plus
`build/version.txt` and `build/manifest.tsv`. Per-species configs are generated
into `reactome/config/` from `reactome.config.template` and are gitignored.

```bash
python3 reactome/reactome.py --species hsa
python3 reactome/reactome.py --version 97     # a specific release
```

Human only for now. Release 097 gives 2310 pathways, 11337 genes and 48811
edges, with 99.8% of gene nodes mapped.

## Gotchas

**The input is `NCBI2Reactome.txt`, not `ReactomePathways.gmt`.** The GMT looks
like the WikiPathways analogue and is not: it is human-only, the bytes served at
the `.gmt` URL are a ZIP archive, and the plain `.gmt` path 403s in the versioned
archive, so it cannot support a pinned build. `NCBI2Reactome.txt` carries Entrez
ids for every species and lands on the existing syscode `L`.

**A (gene, pathway) pair can appear twice**, once `TAS` and once `IEA` — 4066 of
them in human, an 8% edge inflation if written straight out. The pair is emitted
once with the codes aggregated (`TAS`, `IEA`, `TAS+IEA`), so the evidence column
is not a copy of the source column.

**Genes NCBI has not named are rejected, not kept.** WikiPathways falls back to
the Entrez id as a label; here the id is dropped, because the ones `gene_info`
lacks are pathogen-side genes Reactome files under the host species in its
infectious-disease pathways (*Escherichia ruysiae* `skp`, `secY`, `secE`) plus 22
GenBank/RefSeq accessions that are not genes. None map in BridgeDb. The cost is
real and counted: on release 097 this drops 3.4% of human genes and empties 30
pathways, so `max_reject_frac` in `species.tsv` gates it.

**Versions are zero-padded to three digits.** The Zenodo gate compares version
strings with the shell's `\<`, under which `100` sorts before `97`.

**"Lowest level" means lowest-level diagram, not leaf.** 325 of the 2343 human
pathways still have children, so about 9% of pairs attribute a gene to both a
pathway and an ancestor of it. The all-levels file is 3x larger and gives
*Signal Transduction* 2622 genes, which is useless as a neighbour set.

**Non-human species are mostly orthology projections** and are all `off` in
`species.tsv` pending a decision — see the comments there, which record which
species are genuinely curated.
