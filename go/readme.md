# Gene Ontology linkset

GO term–gene linkset for human, built from NCBI gene2go and the GO ontology.
Written as three separate linksets, one per GO aspect.

| | |
| --- | --- |
| Source | https://geneontology.org (annotations via NCBI gene2go, ontology via go-basic.obo) |
| License | CC BY 4.0 — see the [GO citation policy](https://geneontology.org/docs/go-citation-policy/) |
| DOI | not yet deposited |
| Workflow | `.github/workflows/create-linkset-go.yml` |

## What is included

Only experimental evidence codes are kept, so these are curated experimental
annotations rather than the full electronic set. Negated (`NOT`-qualified)
annotations are dropped. gene2go carries Entrez Gene IDs (BridgeDb syscode `L`),
which is what makes the GO linkset line up with the WikiPathways one.

## Ontology propagation (the true-path rule)

gene2go records only the most specific annotation for each gene. By the GO
true-path rule a gene annotated to a term is implicitly annotated to all of that
term's ancestors, so a direct-only linkset leaves high-level terms looking almost
empty — "response to decreased oxygen levels" had one direct human gene but
hundreds once descendants were counted.

`go.py` therefore downloads `go-basic.obo` and, for every direct annotation, also
emits an edge to each ancestor reachable via `is_a` or `part_of`. **`regulates`
is excluded on purpose**: a regulator of process X is not a participant in X, so
propagating across it would be biologically wrong.

Two consequences of propagation:

- After propagation the near-root terms link to most of the genome, so terms with
  more than `MAX_GENES` (500) genes are dropped as uninformative. Set it to
  `None` to keep everything.
- The three aspects are written separately (`bp` / `mf` / `cc`) because
  propagation never crosses namespaces, so each is self-contained.

Each GO term–gene pair is written once. A pair reachable from several direct
annotations has the evidence codes, qualifiers and PubMed IDs of every
contributing annotation merged into one edge, `|`-joined.

## Gene labels

gene2go has no gene symbol, so the readable label comes from NCBI `gene_info`
(GeneID → Symbol). Without it the genes would appear in Cytoscape as bare Entrez
numbers. If a symbol is missing, the Entrez ID is used instead.

## Versioning and attribution

GO's citation policy asks redistributed data to state the release it came from.
`go.py` reads the `data-version` header from `go-basic.obo` and `check_version`
warns when `VERSION` has drifted from the ontology actually downloaded. The
configs are static files, so correcting a drift is a manual edit: bump `VERSION`
in `go.py` and `version=` in all three `go_hsa_*.config` files.

## Not done yet

There is no Zenodo deposit step — the workflow builds and QCs the three
linksets and uploads them as artifacts, but nothing is archived or gets a DOI.

## Adding a species

Append a `(tax_id, code, gene_info_url)` entry to `SPECIES` in `go.py` and add a
matching `go_<code>_<aspect>.config` for each of the three aspects.
