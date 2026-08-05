#!/usr/bin/env python3
"""Quality-control validator for CyTargetLinker XGMML linkset files.

Checks that an XGMML file is well-formed and has the structure every linkset
produced by the linkset-creator carries:

  - root <graph> (XGMML namespace http://www.cs.rpi.edu/XGMML)
  - each <node> has an id, a non-empty `identifiers` list att, a `type` att,
    and a `label` att
  - each <edge> has source/target referencing existing nodes, plus
    `datasource` and `interaction` atts

By default an empty graph (0 nodes or 0 edges) FAILS; pass --allow-empty to
permit it. --min-nodes / --min-edges add explicit thresholds.

--mapped-type TYPE together with --min-mapped-frac F guards against a silent
BridgeDb mapping failure: among nodes whose `type` att equals TYPE, it fails
unless at least fraction F carry more than one identifier (i.e. got expanded to
cross-references). A mapper that failed to load leaves every node with only its
input id, which this catches while a non-empty `identifiers` list alone does not.
The observed fraction is always reported when --mapped-type is given, so the
mapping yield of a new species can be measured before a threshold is chosen.

--labelled-type TYPE reports what fraction of nodes of that type carry a label
that is not simply their own id, i.e. how many got a real name rather than a
bare identifier. It is a probe, not a gate: there is no threshold, and unlike
--mapped-type it stays silent when a file has no nodes of that type.

Usage:
    python3 qc_xgmml.py [--allow-empty] [--min-nodes N] [--min-edges N]
                        [--mapped-type TYPE --min-mapped-frac F]
                        [--labelled-type TYPE] FILE_OR_GLOB...

Each positional argument is treated as a glob (literal paths work too).
Exit code is 0 only if every matched file passes; 1 otherwise.
"""

import argparse
import glob
import sys
import xml.etree.ElementTree as ET

XGMML_NS = "http://www.cs.rpi.edu/XGMML"


def _local(tag):
    """Strip a `{namespace}` prefix from an ElementTree tag."""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _att_index(element):
    """Map immediate child <att> elements by their `name`, plus the list of them.

    Returns (by_name, atts) where by_name maps name -> the <att> element
    (last wins on duplicates) and atts is the ordered list of <att> children.
    """
    atts = [c for c in element if _local(c.tag) == "att"]
    by_name = {a.get("name"): a for a in atts if a.get("name") is not None}
    return by_name, atts


def validate_file(path, allow_empty=False, min_nodes=0, min_edges=0,
                  mapped_type=None, min_mapped_frac=0.0, labelled_type=None):
    """Validate one XGMML file.

    Returns (errors, warnings, infos, n_nodes, n_edges).
    """
    errors = []
    warnings = []
    infos = []
    mapped_total = 0     # nodes whose type == mapped_type
    mapped_ok = 0        # ...of those, with >1 identifier (got cross-references)
    labelled_total = 0   # nodes whose type == labelled_type
    labelled_ok = 0      # ...of those, whose label is not just their id

    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        return ([f"not well-formed XML: {e}"], [], [], 0, 0)
    except OSError as e:
        return ([f"cannot read file: {e}"], [], [], 0, 0)

    if _local(root.tag) != "graph":
        return ([f"root element is <{_local(root.tag)}>, expected <graph>"], [], [], 0, 0)
    if not root.tag.startswith("{" + XGMML_NS + "}"):
        warnings.append(f"root namespace is not the XGMML namespace ({XGMML_NS})")

    graph_atts, _ = _att_index(root)
    for required in ("LinkSet Name", "Creation Date"):
        if required not in graph_atts:
            warnings.append(f"graph is missing the '{required}' attribute")

    nodes = [c for c in root if _local(c.tag) == "node"]
    edges = [c for c in root if _local(c.tag) == "edge"]

    node_ids = set()
    for node in nodes:
        nid = node.get("id")
        if not nid:
            errors.append("a <node> is missing its 'id' attribute")
            continue
        if nid in node_ids:
            warnings.append(f"duplicate node id '{nid}'")
        node_ids.add(nid)

        by_name, _ = _att_index(node)
        ident = by_name.get("identifiers")
        n_ident = 0
        if ident is None:
            errors.append(f"node '{nid}' has no 'identifiers' attribute")
        else:
            values = [c for c in ident if _local(c.tag) == "att"]
            n_ident = len(values)
            if not values:
                errors.append(f"node '{nid}' has an empty 'identifiers' list")
        if "type" not in by_name:
            errors.append(f"node '{nid}' has no 'type' attribute")
        else:
            node_type = by_name["type"].get("value")
            if mapped_type is not None and node_type == mapped_type:
                mapped_total += 1
                if n_ident > 1:
                    mapped_ok += 1
            if labelled_type is not None and node_type == labelled_type:
                labelled_total += 1
                label_att = by_name.get("label")
                label = (label_att.get("value") if label_att is not None
                         else node.get("label"))
                if label and label != nid:
                    labelled_ok += 1
        if "label" not in by_name and node.get("label") is None:
            warnings.append(f"node '{nid}' has no 'label'")

    for edge in edges:
        eid = edge.get("id") or "?"
        src = edge.get("source")
        tgt = edge.get("target")
        if not src or not tgt:
            errors.append(f"edge '{eid}' is missing source/target")
        else:
            if src not in node_ids:
                errors.append(f"edge '{eid}' references missing source node '{src}'")
            if tgt not in node_ids:
                errors.append(f"edge '{eid}' references missing target node '{tgt}'")
        by_name, _ = _att_index(edge)
        if "datasource" not in by_name:
            errors.append(f"edge '{eid}' has no 'datasource' attribute")
        if "interaction" not in by_name:
            warnings.append(f"edge '{eid}' has no 'interaction' attribute")

    n_nodes, n_edges = len(nodes), len(edges)
    if n_nodes == 0 or n_edges == 0:
        if allow_empty:
            warnings.append(f"empty linkset ({n_nodes} nodes, {n_edges} edges)")
        else:
            errors.append(
                f"empty linkset ({n_nodes} nodes, {n_edges} edges); "
                "use --allow-empty to permit"
            )
    if n_nodes < min_nodes:
        errors.append(f"only {n_nodes} nodes (min-nodes={min_nodes})")
    if n_edges < min_edges:
        errors.append(f"only {n_edges} edges (min-edges={min_edges})")

    if mapped_type is not None:
        if mapped_total == 0:
            errors.append(f"no nodes of type '{mapped_type}' to check mapping coverage")
        else:
            frac = mapped_ok / mapped_total
            infos.append(
                f"{mapped_ok}/{mapped_total} ({frac:.1%}) '{mapped_type}' nodes "
                "have cross-references"
            )
            if frac < min_mapped_frac:
                errors.append(
                    f"only {mapped_ok}/{mapped_total} ({frac:.1%}) '{mapped_type}' nodes "
                    f"have cross-references (min-mapped-frac={min_mapped_frac:.0%}); "
                    "BridgeDb mapping may have failed"
                )

    # Report-only: no threshold. How many nodes carry a real name rather than a
    # bare identifier varies by species far too much for one floor to be useful
    # (measured 81% to 100%), and the pipeline that builds the labels already
    # fails earlier, with a better message, when the source is broken.
    if labelled_type is not None and labelled_total:
        frac = labelled_ok / labelled_total
        infos.append(
            f"{labelled_ok}/{labelled_total} ({frac:.1%}) '{labelled_type}' nodes "
            "have a label that is not just their id"
        )

    return (errors, warnings, infos, n_nodes, n_edges)


def main(argv=None):
    parser = argparse.ArgumentParser(description="QC validator for CyTargetLinker XGMML linksets.")
    parser.add_argument("paths", nargs="+", help="XGMML file paths or globs to validate")
    parser.add_argument("--allow-empty", action="store_true", help="treat empty graphs as a warning, not an error")
    parser.add_argument("--min-nodes", type=int, default=0, help="fail if fewer than this many nodes")
    parser.add_argument("--min-edges", type=int, default=0, help="fail if fewer than this many edges")
    parser.add_argument("--mapped-type", default=None,
                        help="node 'type' value to report BridgeDb mapping coverage for (e.g. gene)")
    parser.add_argument("--min-mapped-frac", type=float, default=0.0,
                        help="fail unless this fraction of --mapped-type nodes have >1 identifier")
    parser.add_argument("--labelled-type", default=None,
                        help="node 'type' value to report label coverage for (e.g. gene); "
                             "reports only, never fails")
    args = parser.parse_args(argv)

    files = []
    for pattern in args.paths:
        files.extend(sorted(glob.glob(pattern, recursive=True)))
    # de-duplicate while preserving order
    seen = set()
    files = [f for f in files if not (f in seen or seen.add(f))]

    if not files:
        print(f"ERROR: no files matched: {' '.join(args.paths)}", file=sys.stderr)
        return 1

    failed = 0
    for path in files:
        errors, warnings, infos, n_nodes, n_edges = validate_file(
            path, allow_empty=args.allow_empty, min_nodes=args.min_nodes, min_edges=args.min_edges,
            mapped_type=args.mapped_type, min_mapped_frac=args.min_mapped_frac,
            labelled_type=args.labelled_type
        )
        status = "PASS" if not errors else "FAIL"
        if errors:
            failed += 1
        print(f"{status}  {path}  ({n_nodes} nodes, {n_edges} edges)")
        for i in infos:
            print(f"        INFO: {i}")
        for e in errors:
            print(f"        ERROR: {e}")
        for w in warnings:
            print(f"        WARNING: {w}")

    passed = len(files) - failed
    print(f"\nSummary: {passed} passed, {failed} failed (of {len(files)} file(s))")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
