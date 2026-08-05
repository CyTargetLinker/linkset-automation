import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

URLS = {
    "product1": "https://www.orphadata.com/data/xml/en_product1.xml",
    "product6": "https://www.orphadata.com/data/xml/en_product6.xml",
    "product9_prev": "https://www.orphadata.com/data/xml/en_product9_prev.xml",
    "product9_ages": "https://www.orphadata.com/data/xml/en_product9_ages.xml",
}

OUTPUT_FILE = "input.txt"
PMID_RE = re.compile(r"(\d{6,9})\s*\[?\s*pmid\s*\]?", re.IGNORECASE)

HEADER = [
    "OrphaID", "DiseaseName", "GeneSymbol", "HGNC_ID", "Synonyms", "Source",
    "Prevalence", "Inheritance", "AgeOfOnset",
    "ICD_10", "ICD_11", "OMIM", "UMLS",
]


def fetch(name, url):
    dest = Path(f"{name}.xml")
    if not dest.exists():
        with requests.get(url, stream=True, timeout=180) as r:
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 16):
                    f.write(chunk)
    return dest


def text_of(elem, tag):
    child = elem.find(tag)
    return child.text.strip() if child is not None and child.text else None


def all_texts(elem, path):
    return [c.text.strip() for c in elem.findall(path) if c is not None and c.text]


def parse_cross_references(xml_path):
    data = {}

    for _, elem in ET.iterparse(xml_path, events=("end",)):
        if elem.tag != "Disorder":
            continue

        code = text_of(elem, "OrphaCode")
        if not code:
            elem.clear()
            continue

        refs = {"ICD-10": [], "ICD-11": [], "OMIM": [], "UMLS": []}

        xref_list = elem.find("ExternalReferenceList")
        if xref_list is not None:
            for xref in xref_list.findall("ExternalReference"):
                source = text_of(xref, "Source")
                reference = text_of(xref, "Reference")
                if source in refs and reference:
                    refs[source].append(reference)

        data[code] = {
            "disease_name": text_of(elem, "Name"),
            "synonyms": all_texts(elem, "SynonymList/Synonym"),
            **refs,
        }

        elem.clear()

    return data


def parse_prevalence(xml_path):
    data = {}

    for _, elem in ET.iterparse(xml_path, events=("end",)):
        if elem.tag != "Disorder":
            continue

        code = text_of(elem, "OrphaCode")
        if not code:
            elem.clear()
            continue

        entries = []
        prev_list = elem.find("PrevalenceList")

        if prev_list is not None:
            for prev in prev_list.findall("Prevalence"):
                p_class = text_of(prev, "PrevalenceClass/Name")
                p_geo = text_of(prev, "PrevalenceGeographic/Name")

                if p_class:
                    entries.append({
                        "class": p_class,
                        "geographic": p_geo,
                    })

        if entries:
            data[code] = entries

        elem.clear()

    return data


def summarize_prevalence(entries):
    if not entries:
        return None

    worldwide = [
        e for e in entries
        if (e.get("geographic") or "").lower() == "worldwide"
    ]

    pool = worldwide if worldwide else entries

    classes = []
    for e in pool:
        if e["class"] not in classes:
            classes.append(e["class"])

    return "; ".join(classes)


def parse_natural_history(xml_path):
    data = {}

    for _, elem in ET.iterparse(xml_path, events=("end",)):
        if elem.tag != "Disorder":
            continue

        code = text_of(elem, "OrphaCode")
        if not code:
            elem.clear()
            continue

        inheritance = all_texts(
            elem,
            "TypeOfInheritanceList/TypeOfInheritance/Name"
        )

        onset = all_texts(
            elem,
            "AverageAgeOfOnsetList/AverageAgeOfOnset/Name"
        )

        if inheritance or onset:
            data[code] = {
                "inheritance": inheritance,
                "onset": onset,
            }

        elem.clear()

    return data


def parse_genes(xml_path):
    for _, elem in ET.iterparse(xml_path, events=("end",)):
        if elem.tag != "Disorder":
            continue

        code = text_of(elem, "OrphaCode")
        name = text_of(elem, "Name")
        genes = []
        pmids = []

        assoc_list = elem.find("DisorderGeneAssociationList")

        if assoc_list is not None:
            for assoc in assoc_list.findall("DisorderGeneAssociation"):
                gene_el = assoc.find("Gene")
                if gene_el is None:
                    continue

                symbol = text_of(gene_el, "Symbol")
                xrefs = {}

                xref_list = gene_el.find("ExternalReferenceList")
                if xref_list is not None:
                    for xref in xref_list.findall("ExternalReference"):
                        source = text_of(xref, "Source")
                        reference = text_of(xref, "Reference")

                        if source and reference:
                            xrefs.setdefault(source, []).append(reference)

                if symbol:
                    genes.append({
                        "symbol": symbol,
                        "xrefs": xrefs,
                    })

                source_val = text_of(assoc, "SourceOfValidation")
                if source_val:
                    pmids.extend(PMID_RE.findall(source_val))

        if code and name and genes:
            yield {
                "orpha_code": code,
                "disease_name": name,
                "genes": genes,
                "pmids": sorted(set(pmids), key=int),
            }

        elem.clear()


def build_rows(rec):
    disease_name = rec.get("disease_name") or "-"
    synonyms = "; ".join(rec.get("synonyms") or []) or "-"

    pmids = rec.get("pmids") or []
    source = "; ".join(
        f"PMID:{p} (https://pubmed.ncbi.nlm.nih.gov/{p})"
        for p in pmids
    ) or "-"

    prevalence = summarize_prevalence(rec.get("prevalence")) or "-"
    inheritance = "; ".join(rec.get("inheritance") or []) or "-"
    onset = "; ".join(rec.get("onset") or []) or "-"

    icd10 = "; ".join(rec.get("ICD-10") or []) or "-"

    icd11_codes = rec.get("ICD-11") or []
    icd11 = "; ".join(
        f"{c} (https://icd.who.int/browse/latest-release/mms/search?q={c})"
        for c in icd11_codes
    ) or "-"

    omim_codes = rec.get("OMIM") or []
    omim = "; ".join(
        f"{o} (https://omim.org/entry/{o})"
        for o in omim_codes
    ) or "-"

    umls = "; ".join(rec.get("UMLS") or []) or "-"

    for gene in rec["genes"]:
        hgnc_ids = gene.get("xrefs", {}).get("HGNC", [])
        hgnc_id = "; ".join(hgnc_ids) if hgnc_ids else "-"

        yield [
            f"ORPHA:{rec['orpha_code']}",
            disease_name,
            gene["symbol"],
            hgnc_id,
            synonyms,
            source,
            prevalence,
            inheritance,
            onset,
            icd10,
            icd11,
            omim,
            umls,
        ]


def main():
    paths = {
        name: fetch(name, url)
        for name, url in URLS.items()
    }

    cross_refs = parse_cross_references(paths["product1"])
    prevalence = parse_prevalence(paths["product9_prev"])
    natural_history = parse_natural_history(paths["product9_ages"])
    records = list(parse_genes(paths["product6"]))

    for rec in records:
        code = rec["orpha_code"]

        rec.update(cross_refs.get(code, {}))

        rec["prevalence"] = prevalence.get(code)

        nh = natural_history.get(code, {})
        rec["inheritance"] = nh.get("inheritance")
        rec["onset"] = nh.get("onset")

    with open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(HEADER)

        for rec in records:
            for row in build_rows(rec):
                writer.writerow(row)


if __name__ == "__main__":
    main()