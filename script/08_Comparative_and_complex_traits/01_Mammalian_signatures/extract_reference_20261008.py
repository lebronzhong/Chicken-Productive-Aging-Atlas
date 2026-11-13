#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 08_Comparative_and_complex_traits
# Purpose: Extract reference 20261008.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--expected-reference-genes", type=int)
parser.add_argument("--expected-strict-genes", type=int)
args = parser.parse_args()
import zipfile, xml.etree.ElementTree as ET, json, hashlib
import pandas as pd, numpy as np

root = Path("/path/to/chicken-aging/analysis/conservation")
p = root / "reference/Supplementary_Table_2.xlsx"
ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
with zipfile.ZipFile(p) as z:
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        shared = ["".join(si.itertext()) for si in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    sheets = wb.find("m:sheets", ns)
    rels = {
        r.attrib["Id"]: r.attrib["Target"]
        for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    }
    parts = []
    audit = []
    for idx, sp in zip(range(14, 18), ["Mouse", "Rat", "Macaque", "Human"]):
        sh = sheets[idx]
        rel = sh.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        target = rels[rel]
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        rows = []
        for row in ET.fromstring(z.read(target)).findall(".//m:sheetData/m:row", ns):
            cells = {}
            for c in row:
                letters = "".join(x for x in c.attrib["r"] if x.isalpha())
                col = 0
                for x in letters:
                    col = col * 26 + ord(x) - 64
                v = c.find("m:v", ns)
                v = v.text if v is not None else ""
                if c.attrib.get("t") == "s":
                    v = shared[int(v)]
                elif c.attrib.get("t") == "inlineStr":
                    v = "".join(c.find("m:is", ns).itertext())
                cells[col - 1] = v
            rows.append(cells)
        header = rows[0]
        data = pd.DataFrame([{header[k]: v for k, v in r.items() if k in header} for r in rows[1:]])
        fields = [
            "Entrez.ID",
            "Gene.symbol",
            "Slope",
            "SE",
            "Pearson.corr",
            "P.Value",
            "P.Adjusted",
        ]
        d = data[fields].copy()
        d["Entrez.ID"] = pd.to_numeric(d["Entrez.ID"]).astype(int).astype(str)
        for c in fields[2:]:
            d[c] = pd.to_numeric(d[c])
        assert not d["Entrez.ID"].duplicated().any()
        d.to_csv(root / "reference" / f"{sp}_published_aging_signature.tsv", sep="\t", index=False)
        audit.append(dict(species=sp, sheet=sh.attrib["name"], rows=len(d)))
        parts.append(d.set_index("Entrez.ID").add_prefix(sp + "_"))
m = pd.concat(parts, axis=1, join="inner")
s = m[[f"{sp}_Slope" for sp in ["Mouse", "Rat", "Macaque", "Human"]]]
pvals = m[[f"{sp}_P.Adjusted" for sp in ["Mouse", "Rat", "Macaque", "Human"]]]
m["direction"] = np.where((s > 0).all(axis=1), "up", np.where((s < 0).all(axis=1), "down", "mixed"))
m["strict"] = (pvals < 0.05).all(axis=1) & (m.direction != "mixed")
m.to_csv(root / "reference/four_species_intersection.tsv", sep="\t")
m[m.strict].to_csv(root / "reference/four_species_strict.tsv", sep="\t")
info = dict(
    source="https://www.nature.com/articles/s41586-026-10542-3",
    xlsx_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
    sheets=audit,
    common=len(m),
    strict=int(m.strict.sum()),
)
if args.expected_reference_genes is not None:
    assert len(m) == args.expected_reference_genes
if args.expected_strict_genes is not None:
    assert int(m.strict.sum()) == args.expected_strict_genes
(root / "reference/provenance.json").write_text(json.dumps(info, indent=2))
print(json.dumps(info, indent=2), flush=True)
