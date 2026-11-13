#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 01_Transcriptome
# Purpose: Download 7w.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

import urllib.request, hashlib, json, gzip, subprocess, time
from pathlib import Path

root = Path(__file__).resolve().parent / "resources_7w"
root.mkdir(exist_ok=True)
base = "https://ftp.ensembl.org/pub/release-114/gtf/gallus_gallus_gca016700215v2/"
name = "Gallus_gallus_gca016700215v2.bGalGal1.pat.whiteleghornlayer.GRCg7w.114.gtf.gz"
report = "GCA_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_WZ"
urls = {
    name: base + name,
    "CHECKSUMS": base + "CHECKSUMS",
    "README": base + "README",
    "GRCg7w_assembly_report.txt": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/700/215/"
    + report
    + "/"
    + report
    + "_assembly_report.txt",
}
manifest = {}
for file, url in urls.items():
    target = root / file
    if not target.exists():
        for attempt in range(3):
            try:
                with (
                    urllib.request.urlopen(url, timeout=90) as r,
                    open(str(target) + ".part", "wb") as out,
                ):
                    while block := r.read(1024 * 1024):
                        out.write(block)
                Path(str(target) + ".part").replace(target)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2)
    manifest[file] = {
        "url": url,
        "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "bytes": target.stat().st_size,
    }
    print("DOWNLOADED", file, manifest[file]["bytes"], flush=True)
expected = next(
    l.split()[:2] for l in (root / "CHECKSUMS").read_text().splitlines() if l.split()[-1] == name
)
observed = subprocess.check_output(["sum", str(root / name)], text=True).split()[:2]
assert list(map(int, expected)) == list(map(int, observed)), (expected, observed)
manifest["checksum_verification"] = "PASS: Ensembl BSD sum and block count"
with gzip.open(root / name, "rt") as f:
    header = []
    for l in f:
        if not l.startswith("#"):
            break
        header.append(l.strip())
assert any("GCA_016700215.2" in x for x in header), header
assert any("GRCg7w" in x for x in header), header
manifest["header"] = header
manifest["assembly_report_note"] = (
    "NCBI directory currently names accession .2 GRCg7w_WZ. Validate sequence names and lengths against supplied FASTA before analysis."
)
(root / "download_manifest.json").write_text(json.dumps(manifest, indent=2))
print("DOWNLOAD_VERIFIED", header, flush=True)
