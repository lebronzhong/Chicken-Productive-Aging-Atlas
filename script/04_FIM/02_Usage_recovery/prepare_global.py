#!/usr/bin/env python3
# Chicken Productive Aging Atlas | 04_FIM
# Purpose: Prepare global.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, shutil, gc, re
from scipy import sparse

r = Path(__file__).resolve().parents[1]
b = r.parent
old = b / "count_DTU"
psiold = b / "isoform_PSI"
mp = pd.read_csv(old / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t")
mp.to_csv(r / "inputs/canonical_transcript_mapping.tsv.gz", sep="\t", index=False)
src = {x: y.transcript_id for y in mp.itertuples() for x in y.source_ids.split(",")}
pos = {x: i for i, x in enumerate(mp.transcript_id)}
par = mp.set_index("transcript_id").gene_id
srcset = set(src)
meta = pd.read_csv(psiold / "inputs/sample_manifest.tsv", sep="\t")
tissues = sorted(meta.tissue.unique())
manifest = []
audit = []


def readfull(p):
    d = pd.read_csv(p, sep="\t", index_col=0)
    assert d.index.is_unique
    audit.append(
        dict(
            path=str(p),
            rows=len(d),
            columns=len(d.columns),
            unmapped_rows=int((~d.index.isin(src)).sum()),
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        )
    )
    d = d[d.index.isin(src)]
    present = set(d.index)
    assert set(src) <= present, (str(p), len(set(src) - present))
    d.index = d.index.map(src)
    d = d.groupby(level=0, sort=False).sum().loc[mp.transcript_id]
    assert np.isfinite(d.to_numpy()).all()
    assert d.to_numpy().min() >= 0
    return d


if not (r / "inputs/COUNTS_READY").exists():
    m = readfull(b / "12.20240919_TPM.count/count/M_merged_transcript_grouped_counts.tsv")
    control = readfull(b / "12.20240919_TPM.count/count/100wk_merged_transcript_grouped_counts.tsv")
    validations = []
    for t in tissues:
        ss = meta[meta.tissue.eq(t)].copy()
        ss["group"] = ss.age.astype(str)
        ss["animal_id"] = ss.individual
        ss["sample_id"] = t + "::" + ss["sample"]
        ss["treatment_stage"] = "natural"
        ss["identity_basis"] = "Sample manifest used for natural-age isoform-usage analysis"
        mcols = [f"{t}M-{i}" for i in [1, 3, 4, 5, 6]]
        assert set(mcols) <= set(m.columns)
        more = pd.DataFrame(
            [
                dict(
                    tissue=t,
                    sample=f"M{i}",
                    source_sample=f"{t}M-{i}",
                    psi_sample=f"{t}M-{i}",
                    age=100,
                    individual=f"M_{i}",
                    animal_id=f"M_{i}",
                    sample_id=f"{t}::M{i}",
                    batch="UNKNOWN_FIM",
                    group="FIM",
                    treatment_stage="post_fasting_resumed_laying",
                    identity_basis="FIM sample labels and IsoQuant input list; cross-tissue animal identity is not independently assessed by this script",
                )
                for i in [1, 3, 4, 5, 6]
            ]
        )
        sm = pd.concat([ss, more], ignore_index=True)
        manifest.append(sm)
        out = r / "inputs" / t
        out.mkdir(exist_ok=True)
        x = pd.read_csv(
            old / f"inputs/validated_fig3/{t}/tx_counts.tsv.gz", sep="\t", index_col=0
        ).loc[mp.transcript_id, ss["sample"]]
        cc = control[[f"{t}100-{i}" for i in [1, 2, 3, 4, 5]]].to_numpy()
        difference = np.max(np.abs(x[ss.loc[ss.age.eq(100), "sample"]].to_numpy() - cc))
        assert difference < 1e-6, (t, difference)
        xx = np.concatenate([x.to_numpy(), m[mcols].to_numpy()], axis=1)
        assert np.max(abs(xx - np.round(xx))) < 1e-5
        np.savez_compressed(out / "counts.npz", counts=np.round(xx).astype(np.int32))
        sm.to_csv(out / "samples.tsv", sep="\t", index=False)
        validations.append(
            dict(
                tissue=t,
                canonical_transcripts=len(mp),
                control_max_abs_difference=float(difference),
                M_samples=5,
                M_count_sum=int(xx[:, -5:].sum()),
            )
        )
        print("COUNTS_VALIDATED", t, flush=True)
    pd.concat(manifest).to_csv(r / "inputs/sample_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(validations).to_csv(r / "audit/full_count_validation.tsv", sep="\t", index=False)
    (r / "inputs/COUNTS_READY").write_text("SUCCESS")
    del m, control, x, xx
    gc.collect()
if not (r / "inputs/TPM_READY").exists():
    p = b / "12.20240919_TPM.count/TPM/final_merged_transcript_grouped_tpm.tsv"
    d = readfull(p)
    for t in tissues:
        sm = pd.read_csv(r / f"inputs/{t}/samples.tsv", sep="\t")
        cols = []
        for row in sm.itertuples():
            opts = [v for v in [row.psi_sample, row.source_sample] if v in d.columns]
            assert opts, (t, row.sample)
            cols.append(opts[0])
        assert len(set(cols)) == 25
        np.savez_compressed(r / f"inputs/{t}/tpm.npz", TPM=d[cols].to_numpy())
        sm["TPM_source_sample"] = cols
        sm.to_csv(r / f"inputs/{t}/samples.tsv", sep="\t", index=False)
        print("TPM_PREPARED", t, flush=True)
    (r / "inputs/TPM_READY").write_text("SUCCESS")
    del d
    gc.collect()
pd.DataFrame(audit).to_csv(r / "audit/global_matrix_provenance.tsv", sep="\t", index=False)
if not (r / "inputs/EVENTS_READY").exists():
    p = b / "11_20240903.SUPPA/rename_Swiss/9.9_finalgtf_annotated_rename_AS.events.ioe"
    e = pd.read_csv(p, sep="\t")
    e = e[e.event_id.ne("event_id")].drop_duplicates()
    rows = []
    for z in e.itertuples():
        inc = set(z.alternative_transcripts.split(","))
        tot = set(z.total_transcripts.split(","))
        reason = "valid"
        if not inc <= tot:
            reason = "inclusion_not_subset"
        if not tot <= srcset:
            reason = "unmapped_denominator_member"
        ii = {src[x] for x in inc if x in src}
        tt = {src[x] for x in tot if x in src}
        excluded = {src[x] for x in tot - inc if x in src}
        if ii & excluded:
            reason = "canonical_inclusion_exclusion_conflict"
        parents = {par[x] for x in tt}
        if len(parents) != 1 or (parents and next(iter(parents)) != z.gene_id):
            reason = "parent_mapping_conflict"
        coord = z.event_id.split(";", 1)[-1]
        rows.append(
            dict(
                event_id=coord,
                source_event_id=z.event_id,
                gene_id=z.gene_id,
                event_type=coord.split(":")[0],
                coordinates=coord,
                scope="terminal_exon" if coord.startswith(("AF:", "AL:")) else "internal_splicing",
                inclusion_transcripts=",".join(sorted(ii)),
                total_transcripts=",".join(sorted(tt)),
                definition_status=reason,
            )
        )
    ed = pd.DataFrame(rows)
    ed["signature"] = ed.inclusion_transcripts + "|" + ed.total_transcripts + "|" + ed.gene_id
    conflicts = set(ed.groupby("event_id").signature.nunique().loc[lambda x: x > 1].index)
    ed.loc[ed.event_id.isin(conflicts), "definition_status"] = (
        "conflicting_same_coordinate_definition"
    )
    ed.to_csv(r / "inputs/event_definition_audit.tsv.gz", sep="\t", index=False)
    ed = (
        ed[ed.definition_status.eq("valid")]
        .drop_duplicates("event_id")
        .drop(columns="signature")
        .reset_index(drop=True)
    )
    for field, name in [("inclusion_transcripts", "inclusion"), ("total_transcripts", "total")]:
        ri = []
        ci = []
        for i, x in enumerate(ed[field]):
            ids = x.split(",")
            ri.extend([i] * len(ids))
            ci.extend(pos[t] for t in ids)
        mat = sparse.csr_matrix(
            (np.ones(len(ri), dtype=np.float64), (ri, ci)), shape=(len(ed), len(mp))
        )
        sparse.save_npz(r / f"inputs/event_{name}_matrix.npz", mat)
    ed.to_csv(r / "inputs/event_definitions.tsv.gz", sep="\t", index=False)
    (r / "inputs/EVENTS_READY").write_text("SUCCESS")
    print("EVENTS_DEFINED", len(ed), "of", len(rows), flush=True)
(r / "inputs/GLOBAL_READY").write_text("SUCCESS")
