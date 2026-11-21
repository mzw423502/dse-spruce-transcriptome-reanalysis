#!/usr/bin/env python3
"""Read-only descriptive summary for the fixed original-391 selection.

No model fitting is performed. Inputs are derived TSVs produced by the approved
classification extension. The resulting q-status counts are exploratory,
conditional on the original selection, and are not equivalence tests.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import sys


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def finite_float(x):
    try:
        value = float(x)
    except (TypeError, ValueError):
        return math.nan
    return value if math.isfinite(value) else math.nan


def support_ids(path):
    return {
        row["feature_id"]
        for row in read_tsv(path)
        if row.get("retained_supported_original_pattern") == "TRUE"
    }


def summarize(ids, label, p0, p1, d0, d1):
    rows = []
    strata = [
        ("all", ids),
        ("P0_positive_P1_negative", {fid for fid in ids if p0[fid] > 0 and p1[fid] < 0}),
        ("P0_negative_P1_positive", {fid for fid in ids if p0[fid] < 0 and p1[fid] > 0}),
    ]
    for stratum, subset in strata:
        q0 = {fid: d0[fid] for fid in subset}
        q1 = {fid: d1[fid] for fid in subset}
        valid = {fid for fid in subset if math.isfinite(q0[fid]) and math.isfinite(q1[fid])}
        d0sig = {fid for fid in valid if q0[fid] < 0.05}
        d1sig = {fid for fid in valid if q1[fid] < 0.05}
        rows.append({
            "set": label,
            "P0_P1_direction_stratum": stratum,
            "n": len(subset),
            "q_pair_finite": len(valid),
            "q_pair_nonfinite": len(subset - valid),
            "D0_q_lt_0_05": len(d0sig),
            "D1_q_lt_0_05": len(d1sig),
            "both_D0_D1_q_lt_0_05": len(d0sig & d1sig),
            "D0_only_q_lt_0_05": len(d0sig - d1sig),
            "D1_only_q_lt_0_05": len(d1sig - d0sig),
            "neither_D0_D1_q_lt_0_05": len(valid - (d0sig | d1sig)),
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True)
    args = parser.parse_args()
    directory = args.analysis_dir
    files = {
        "five_contrasts": os.path.join(directory, "primary391_five_contrasts.tsv"),
        "allN": os.path.join(directory, "primary391_allN_classification.tsv"),
        "competition": os.path.join(directory, "primary391_competition_classification.tsv"),
        "edgeR_QL": os.path.join(directory, "primary391_edgeR_QL_classification.tsv"),
    }
    for path in files.values():
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    source = read_tsv(files["five_contrasts"])
    by_contrast = {}
    for row in source:
        by_contrast.setdefault(row["contrast"], {})[row["feature_id"]] = row
    required = ["DSE_at_pathogen0", "DSE_at_pathogen1", "pathogen_at_DSE0", "pathogen_at_DSE1", "interaction"]
    if set(required) - set(by_contrast):
        raise ValueError("missing contrasts: " + repr(set(required) - set(by_contrast)))
    primary = set(by_contrast["interaction"])
    if len(primary) != 391 or any(set(by_contrast[c]) != primary for c in required):
        raise ValueError("five-contrast table is not a complete 391-unit rectangle")
    triple = support_ids(files["allN"]) & support_ids(files["competition"]) & support_ids(files["edgeR_QL"])
    if len(triple) != 127 or not triple <= primary:
        raise ValueError("unexpected three-route support intersection")
    p0 = {fid: finite_float(by_contrast["pathogen_at_DSE0"][fid]["log2FoldChange"]) for fid in primary}
    p1 = {fid: finite_float(by_contrast["pathogen_at_DSE1"][fid]["log2FoldChange"]) for fid in primary}
    d0 = {fid: finite_float(by_contrast["DSE_at_pathogen0"][fid]["padj_simple_four_fixed"]) for fid in primary}
    d1 = {fid: finite_float(by_contrast["DSE_at_pathogen1"][fid]["padj_simple_four_fixed"]) for fid in primary}
    if any(not math.isfinite(p0[fid]) or not math.isfinite(p1[fid]) for fid in primary):
        raise ValueError("nonfinite pathogen effects in primary 391")
    rows = summarize(primary, "original_DESeq2_391", p0, p1, d0, d1)
    rows.extend(summarize(triple, "three_route_supported_127", p0, p1, d0, d1))
    out_tsv = os.path.join(directory, "baseline_context_summary.tsv")
    columns = list(rows[0])
    with open(out_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    script_sha = sha256(os.path.realpath(__file__))
    report = {
        "status": "COMPLETE",
        "scope": "Read-only descriptive counts from the existing primary391 five-contrast export and three existing route-classification exports; no model refit.",
        "exploratory_status": "Post hoc descriptive context only; conditional on the original 391-unit DESeq2 selection and the later route-specific support intersection. It is not an independently designed validation analysis.",
        "definitions": {
            "D0": "DSE_at_pathogen0, evaluated with the already exported fixed-four-simple-effect BH q value.",
            "D1": "DSE_at_pathogen1, evaluated with the already exported fixed-four-simple-effect BH q value.",
            "direction_strata": "P0 and P1 use the already exported DESeq2 pathogen simple-effect coefficients.",
            "three_route_supported_127": "Intersection of allN and competition saved original-DESeq2 support with edgeR's distinct QL-analogue support; it does not mean the classifiers are identical.",
            "non_significance": "Neither D0 nor D1 reaching q<0.05 is not evidence of equivalence or absence of a baseline DSE effect."
        },
        "inputs": [{"path": path, "sha256": sha256(path)} for path in files.values()],
        "script": {"path": os.path.realpath(__file__), "sha256": script_sha},
        "rows": rows,
    }
    with open(os.path.join(directory, "baseline_context_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(json.dumps({"status": "COMPLETE", "rows": rows}, ensure_ascii=False))


if __name__ == "__main__":
    main()
