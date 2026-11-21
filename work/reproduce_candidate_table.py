"""Read-only cross-route and descriptive-expression extraction for 39 candidates.

This portable entry point performs the same extraction as
review/build_candidate_crossroute_audit.py. No fit, hypothesis test, or
reclassification is run here. Group means merely summarize the already
supplied 12-pool normalized-count file and must not be substituted for model
estimates.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("analysis_dir", type=Path)
parser.add_argument("source_audit_tsv", type=Path)
parser.add_argument("output_tsv", type=Path)
args = parser.parse_args()

AUDIT = args.source_audit_tsv
ANALYSIS = args.analysis_dir
OUT = args.output_tsv
OUT.parent.mkdir(parents=True, exist_ok=True)

def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

source_first = {}
for row in read_tsv(AUDIT):
    source_first.setdefault(row["unit_id"], row)
routes = {}
for route in ("allN", "competition", "edgeR_QL"):
    routes[route] = {
        row["feature_id"]: row
        for row in read_tsv(ANALYSIS / f"primary391_{route}_classification.tsv")
    }
contrasts = defaultdict(dict)
for row in read_tsv(ANALYSIS / "primary391_five_contrasts.tsv"):
    if row["feature_id"] in source_first:
        contrasts[row["feature_id"]][row["contrast"]] = row
values = defaultdict(lambda: defaultdict(list))
for row in read_tsv(ANALYSIS / "primary391_expression_12pools.tsv"):
    if row["feature_id"] in source_first:
        values[row["feature_id"]][row["group"]].append(float(row["normalized_count"]))

out_rows=[]
for unit_id, src in sorted(source_first.items()):
    r = dict(src)
    for route, tbl in routes.items():
        x=tbl[unit_id]
        r[f"{route}_status"] = x["status"]
        r[f"{route}_full_support"] = x["retained_supported_original_pattern"]
    for short, contrast in {
        "D0":"DSE_at_pathogen0", "D1":"DSE_at_pathogen1",
        "P0":"pathogen_at_DSE0", "P1":"pathogen_at_DSE1", "I":"interaction"
    }.items():
        x=contrasts[unit_id][contrast]
        r[f"{short}_log2FC_frozen"] = x["log2FoldChange"]
        r[f"{short}_q_frozen"] = x["padj_simple_four_fixed"] if short != "I" else x["padj_contrast_fixed"]
    for group in ("Pa", "PaPs", "PaHp", "PaPsHp"):
        x=values[unit_id][group]
        r[f"{group}_norm_count_mean_descriptive"] = f"{sum(x)/len(x):.6f}"
        r[f"{group}_norm_count_3pools"] = ";".join(f"{v:.6f}" for v in x)
    r["expression_audit_note"] = "Descriptive normalized counts from three pooled biological replicates per group; group means are not model contrasts."
    out_rows.append(r)

with OUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(out_rows[0]), delimiter="\t")
    writer.writeheader(); writer.writerows(out_rows)
print(f"wrote {len(out_rows)} candidates to {OUT}")
