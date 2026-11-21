"""Read-only extraction of frozen functional-source and response evidence.

This portable entry point performs the same join as
review/extract_functional_source_audit.py. It does not calculate new
statistics. It joins already frozen tables to make the provenance of the 39
source-linked direction-opposite units auditable.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("frozen_base", type=Path)
parser.add_argument("output_reviewdir", type=Path)
args = parser.parse_args()

ROOT = args.frozen_base
OUT = args.output_reviewdir
OUT.mkdir(parents=True, exist_ok=True)
LEDGER = ROOT / "functional_source_ledger"
NUM = ROOT / "numerical"

def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

opposite = read_tsv(NUM / "annotation_response_opposite.tsv")
unit_ledger = {r["unit_id"]: r for r in read_tsv(LEDGER / "all64238_unit_evidence_ledger.tsv")}
selected = [
    row for row in opposite
    if int(unit_ledger[row["unit_id"]]["source_linked_hypothesis_pairs_n"]) > 0
]
selected_ids = {r["unit_id"] for r in selected}
chains = [
    r for r in read_tsv(LEDGER / "all27930_MA_source_chain_ledger.tsv")
    if r["unit_id"] in selected_ids
]
chain_by_pair = {(r["unit_id"], r["source_accession_version"]): r for r in chains}
accessions = {r["source_accession_version"] for r in chains}
source_records = {
    r["source_accession_version"]: r
    for r in read_tsv(LEDGER / "all3908_source_record_summary.tsv")
    if r["source_accession_version"] in accessions
}
relations = [
    r for r in read_tsv(LEDGER / "source_GO_relation_ledger.tsv")
    if r["source_accession_version"] in accessions
]

# GO annotations are entry keyed, while relations provide required source entries.
ann_by_entry = defaultdict(list)
with (LEDGER / "normalized_source_GO_annotations.jsonl").open(encoding="utf-8") as handle:
    for line in handle:
        row = json.loads(line)
        ann_by_entry[row["entry"]].append(row)

rel_by_acc = defaultdict(list)
for rel in relations:
    rel_by_acc[rel["source_accession_version"]].append(rel)

def evidence_for_relation(rel):
    entries = [x for x in rel["required_entries"].split(";") if x]
    support = []
    for entry in entries:
        for ann in ann_by_entry.get(entry, []):
            if ann["go_id"] == rel["GO"] and not ann["is_NOT"]:
                support.append({
                    "entry": entry,
                    "GO_evidence": ann["GO_evidence"],
                    "reference": ann["reference"],
                    "annotation_id": ann["annotation_id"],
                    "taxon_id": ann["taxon_id"],
                    "term_name": ann["term_name"],
                })
    return support

rows = []
for op in selected:
    unit = unit_ledger[op["unit_id"]]
    for acc in unit["source_accessions"].split(";"):
        if not acc:
            continue
        chain = chain_by_pair.get((op["unit_id"], acc), {})
        record = source_records.get(acc, {})
        acc_relations = rel_by_acc.get(acc, [])
        if not acc_relations:
            rows.append({
                "unit_id": op["unit_id"], "representative_gene": op["representative_gene"],
                "host_annotation_description": op["annotation_descriptions"],
                "source_accession": acc, "source_record_state": record.get("state", "MISSING"),
                "source_GO_chain_state": chain.get("source_GO_candidate_state", "MISSING"),
                "source_GO": "", "source_GO_name": "", "GO_evidence": "", "GO_reference": "",
                "local_identity_pct": chain.get("local_identity_pct_paired", ""),
                "local_host_coverage_pct": chain.get("local_host_paired_residue_coverage_pct", ""),
                "local_source_coverage_pct": chain.get("local_source_paired_residue_coverage_pct", ""),
                "P0_log2FC": op["log2FoldChange_P0"], "P0_q": op["padj_simple_four_fixed_P0"],
                "P1_log2FC": op["log2FoldChange_P1"], "P1_q": op["padj_simple_four_fixed_P1"],
                "I_log2FC": op["log2FoldChange_I"], "I_q": op["padj_I"],
                "audit_note": "Source-linked unit has no retained direct GO relation row."
            })
        for rel in acc_relations:
            supports = evidence_for_relation(rel)
            evidence = ";".join(sorted({x["GO_evidence"] for x in supports}))
            refs = ";".join(sorted({x["reference"] for x in supports}))
            rows.append({
                "unit_id": op["unit_id"], "representative_gene": op["representative_gene"],
                "host_annotation_description": op["annotation_descriptions"],
                "source_accession": acc, "source_record_state": record.get("state", "MISSING"),
                "source_GO_chain_state": chain.get("source_GO_candidate_state", "MISSING"),
                "source_GO": rel["GO"], "source_GO_name": ";".join(sorted({x["term_name"] for x in supports})),
                "GO_evidence": evidence, "GO_reference": refs,
                "local_identity_pct": chain.get("local_identity_pct_paired", ""),
                "local_host_coverage_pct": chain.get("local_host_paired_residue_coverage_pct", ""),
                "local_source_coverage_pct": chain.get("local_source_paired_residue_coverage_pct", ""),
                "P0_log2FC": op["log2FoldChange_P0"], "P0_q": op["padj_simple_four_fixed_P0"],
                "P1_log2FC": op["log2FoldChange_P1"], "P1_q": op["padj_simple_four_fixed_P1"],
                "I_log2FC": op["log2FoldChange_I"], "I_q": op["padj_I"],
                "audit_note": "GO applies to required source entry/entries; it is not a host GO membership."
            })

fieldnames = list(rows[0])
with (OUT / "opposite39_source_function_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader(); writer.writerows(rows)

summary = {
    "scope": "Read-only join of frozen numerical and functional-source-ledger files; no new statistics or alignment.",
    "opposite_nonrisk_units": len(opposite),
    "source_linked_hypothesis_units_within_opposite": len(selected),
    "source_accession_pairs": len(chains),
    "direct_source_GO_relation_rows": len(rows),
    "source_record_states": Counter(r["source_record_state"] for r in rows),
    "GO_evidence_codes": Counter(code for r in rows for code in r["GO_evidence"].split(";") if code),
    "source_chain_states": Counter(r["source_GO_chain_state"] for r in rows),
    "input_files": [
        "numerical/annotation_response_opposite.tsv",
        "functional_source_ledger/all64238_unit_evidence_ledger.tsv",
        "functional_source_ledger/all27930_MA_source_chain_ledger.tsv",
        "functional_source_ledger/all3908_source_record_summary.tsv",
        "functional_source_ledger/source_GO_relation_ledger.tsv",
        "functional_source_ledger/normalized_source_GO_annotations.jsonl"
    ],
    "hard_boundary": "The ledger explicitly creates zero host_GO_memberships and zero module_memberships. Source-linked GO is evidence about a required source protein entry only, not direct annotation, orthology, or validated function of the spruce sequence-equivalence unit."
}
(OUT / "opposite39_source_function_audit_summary.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(summary, ensure_ascii=False, indent=2))
