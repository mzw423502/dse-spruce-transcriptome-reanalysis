#!/usr/bin/env python3
"""Portable rerun of the six exploratory, outside-GA diagnostic records."""
import argparse, csv, hashlib, json, sys
from pathlib import Path

TARGETS = {
    "PF00847": ["MA_81979g0010|host|PaEq_960a1016a6ff6475", "MA_81979g0010|frozen_source|XP_016682358.1"],
    "PF02798": ["MA_10426104g0010|host|PaEq_6c314d0bda685b5e", "MA_10426104g0010|frozen_source|XP_010923805.1"],
    "PF13417": ["MA_10426104g0010|host|PaEq_6c314d0bda685b5e", "MA_10426104g0010|frozen_source|XP_010923805.1"],
}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

parser = argparse.ArgumentParser()
parser.add_argument("--root", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()
root, output = Path(args.root), Path(args.output)
output.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(root / "tools/pyhmmer"))
import pyhmmer

alphabet = pyhmmer.easel.Alphabet.amino()
with pyhmmer.easel.SequenceFile(root / "inputs/candidate_host_and_frozen_source_peptides.faa", digital=True, alphabet=alphabet) as sf:
    sequences = sf.read_block()
rows = []
for acc, wanted in TARGETS.items():
    hmm = next(pyhmmer.plan7.HMMFile(root / "models" / (acc + ".hmm")))
    top = next(pyhmmer.hmmer.hmmsearch([hmm], sequences, cpus=1, E=1000, domE=1000, incE=1000, incdomE=1000))
    hits = {str(hit.name): hit for hit in top.reported}
    for name in wanted:
        hit = hits.get(name)
        row = {"profile_accession": str(hmm.accession), "profile_name": str(hmm.name), "sequence_header": name,
               "sequence_ga_score": hmm.cutoffs.gathering[0], "domain_ga_score": hmm.cutoffs.gathering[1],
               "diagnostic_search_E": 1000, "diagnostic_domain_E": 1000, "formal_GA_call": "not_included"}
        if hit:
            domains = list(hit.domains.reported)
            row.update({"full_sequence_score": hit.score, "full_sequence_evalue": hit.evalue,
                        "reported_domain_count": len(domains), "included_domain_count": len(hit.domains.included)})
            if domains:
                dom = domains[0]
                ok = hit.score >= hmm.cutoffs.gathering[0] and dom.score >= hmm.cutoffs.gathering[1]
                row.update({"formal_GA_call": "included" if ok else "not_included", "best_domain_score": dom.score,
                            "best_domain_i_evalue": dom.i_evalue, "hmm_from": dom.alignment.hmm_from,
                            "hmm_to": dom.alignment.hmm_to, "ali_from": dom.alignment.target_from,
                            "ali_to": dom.alignment.target_to})
        for key in ("full_sequence_score", "full_sequence_evalue", "reported_domain_count", "included_domain_count", "best_domain_score", "best_domain_i_evalue", "hmm_from", "hmm_to", "ali_from", "ali_to"):
            row.setdefault(key, "")
        rows.append(row)
fields = list(rows[0])
table = output / "rerun_ga_outside_diagnostic.tsv"
with table.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader(); writer.writerows(rows)
receipt = {"status": "COMPLETE_EXPLORATORY_DIAGNOSTIC_RERUN", "rows": len(rows), "cpus": 1,
           "parameters": {"E": 1000, "domE": 1000, "incE": 1000, "incdomE": 1000},
           "output_sha256": sha256(table), "limit": "Exploratory diagnostic only; it does not alter formal Pfam-GA classification."}
(output / "rerun_diagnostic_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt, sort_keys=True))
