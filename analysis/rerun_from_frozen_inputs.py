#!/usr/bin/env python3
import argparse,csv,hashlib,sys
from pathlib import Path
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();root=Path(a.root);out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(root/'tools/pyhmmer'));import pyhmmer
models=[]
for r in csv.DictReader((root/'models/model_manifest.tsv').open(),delimiter='\t'):
 h=next(pyhmmer.plan7.HMMFile(root/'models'/(r['accession']+'.hmm')));models.append(h)
alpha=pyhmmer.easel.Alphabet.amino()
with pyhmmer.easel.SequenceFile(root/'inputs/candidate_host_and_frozen_source_peptides.faa',digital=True,alphabet=alpha) as sf: seqs=sf.read_block()
rows=[]
for top in pyhmmer.hmmer.hmmsearch(models,seqs,cpus=1,bit_cutoffs='gathering'):
 for hit in top:
  for idx,d in enumerate(hit.domains.included,1):
   sg,dg=top.query.cutoffs.gathering
   rows.append({'profile_accession':str(top.query.accession),'profile_name':str(top.query.name),'sequence_header':str(hit.name),'protein_length':hit.length,'full_sequence_score':hit.score,'sequence_ga_score':sg,'domain_score':d.score,'domain_ga_score':dg,'full_sequence_evalue':hit.evalue,'independent_domain_evalue':d.i_evalue,'hmm_from':d.alignment.hmm_from,'hmm_to':d.alignment.hmm_to,'ali_from':d.alignment.target_from,'ali_to':d.alignment.target_to,'hmm_coverage_fraction':round((d.alignment.hmm_to-d.alignment.hmm_from+1)/top.query.M,6),'protein_coverage_fraction':round((d.alignment.target_to-d.alignment.target_from+1)/hit.length,6),'ga_pass':str(hit.score>=sg and d.score>=dg).lower()})
fields=list(rows[0]);p=out/'rerun_domain_hits.tsv'
with p.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
print(len(rows),sha(p))
