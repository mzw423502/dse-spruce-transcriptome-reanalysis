from pathlib import Path
import json,csv,hashlib,datetime,os,time,resource,traceback,itertools,sys
import Bio,numpy
from Bio import Align
from Bio.Align import substitution_matrices
O=Path(__file__).parent;S=O.parents[2];ENV=S/'04_environment/sequence_pairwise_20261028_1212';M=substitution_matrices.load('BLOSUM62');CAP=32
now=lambda:datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
sha=lambda b:hashlib.sha256(b).hexdigest()
def meta(p):
 p=Path(p);return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p.read_bytes())}
def fasta(p):
 d={};name=None;buf=[]
 for line in Path(p).read_text().splitlines():
  if line.startswith('>'):
   if name is not None:assert name not in d;d[name]=''.join(buf)
   name=line[1:].split()[0];buf=[]
  else:buf.append(line)
 if name is not None:assert name not in d;d[name]=''.join(buf)
 return d
def readrows(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def aligner(mode):
 a=Align.PairwiseAligner();a.mode=mode;a.substitution_matrix=M;a.open_gap_score=-11;a.extend_gap_score=-1;a.epsilon=1e-6;return a
ALIGN={k:aligner(k) for k in ['local','global']}
def analyze(h,s,mode):
 base={'mode':mode,'host_length':len(h),'source_length':len(s),'host_sha256':sha(h.encode()),'source_sha256':sha(s.encode()),'host_X':h.count('X'),'source_X':s.count('X'),'host_internal_stop':h.count('*'),'source_internal_stop':s.count('*'),'host_unsupported':''.join(sorted(set(h)-set(M.alphabet))),'source_unsupported':''.join(sorted(set(s)-set(M.alphabet)))}
 if not h or not s:return {**base,'status':'EMPTY_SEQUENCE'}
 if len(h)*len(s)>200000000:return {**base,'status':'PAIR_CELL_RESOURCE_LIMIT'}
 if base['host_unsupported'] or base['source_unsupported']:return {**base,'status':'UNSUPPORTED_RESIDUE'}
 als=ALIGN[mode].align(h,s)
 try:nt=len(als);overflow=False
 except OverflowError:nt=None;overflow=True
 candidates=list(itertools.islice(als,CAP))
 if not candidates:return {**base,'status':'NO_LOCAL_ALIGNMENT' if mode=='local' else 'NO_GLOBAL_ALIGNMENT','tie_count':nt,'tie_count_overflow':overflow}
 al=min(candidates,key=lambda x:tuple(x.coordinates.flatten().tolist()));coords=al.coordinates.tolist();hs,ss=al[0],al[1];start_h,end_h=int(coords[0][0]),int(coords[0][-1]);start_s,end_s=int(coords[1][0]),int(coords[1][-1]);assert len(hs)==len(ss) and hs.replace('-','')==h[start_h:end_h] and ss.replace('-','')==s[start_s:end_s]
 ident=paired=positive=unambiguous_ident=0;ops=[];calcscore=0.;gaprun=None;gaplen=0;hpos=start_h;spos=start_s;paired_h=[];paired_s=[]
 def flush_gap():return -11-(gaplen-1) if gaplen else 0
 for hc,sc in zip(hs,ss):
  assert not(hc=='-' and sc=='-')
  op='I' if hc=='-' else 'D' if sc=='-' else '=' if hc==sc else 'X'
  if op in ['I','D']:
   if gaprun!=op:calcscore+=flush_gap();gaplen=0;gaprun=op
   gaplen+=1
  else:
   calcscore+=flush_gap();gaprun=None;gaplen=0;calcscore+=float(M[hc,sc]);paired+=1;ident+=int(hc==sc);positive+=int(M[hc,sc]>0);unambiguous_ident+=int(hc==sc and hc in 'ACDEFGHIKLMNPQRSTVWY');paired_h.append(hpos);paired_s.append(spos)
  if ops and ops[-1][1]==op:ops[-1][0]+=1
  else:ops.append([1,op])
  hpos+=int(hc!='-');spos+=int(sc!='-')
 calcscore+=flush_gap();assert abs(calcscore-al.score)<1e-6 and hpos==end_h and spos==end_s
 nc=len(hs);pct=lambda n,d:100*n/d if d else None
 result={**base,'status':'ALIGNED','score':float(al.score),'score_recomputed':calcscore,'alignment_columns':nc,'paired_residues':paired,'identical_residues':ident,'unambiguous_identical_residues':unambiguous_ident,'positive_BLOSUM62_residues':positive,'identity_pct_columns':pct(ident,nc),'identity_pct_paired':pct(ident,paired),'positive_pct_paired':pct(positive,paired),'host_gap_columns':hs.count('-'),'source_gap_columns':ss.count('-'),'gap_runs':sum(op in ['I','D'] for n,op in ops),'host_start0':start_h,'host_end0_exclusive':end_h,'source_start0':start_s,'source_end0_exclusive':end_s,'host_span_coverage_pct':pct(end_h-start_h,len(h)),'source_span_coverage_pct':pct(end_s-start_s,len(s)),'host_paired_residue_coverage_pct':pct(paired,len(h)),'source_paired_residue_coverage_pct':pct(paired,len(s)),'host_outside_alignment_left':start_h,'host_outside_alignment_right':len(h)-end_h,'source_outside_alignment_left':start_s,'source_outside_alignment_right':len(s)-end_s,'host_unpaired_left_to_first_pair':paired_h[0] if paired_h else len(h),'host_unpaired_right_from_last_pair':len(h)-paired_h[-1]-1 if paired_h else len(h),'source_unpaired_left_to_first_pair':paired_s[0] if paired_s else len(s),'source_unpaired_right_from_last_pair':len(s)-paired_s[-1]-1 if paired_s else len(s),'coordinates':coords,'CIGAR':''.join(str(n)+op for n,op in ops),'host_aligned':hs,'source_aligned':ss,'tie_count':nt,'tie_count_overflow':overflow,'ties_considered':len(candidates),'tie_enumeration_capped':overflow or nt>CAP,'chosen_rule':'lexicographically smallest flattened coordinates among first <=32 optimal alignments in pinned library order','considered_coordinates_sha256':sha(json.dumps([x.coordinates.tolist() for x in candidates],separators=(',',':')).encode()),'algorithm':ALIGN[mode].algorithm,'reconstruction_and_score_check':'PASS'}
 return result

def synthetic_qa():
 tests=[];s='ACDEFGHIKLMNPQRSTVWY'
 for mode in ALIGN:
  r=analyze(s,s,mode);assert r['paired_residues']==20 and r['identical_residues']==20 and r['CIGAR']=='20=' and r['score']==sum(M[x,x] for x in s);tests.append(mode+'_identity')
 r=analyze('WWW'+s+'WWW',s,'local');assert r['host_start0']==3 and r['host_end0_exclusive']==23 and r['source_paired_residue_coverage_pct']==100 and r['host_unpaired_left_to_first_pair']==3;tests.append('local_subsequence_ends_coverage')
 r=analyze('WWW'+s+'WWW',s,'global');assert r['CIGAR']=='3D20=3D' and r['score']==sum(M[x,x] for x in s)-26 and r['host_unpaired_left_to_first_pair']==3 and r['host_unpaired_right_from_last_pair']==3;tests.append('global_terminal_gaps_length3_penalty13_each')
 r=analyze('ACDEFGHIKLMNPQRSTVWY','ACDEFGHINPQRSTVWY','global');assert r['CIGAR']=='8=3D9=' and r['score_recomputed']==r['score'];tests.append('internal_gap_length3_CIGAR_score_reconstruction')
 r=analyze('AAAA','AA','global');assert r['tie_count']>1 and r==analyze('AAAA','AA','global');tests.append('repeat_tie_stable_serialization')
 r=analyze('AAAA','WWWW','local');assert r['status']=='NO_LOCAL_ALIGNMENT';tests.append('no_positive_local_explicit')
 r=analyze('ACDX','ACDX','global');assert r['host_X']==1 and r['identical_residues']==4 and r['unambiguous_identical_residues']==3;tests.append('X_separate_identity')
 assert analyze('ACDU','ACDU','global')['status']=='UNSUPPORTED_RESIDUE';assert analyze('','A','local')['status']=='EMPTY_SEQUENCE';tests.append('unsupported_and_empty_explicit')
 return {'status':'PASS','tests':tests}

def run(mode):
 dest=O/('pilot100' if mode=='pilot' else 'all_fixed_pairs');dest.mkdir(exist_ok=False);started=time.monotonic();rec={'started_bjt':now(),'pid':os.getpid(),'process_start_ticks':Path('/proc/self/stat').read_text().split()[21],'script':meta(__file__),'mode':mode,'exit':None};(dest/'execution.reservation').open('x').write(json.dumps(rec))
 try:
  assert json.loads((O/'export_receipt.json').read_text())['exit']==0;qa=synthetic_qa();(dest/'synthetic_QA.json').open('x').write(json.dumps(qa,indent=2))
  inputfile=O/('pilot100_fixed_pairs.tsv' if mode=='pilot' else 'unique_sequence_pair_manifest.tsv');pairs=readrows(inputfile);hosts=fasta(O/'host_units_current_peptides.faa');sources=fasta(O/'source_accession_version_proteins.faa')
  ins=[meta(p) for p in [inputfile,O/'host_units_current_peptides.faa',O/'source_accession_version_proteins.faa',ENV/'functional_QA.json',ENV/'installation_receipt.json']]
  protocol={'frozen_before_alignments_bjt':now(),'input_files':ins,'scope':'Current existing legacy candidate pairs only, sequence evidence audit, no database search, no orthology/function equivalence or GO transfer.','Bio':Bio.__version__,'numpy':numpy.__version__,'python':sys.version,'matrix':'BLOSUM62','matrix_text_sha256':sha(str(M).encode()),'modes':['local','global'],'gap_open':-11,'gap_extend':-1,'gap_formula':'open+(n-1)*extend for n>=1; terminal gaps use same scores as internal; all insertion/deletion sides identical.','epsilon':1e-6,'algorithms':{k:a.algorithm for k,a in ALIGN.items()},'coordinates':'0-based half-open target=host, query=source','CIGAR':'= identical chars; X unequal nongap; D host residue/source gap; I host gap/source residue; no soft-clips in CIGAR, outside alignment recorded separately.','metrics':'Identity numerator literal char equality; separate canonical20-aa identity count excludes BZX*. Columns denominator includes gaps; paired denominator excludes gaps. Positive=BLOSUM62>0 among paired residues. Residue coverage=paired residues/full sequence length. Span coverage=(end-start)/full length and can be100% in global despite large gaps; must not substitute span for residue coverage.','terminal_metrics':'outside_alignment_left/right counts local excluded residues. unpaired_left_to_first_pair/right_from_last_pair counts all residues before first or after last actual paired residue, including terminal global gaps.','ties':'len(PairwiseAlignments) exact unless OverflowError; enumerate <=32 optimal outputs; choose lexicographically smallest flattened coordinates in that subset; capped is not canonical among all alignments and does not establish unique biological correspondence.','errors':'Empty or unsupported alphabet explicit status, no residue replacement. No positive local ->NO_LOCAL_ALIGNMENT. Each ALIGNED row must reconstruct both source substrings and hand-recompute gap/matrix score.','resources':{'threads':1,'max_memory_GiB':3,'per_pair_max_DP_cells':200000000,'full_phase_max_wall_seconds':3600},'functional_memberships_created':0}
  (dest/'protocol.json').open('x').write(json.dumps(protocol,indent=2));status={};aligned=0
  with (dest/'alignments.jsonl').open('x') as out:
   for idx,p in enumerate(pairs):
    h=hosts[p['unit_id']];s=sources[p['source_accession_version']];assert sha(h.encode())==p['host_sha256'] and sha(s.encode())==p['source_sha256'];assert len(h)==int(p['host_length']) and len(s)==int(p['source_length'])
    for algmode in ALIGN:
     result=analyze(h,s,algmode);row={'sequence_pair_id':p['sequence_pair_id'],'unit_ids':p.get('unit_ids',p['unit_id']),'unit_source_pair_ids':p.get('unit_source_pair_ids',p['pair_id']),'source_accession_version':p['source_accession_version'],**result};out.write(json.dumps(row,separators=(',',':'))+'\n');status[result['status']]=status.get(result['status'],0)+1;aligned+=int(result['status']=='ALIGNED')
    if idx%100==0 or idx+1==len(pairs):
     out.flush();tmp=dest/'progress.json.tmp';tmp.write_text(json.dumps({'beijing_time':now(),'pid':os.getpid(),'pairs_done':idx+1,'pairs_total':len(pairs),'status_counts':status,'wall_seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}));tmp.replace(dest/'progress.json')
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<3*1024**2,'MEMORY_LIMIT';assert time.monotonic()-started<3600,'WALL_LIMIT'
  for m in ins:assert meta(m['path'])['sha256']==m['sha256']
  rec.update({'status':'COMPLETED_SEQUENCE_AUDIT_NOT_FUNCTION_TRANSFER','exit':0,'pairs':len(pairs),'mode_rows':2*len(pairs),'status_counts':status,'input_hashes_unchanged':True,'output':meta(dest/'alignments.jsonl'),'protocol':meta(dest/'protocol.json')})
 except Exception as e:rec.update({'status':'FAILED','exit':1,'error':repr(e),'traceback':traceback.format_exc()})
 rec.update({'finished_bjt':now(),'wall_seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss});(dest/'execution_receipt.json').open('x').write(json.dumps(rec,indent=2));(dest/'exit.code').write_text(str(rec['exit'])+'\n');print(json.dumps(rec,indent=2));return rec['exit']
if __name__=='__main__':
 if sys.argv[1]=='qa':print(json.dumps(synthetic_qa()));sys.exit(0)
 assert sys.argv[1] in ['pilot','full'];sys.exit(run(sys.argv[1]))
