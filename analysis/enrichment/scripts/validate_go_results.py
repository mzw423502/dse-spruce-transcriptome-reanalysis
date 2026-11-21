#!/usr/bin/env python3
"""Independent numerical validation of frozen GO-BP enrichment outputs.

This validator never changes the frozen R results.  It reconstructs the stated
universe and tests directly from the input tables, checks all 58 Fisher tables
and BH values, obtains deterministic conditional-hypergeometric p-values for
the same 5x5 strata, and exports transparent membership/context tables.
"""
import argparse, csv, hashlib, json, math, os, sys
from collections import defaultdict
import numpy as np

ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
CANDIDATES = {
    "MA_81979g0010": "ERF017-like",
    "MA_71687g0010": "NHL13-like",
    "MA_595526g0010": "CSLD3-like",
    "MA_10426104g0010": "GST-like",
    "MA_10432086g0010": "XTH23-like",
    "MA_5735g0010": "CHS4-like",
}

def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def write_tsv(path, rows, fields):
    with open(path, "w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def flag(x):
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}

def bh(p):
    n = len(p); order = sorted(range(n), key=lambda i: p[i]); out = [0.0] * n; run = 1.0
    for rank0 in range(n - 1, -1, -1):
        i = order[rank0]; run = min(run, p[i] * n / (rank0 + 1)); out[i] = run
    return out

def logchoose(n, k):
    if k < 0 or k > n: return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def hg_logpmf(x, N, K, n):
    return logchoose(K, x) + logchoose(N - K, n - x) - logchoose(N, n)

def fisher_one_greater(a, b, c, d):
    N, K, n = a + b + c + d, a + c, a + b
    lo, hi = max(0, n - (N - K)), min(n, K)
    lps = [hg_logpmf(x, N, K, n) for x in range(lo, hi + 1)]
    mx = max(lps); probs = [math.exp(z - mx) for z in lps]
    return min(1.0, math.exp(mx) * sum(v for x, v in zip(range(lo, hi + 1), probs) if x >= a))

def fisher_two_sided(a, b, c, d):
    N, K, n = a + b + c + d, a + c, a + b
    lo, hi = max(0, n - (N - K)), min(n, K)
    lps = [hg_logpmf(x, N, K, n) for x in range(lo, hi + 1)]
    lpobs = hg_logpmf(a, N, K, n)
    mx = max(lps)
    return min(1.0, math.exp(mx) * sum(math.exp(z - mx) for z in lps if z <= lpobs + 1e-12))

def rank_quintile(x):
    n = len(x); ans = [None] * n
    order = sorted(range(n), key=lambda i: x[i]); left = 0
    while left < n:
        right = left
        while right + 1 < n and x[order[right + 1]] == x[order[left]]: right += 1
        rank = ((left + 1) + (right + 1)) / 2.0
        q = min(5, int(math.floor(rank / (n + 1) * 5)) + 1)
        for z in range(left, right + 1): ans[order[z]] = q
        left = right + 1
    return ans

def hypergeom_vector(N, K, n):
    lo, hi = max(0, n - (N - K)), min(n, K)
    x = np.arange(lo, hi + 1)
    lp = np.array([hg_logpmf(int(v), N, K, n) for v in x])
    pr = np.exp(lp - np.max(lp)); pr /= pr.sum()
    return lo, pr

def convolved_counts(strata):
    """Return exact pmf of summed hypergeometrics as a vector indexed by count."""
    dist = np.array([1.0])
    offset = 0
    for N, K, n in strata:
        lo, pr = hypergeom_vector(N, K, n)
        dist = np.convolve(dist, pr); offset += lo
    dist /= dist.sum()
    return offset, dist

def descriptive(values):
    vals = sorted(values)
    if not vals: return {"n": 0, "min": None, "q25": None, "median": None, "mean": None, "q75": None, "max": None}
    def q(p):
        k = (len(vals) - 1) * p; lo = int(math.floor(k)); hi = int(math.ceil(k))
        return vals[lo] if lo == hi else vals[lo] + (vals[hi] - vals[lo]) * (k - lo)
    return {"n": len(vals), "min": vals[0], "q25": q(.25), "median": q(.5), "mean": sum(vals)/len(vals), "q75": q(.75), "max": vals[-1]}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def main(root):
    inp, result, val = os.path.join(root, "input"), os.path.join(root, "run", "results"), os.path.join(root, "run", "validation")
    os.makedirs(val, exist_ok=True)
    universe = read_tsv(os.path.join(inp, "universe.tsv"))
    memberships = read_tsv(os.path.join(inp, "memberships.tsv"))
    terms = read_tsv(os.path.join(inp, "terms.tsv"))
    exclusions = read_tsv(os.path.join(inp, "term_exclusions.tsv"))
    identity = read_tsv(os.path.join(inp, "selected_identity.tsv"))
    primary = read_tsv(os.path.join(result, "02_primary_all_terms.tsv"))
    if len(universe) != 41139: raise RuntimeError("Frozen nominal universe mismatch")
    byid = {r["unit_id"]: r for r in universe}
    terms_by = {r["term_id"]: r for r in terms}
    excluded = {r["term_id"] for r in exclusions}
    valid_bp = {t for t, r in terms_by.items() if r["namespace"].strip().lower() in {"biological_process", "biological process", "bp"} and t not in ROOTS and t not in excluded}
    member_pairs = {(r["unit_id"], r["term_id"]) for r in memberships if r["unit_id"] in byid and r["term_id"] in valid_bp}
    U_ids = sorted({u for u, t in member_pairs})
    U_set = set(U_ids)
    if len(U_ids) != 25482: raise RuntimeError(f"Effective U mismatch: {len(U_ids)} != 25482")
    uidx = {u:i for i,u in enumerate(U_ids)}
    opp = {u for u in U_ids if flag(byid[u]["opposite"])}
    A = {u for u in U_ids if flag(byid[u]["opposite"]) and byid[u]["direction"] == "P0up_P1down"}
    B = {u for u in U_ids if flag(byid[u]["opposite"]) and byid[u]["direction"] == "P0down_P1up"}
    strict = {u for u in opp if flag(byid[u]["strict127"])}
    sA, sB = A & strict, B & strict
    counts = {"U":len(U_ids), "overall391_annotated":len(opp), "A_annotated":len(A), "B_annotated":len(B), "strict127_annotated":len(strict), "strictA_annotated":len(sA), "strictB_annotated":len(sB)}
    expected_counts = {"U":25482,"overall391_annotated":213,"A_annotated":184,"B_annotated":29,"strict127_annotated":79,"strictA_annotated":67,"strictB_annotated":12}
    if counts != expected_counts: raise RuntimeError(f"Selection coverage mismatch: {counts}")
    members_by_term = defaultdict(set); terms_by_unit = defaultdict(set)
    for u,t in member_pairs: members_by_term[t].add(u); terms_by_unit[u].add(t)
    tested = sorted(t for t, us in members_by_term.items() if 5 <= len(us) <= math.floor(.95*len(U_ids)))
    if len(tested) != 29: raise RuntimeError(f"Tested term mismatch: {len(tested)} != 29")

    # Independent reconstruction of 58 Fisher tables and exact p values.
    expected = []
    for t in tested:
        us = members_by_term[t]
        a = len(opp & us); b = len(opp) - a; c = len((U_set-opp) & us); d = len(U_set-opp)-c
        expected.append({"term_id":t,"analysis":"A_overall391_vs_annotated_background","selected_n":len(opp),"selected_term_n":a,"background_n":len(U_ids),"background_term_n":len(us),"p_exact_independent":fisher_one_greater(a,b,c,d)})
        a = len(A & us); b = len(A)-a; c = len(B & us); d = len(B)-c
        expected.append({"term_id":t,"analysis":"B_directionA_vs_directionB_within391","selected_n":len(A),"selected_term_n":a,"background_n":len(opp),"background_term_n":a+c,"p_exact_independent":fisher_two_sided(a,b,c,d)})
    for r,q in zip(expected,bh([r["p_exact_independent"] for r in expected])): r["q_exact_independent"] = q
    reported = {(r["term_id"],r["analysis"]):r for r in primary}
    if len(reported) != 58 or len(primary) != 58: raise RuntimeError("Primary output must contain exactly 58 unique rows")
    check_rows=[]; max_p=max_q=0.0
    for r in expected:
        rr=reported.get((r["term_id"],r["analysis"]))
        if rr is None: raise RuntimeError("Missing output row: "+repr((r["term_id"],r["analysis"])))
        for key in ("selected_n","selected_term_n","background_n","background_term_n"):
            r["reported_"+key]=int(float(rr[key])); r["match_"+key]=r[key]==r["reported_"+key]
        r["reported_p_value"] = float(rr["p_value"]); r["reported_q_combined"] = float(rr["q_combined"])
        r["abs_p_difference"] = abs(r["p_exact_independent"]-r["reported_p_value"])
        r["abs_q_difference"] = abs(r["q_exact_independent"]-r["reported_q_combined"])
        max_p=max(max_p,r["abs_p_difference"]); max_q=max(max_q,r["abs_q_difference"]); check_rows.append(r)
    write_tsv(os.path.join(val,"01_fisher_BH_independent_validation.tsv"),check_rows,list(check_rows[0]))

    # Same input-defined 5x5 rank strata and exact conditional hypergeometric QA.
    lengths=[float(byid[u]["length"]) for u in U_ids]; means=[float(byid[u]["baseMean"]) for u in U_ids]
    lq,bq=rank_quintile(lengths),rank_quintile(means)
    strata=defaultdict(list)
    for i,u in enumerate(U_ids): strata[(lq[i],bq[i])].append(u)
    strata_keys=sorted(strata)
    exact_mc=[]
    for t in tested:
        us=members_by_term[t]
        # A right-tail conditional count distribution, preserving each stratum's 391 count.
        pieces=[]
        for s in strata_keys:
            z=strata[s]; pieces.append((len(z),sum(u in us for u in z),sum(u in opp for u in z)))
        off,dist=convolved_counts(pieces); obs=len(opp&us)
        pA=float(dist[max(0,obs-off):].sum()) if obs-off < len(dist) else 0.0
        exact_mc.append({"term_id":t,"analysis":"A_overall391_vs_annotated_background","conditional_exact_p":min(1.0,pA),"observed_statistic":obs,"conditional_expectation":sum(n*k/N for N,k,n in pieces)})
        # B two-sided distribution centered on its analytic conditional expectation.
        piecesB=[]; expx=0.0; totalK=0
        for s in strata_keys:
            z=[u for u in strata[s] if u in opp]; N=len(z); K=sum(u in us for u in z); n=sum(u in A for u in z)
            if N:
                piecesB.append((N,K,n)); expx += n*K/N
            totalK += K
        offB,distB=convolved_counts(piecesB); nA=len(A); nB=len(B); obsA=len(A&us)
        observed = obsA/nA - (totalK-obsA)/nB
        expectation = expx/nA - (totalK-expx)/nB
        xvals=np.arange(offB,offB+len(distB)); dv=xvals/nA-(totalK-xvals)/nB
        pB=float(distB[np.abs(dv-expectation) >= abs(observed-expectation)-1e-12].sum())
        exact_mc.append({"term_id":t,"analysis":"B_directionA_vs_directionB_within391","conditional_exact_p":min(1.0,pB),"observed_statistic":observed,"conditional_expectation":expectation})
    for r,q in zip(exact_mc,bh([r["conditional_exact_p"] for r in exact_mc])): r["conditional_exact_q_BH_2T"] = q
    for r in exact_mc:
        rr=reported[(r["term_id"],r["analysis"])]
        r["reported_MC_p"] = float(rr["p_empirical"]); r["reported_MC_q_2T"] = float(rr["q_perm_combined"])
        r["MC_minus_exact"] = r["reported_MC_p"]-r["conditional_exact_p"]
        r["same_q05_call"] = (r["reported_MC_q_2T"] < .05)==(r["conditional_exact_q_BH_2T"] < .05)
    write_tsv(os.path.join(val,"02_conditional_exact_QA.tsv"),exact_mc,list(exact_mc[0]))

    # Membership and candidate exports.
    ids_by_rep={r.get("representative_gene",""):r for r in identity}
    all391=[]; coverage=[]
    for u in sorted(u for u,r in byid.items() if flag(r["opposite"])):
        valid_terms=sorted(terms_by_unit.get(u,set())); coverage.append({"unit_id":u,"direction":byid[u]["direction"],"strict127":byid[u]["strict127"],"has_valid_BP":bool(valid_terms),"valid_BP_n":len(valid_terms)})
        ident=ids_by_rep.get(u,{}) # no representative key here; identity is separately mapped below
        for t in valid_terms:
            all391.append({"unit_id":u,"direction":byid[u]["direction"],"strict127":byid[u]["strict127"],"baseMean":byid[u]["baseMean"],"length":byid[u]["length"],"term_id":t,"term_name":terms_by[t]["name"],"tested_term":t in tested})
    write_tsv(os.path.join(val,"03_all391_valid_BP_memberships.tsv"),all391,["unit_id","direction","strict127","baseMean","length","term_id","term_name","tested_term"])
    write_tsv(os.path.join(val,"03b_all391_valid_BP_coverage.tsv"),coverage,["unit_id","direction","strict127","has_valid_BP","valid_BP_n"])

    sig=[r["term_id"] for r in expected if r["analysis"].startswith("A_") and r["q_exact_independent"]<.05 and r["selected_term_n"]>=2]
    ident_by_unit={r["unit_id"]:r for r in identity}
    sigmem=[]
    for t in sig:
        for u in sorted(opp & members_by_term[t]):
            ident=ident_by_unit.get(u,{})
            sigmem.append({"term_id":t,"term_name":terms_by[t]["name"],"unit_id":u,"direction":byid[u]["direction"],"strict127":byid[u]["strict127"],"baseMean":byid[u]["baseMean"],"length":byid[u]["length"],"representative_gene":ident.get("representative_gene",""),"annotation_description":ident.get("annotation_description","")})
    write_tsv(os.path.join(val,"04_significant_A_term_members.tsv"),sigmem,["term_id","term_name","unit_id","direction","strict127","baseMean","length","representative_gene","annotation_description"])
    jac=[]
    for i,t1 in enumerate(sig):
        for t2 in sig[i+1:]:
            for scope,base in (("U",U_set),("overall391",opp)):
                a=members_by_term[t1]&base; b=members_by_term[t2]&base
                jac.append({"scope":scope,"term_1":t1,"term_1_name":terms_by[t1]["name"],"term_2":t2,"term_2_name":terms_by[t2]["name"],"intersection_n":len(a&b),"union_n":len(a|b),"jaccard":len(a&b)/len(a|b) if a|b else 0.0})
    write_tsv(os.path.join(val,"05_significant_term_Jaccard.tsv"),jac,["scope","term_1","term_1_name","term_2","term_2_name","intersection_n","union_n","jaccard"])
    cand=[]
    for u,ident in ident_by_unit.items():
        rep=ident.get("representative_gene","")
        if rep in CANDIDATES:
            overlaps=sorted(terms_by_unit.get(u,set()) & set(sig))
            cand.append({"candidate":CANDIDATES[rep],"unit_id":u,"representative_gene":rep,"in_opposite391":flag(byid[u]["opposite"]),"direction":byid[u]["direction"],"strict127":byid[u]["strict127"],"significant_term_membership":bool(overlaps),"term_ids":";".join(overlaps),"term_names":";".join(terms_by[t]["name"] for t in overlaps)})
    found={r["representative_gene"] for r in cand}
    for rep,name in CANDIDATES.items():
        if rep not in found: cand.append({"candidate":name,"unit_id":"","representative_gene":rep,"in_opposite391":False,"direction":"","strict127":"","significant_term_membership":False,"term_ids":"","term_names":""})
    write_tsv(os.path.join(val,"06_six_candidate_significant_term_membership.tsv"),cand,["candidate","unit_id","representative_gene","in_opposite391","direction","strict127","significant_term_membership","term_ids","term_names"])

    dist=[]
    groups={"U":U_set,"overall391":opp,"directionA":A,"directionB":B,"strict127":strict,"strictA":sA,"strictB":sB}
    for group,us in groups.items():
        for measure in ("length","baseMean"):
            row={"group":group,"measure":measure}; row.update(descriptive([float(byid[u][measure]) for u in us])); dist.append(row)
    write_tsv(os.path.join(val,"07_length_baseMean_descriptive.tsv"),dist,["group","measure","n","min","q25","median","mean","q75","max"])
    hashes={"inputs":{f:sha256(os.path.join(inp,f)) for f in ("universe.tsv","memberships.tsv","terms.tsv","term_exclusions.tsv","selected_identity.tsv")},"primary_output":sha256(os.path.join(result,"02_primary_all_terms.tsv"))}
    exact_cross=[r for r in exact_mc if not r["same_q05_call"]]
    summary={"validator":"independent_python_exact_conditional_v1","nominal_counts":counts,"tested_terms":len(tested),"primary_rows":len(primary),"fisher_max_abs_p_difference":max_p,"fisher_max_abs_BH_difference":max_q,"all_contingency_counts_match":all(all(r["match_"+k] for k in ("selected_n","selected_term_n","background_n","background_term_n")) for r in check_rows),"conditional_exact_vs_MC_q05_disagreements":len(exact_cross),"significant_A_terms":sig,"significant_A_term_names":[terms_by[t]["name"] for t in sig],"hashes":hashes,"notes":["Exact conditional QA uses the same fixed 5x5 rank strata and selection sets, but deterministic finite hypergeometric convolution rather than the frozen 10,000-draw Monte Carlo.","No RNA model was refit; no frozen R result or protocol file was modified.","B rows with background_term_n=0 have no informative odds ratio; original output remains untouched and downstream interpretation should render OR as NA for those rows."]}
    with open(os.path.join(val,"00_validation_summary.json"),"w",encoding="utf-8") as h: json.dump(summary,h,ensure_ascii=False,indent=2)
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__ == "__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("root"); args=ap.parse_args(); main(args.root)
