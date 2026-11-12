
"""Re-computation of the atrial-fibrillation (AF) positive-control analyses without R (independent Python re-implementation).

Purpose: pipeline v9 cached an AF lead-SNP list that was built with a p-value string comparison (130 leads); the corrected rule gives 132.
This script (a) reproduces the v9 table with the v9 rule (validation of the re-implementation) and (b) recomputes the AF -> stroke estimates
with the corrected lead list.  Estimators follow TwoSampleMR: multiplicative random-effects IVW (normal-approximation p-value), MR-Egger
(t, n-2 df, standard error scaled by 1/min(1, sigma)) and the weighted median (bootstrap SE, 1,000 draws, fixed seed).
Usage: python af_update_python.py <project_root> <raw_dir> <out_csv>
"""
import sys, os, json, numpy as np, pandas as pd
from scipy.stats import norm, t as tdist, chi2

root, raw, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def scan(path, key, ids, usecols):
    out = []
    for ch in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=2_000_000, dtype=str, compression="gzip"):
        s = ch[ch[key].isin(ids)]
        if len(s): out.append(s)
    return pd.concat(out)

# ---------------- AF lead SNPs
af_cols = ["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "hm_effect_allele_frequency", "standard_error", "p_value"]
parts = []
for ch in pd.read_csv(os.path.join(raw, "GCST006414.h.tsv.gz"), sep="\t", usecols=af_cols, chunksize=2_000_000, dtype=str, compression="gzip"):
    p = pd.to_numeric(ch.p_value, errors="coerce")
    g = ch[(p < 5e-8)]
    if len(g): parts.append(g.assign(p_num=p[g.index]))
allr = pd.concat(parts); allr = allr[allr.hm_rsid.fillna("").str.startswith("rs") & allr.hm_beta.notna() & allr.standard_error.notna()].copy()
allr["chr"] = pd.to_numeric(allr.hm_chrom, errors="coerce"); allr = allr[allr.chr.notna()].copy(); allr["chr"] = allr.chr.astype(int)
allr["pos"] = allr.hm_pos.astype(float); allr["beta"] = allr.hm_beta.astype(float); allr["se"] = allr.standard_error.astype(float)
allr["z2"] = (allr.beta / allr.se) ** 2

def prune(a):
    a = a.sort_values("z2", ascending=False, kind="mergesort").reset_index(drop=True)
    removed = np.zeros(len(a), bool); lead = []
    for i in range(len(a)):
        if removed[i]: continue
        lead.append(i); removed |= (a.chr.values == a.chr.values[i]) & (np.abs(a.pos.values - a.pos.values[i]) <= 500000)
    return a.iloc[lead].drop_duplicates("hm_rsid")

sets = {"corrected": prune(allr[allr.p_num < 5e-8])}
sets["v9_cached_rule"] = prune(allr[allr.p_value < "5e-08"])   # numeric p < 5e-8 AND (lexicographic) p string < "5e-08": the v9 cache bug
print({k: len(v) for k, v in sets.items()}, flush=True)

# ---------------- outcomes
ids = set().union(*[set(v.hm_rsid) for v in sets.values()])
OUT_COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"]
out = {}
for k, f in (("CES", "GCST90104541"), ("LAS", "GCST90104542"), ("SVS", "GCST90104543")):
    o = scan(os.path.join(raw, f + ".h.tsv.gz"), "rsid", ids, OUT_COLS).drop_duplicates("rsid")
    out[k] = pd.DataFrame(dict(SNP=o.rsid, ea_o=o.effect_allele, oa_o=o.other_allele, beta_o=o.beta, se_o=o.standard_error, eaf_o=o.effect_allele_frequency))
o = scan(os.path.join(raw, "GCST006910.h.tsv.gz"), "hm_rsid", ids, ["hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "hm_effect_allele_frequency"]).drop_duplicates("hm_rsid")
out["MEGA_CES"] = pd.DataFrame(dict(SNP=o.hm_rsid, ea_o=o.hm_effect_allele, oa_o=o.hm_other_allele, beta_o=o.hm_beta, se_o=o.standard_error, eaf_o=o.hm_effect_allele_frequency))

def harmonise(exp, out, tol=0.08):
    m = exp.merge(out, on="SNP", how="inner"); rows = []
    for r in m.itertuples(index=False):
        a1e, a2e, a1o, a2o = r.ea_e, r.oa_e, str(r.ea_o).upper(), str(r.oa_o).upper()
        bo = float(r.beta_o); eo = float(r.eaf_o) if r.eaf_o not in (None, "NA", "") and pd.notna(r.eaf_o) else np.nan
        ee = r.eaf_e; pal = (a1e, a2e) in (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")); keep, flip = True, False
        if not pal:
            if (a1e, a2e) == (a1o, a2o) or (COMP.get(a1e), COMP.get(a2e)) == (a1o, a2o): pass
            elif (a1e, a2e) == (a2o, a1o) or (COMP.get(a1e), COMP.get(a2e)) == (a2o, a1o): flip = True
            else: keep = False
        elif np.isnan(eo) or abs(ee - 0.5) < tol or abs(eo - 0.5) < tol: keep = False
        else: flip = (ee - 0.5) * (eo - 0.5) < 0
        rows.append(dict(SNP=r.SNP, keep=keep, bx=r.beta_e, sex=r.se_e, by=-bo if flip else bo, sey=float(r.se_o)))
    h = pd.DataFrame(rows); return h[h.keep].reset_index(drop=True)

def ivw(bx, by, sey):
    w = 1 / sey**2; n = len(bx); b = np.sum(w * bx * by) / np.sum(w * bx**2); res = by - b * bx
    q = np.sum(w * res**2); sigma = np.sqrt(q / (n - 1)); se = np.sqrt(1 / np.sum(w * bx**2)) * sigma / min(1, sigma)
    return b, se, 2 * norm.sf(abs(b / se)), chi2.sf(q, n - 1)

def egger(bx, by, sey):
    sgn = np.where(bx < 0, -1, 1); bx2, by2 = bx * sgn, by * sgn; w = 1 / sey**2; n = len(bx)
    X = np.column_stack([np.ones(n), bx2]); W = np.diag(w) if n < 400 else None
    XtWX = X.T @ (X * w[:, None]); beta = np.linalg.solve(XtWX, X.T @ (w * by2)); res = by2 - X @ beta
    sigma = np.sqrt(np.sum(w * res**2) / (n - 2)); cov = np.linalg.inv(XtWX) * sigma**2
    se = np.sqrt(np.diag(cov)) / min(1, sigma)
    p = 2 * tdist.sf(np.abs(beta / se), n - 2)
    return beta[1], se[1], p[1], beta[0], se[0], p[0]

def wmed(b, w):
    o = np.argsort(b); bo, wo = b[o], w[o]; cs = np.cumsum(wo) - 0.5 * wo; cs = cs / wo.sum()
    below = np.max(np.where(cs < 0.5)[0]); return bo[below] + (bo[below + 1] - bo[below]) * (0.5 - cs[below]) / (cs[below + 1] - cs[below])

def weighted_median(bx, by, sex, sey, nboot=1000, seed=20261105):
    biv = by / bx; vb = sey**2 / bx**2 + by**2 * sex**2 / bx**4; w = 1 / vb; b = wmed(biv, w)
    rng = np.random.default_rng(seed); meds = np.empty(nboot)
    for i in range(nboot):
        bxb = rng.normal(bx, sex); byb = rng.normal(by, sey); meds[i] = wmed(byb / bxb, w)
    se = meds.std(ddof=1); return b, se, 2 * norm.sf(abs(b / se))

rows = []
for sname, lead in sets.items():
    afx = pd.DataFrame(dict(SNP=lead.hm_rsid, ea_e=lead.hm_effect_allele, oa_e=lead.hm_other_allele, eaf_e=lead.hm_effect_allele_frequency.astype(float), beta_e=lead.beta, se_e=lead.se))
    for k in ("CES", "LAS", "SVS", "MEGA_CES"):
        hh = harmonise(afx, out[k])
        bx, by, sex, sey = hh.bx.values, hh.by.values, hh.sex.values, hh.sey.values
        b, se, p, qp = ivw(bx, by, sey); eb, ese, ep, ic, ice, icp = egger(bx, by, sey); wb, wse, wp = weighted_median(bx, by, sex, sey)
        for meth, bb, ss, pp in (("Inverse variance weighted", b, se, p), ("MR Egger", eb, ese, ep), ("Weighted median", wb, wse, wp)):
            rows.append(dict(lead_set=sname, n_lead=len(lead), outcome=k, method=meth, n_snps=len(hh), beta=bb, se=ss, OR=np.exp(bb), lci=np.exp(bb - 1.96 * ss), uci=np.exp(bb + 1.96 * ss), p=pp, Q_P=qp, Egger_intercept_P=icp))
pd.DataFrame(rows).to_csv(out_csv, index=False)
sets["corrected"][["hm_rsid", "chr", "pos"]].to_csv(out_csv.replace(".csv", "_leads_corrected.csv"), index=False)
sets["v9_cached_rule"][["hm_rsid", "chr", "pos"]].to_csv(out_csv.replace(".csv", "_leads_v9rule.csv"), index=False)
print("done", flush=True)
