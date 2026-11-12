"""Shared functions for the OpenGWAS-API analyses (reverse MR, GlycA + AF MVMR): API access with a local cache, allele alignment,
univariable estimators (TwoSampleMR logic) and multivariable estimators (TwoSampleMR mv_multiple; MVMR strength_mvmr / pleiotropy_mvmr with gencov = 0).
The OpenGWAS token is read from the environment variable OPENGWAS_JWT only when a request is sent (never written to disk)."""
import os, sys, json, time, re, hashlib
import numpy as np, pandas as pd, requests
from scipy.stats import norm, t as tdist, chi2

API = "https://api.opengwas.io/api/"
CFG = dict(cache=None, raw=None)
def init(cache_dir, raw_dir):
    CFG["cache"], CFG["raw"] = cache_dir, raw_dir
    os.makedirs(cache_dir, exist_ok=True)
def _hdr():
    tok = os.environ.get("OPENGWAS_JWT") or sys.exit("set OPENGWAS_JWT (or run with a populated cache)")
    return {"Authorization": "Bearer " + tok, "User-Agent": "glyca-mr-check"}
SEED, COMP = 20261105, {"A": "T", "T": "A", "C": "G", "G": "C"}
COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"]
N_GLYCA = 115078

def post(path, body, tries=4):
    for k in range(tries):
        r = requests.post(API + path, headers=_hdr(), json=body, timeout=300)
        if r.status_code == 200: return r.json()
        if r.status_code in (429, 502, 503, 504): time.sleep(5 * (k + 1)); continue
        raise RuntimeError(f"{path}: HTTP {r.status_code} {r.text[:200]}")
    raise RuntimeError(path + ": retries exhausted")

def cached(name, fn):
    p = os.path.join(CFG["cache"], name + ".json")
    if os.path.exists(p): return json.load(open(p))
    v = fn(); json.dump(v, open(p, "w")); return v

def tophits(gid):
    return pd.DataFrame(cached("tophits_" + gid, lambda: post("tophits", {"id": [gid], "pval": 5e-8, "preclumped": 1, "clump": 1, "r2": 0.001, "kb": 10000, "pop": "EUR"})))

def _assoc_raw(gid, snps, proxies, chunk=64):
    def run():
        out = []
        for i in range(0, len(snps), chunk):
            out += post("associations", {"id": [gid], "variant": list(snps[i:i + chunk]), "proxies": int(proxies), "r2": 0.8, "align_alleles": 1, "palindromes": 1, "maf_threshold": 0.3})
        return out
    key = "assoc_" + gid + ("_px_" if proxies else "_") + hashlib.md5(",".join(sorted(snps)).encode()).hexdigest()[:12]
    return pd.DataFrame(cached(key, run))

def assoc(gid, snps, chunk=64, proxies=True):
    """POST /associations (N(id) * N(variant) <= 64 per request). Exact rsID matches first; variants absent from the dataset are then looked up with LD proxies
    (r2 >= 0.8, EUR; proxy alleles aligned to the target variant by the API), as TwoSampleMR::extract_outcome_data does. Result indexed by the target rsID."""
    snps = list(snps)
    d = _assoc_raw(gid, snps, False, chunk)
    if len(d): d = d.drop_duplicates("rsid")
    if proxies:
        missing = [s for s in snps if len(d) == 0 or s not in set(d.rsid)]
        if missing:
            px = _assoc_raw(gid, missing, True, chunk)
            if len(px):
                if "target_snp" in px.columns: px["rsid"] = px["target_snp"].where(px["target_snp"].notna(), px["rsid"])
                px = px[px.rsid.isin(missing)].drop_duplicates("rsid"); d = pd.concat([d, px], ignore_index=True)
    if len(d): d = d.drop_duplicates("rsid").set_index("rsid", drop=False)
    return d

def clump(rsids, pvals, name):
    pv = [max(float(p), 1e-300) for p in pvals]
    return set(cached("clump_" + name, lambda: post("ld/clump", {"rsid": list(rsids), "pval": pv, "pthresh": 1, "r2": 0.001, "kb": 10000, "pop": "EUR"})))

def pal(a, b): return (a, b) in (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"))
def align(ref, x, tol=0.08):
    """Align dataset x to the reference alleles/frequency (TwoSampleMR harmonise_data, action = 2). ref, x: dicts with ea, nea, eaf, beta.
    Returns the sign (+1 / -1) to apply to x['beta'], or None if the SNP must be dropped."""
    a1, a2, b1, b2 = ref["ea"], ref["nea"], x["ea"], x["nea"]
    if None in (a1, a2, b1, b2): return None
    if not pal(a1, a2):
        if (a1, a2) == (b1, b2) or (COMP[a1], COMP[a2]) == (b1, b2): return 1
        if (a1, a2) == (b2, b1) or (COMP[a1], COMP[a2]) == (b2, b1): return -1
        return None
    e1, e2 = ref["eaf"], x["eaf"]
    if e1 is None or e2 is None or pd.isna(e1) or pd.isna(e2) or abs(e1 - 0.5) < tol or abs(e2 - 0.5) < tol: return None
    return 1 if (e1 - 0.5) * (e2 - 0.5) > 0 else -1

# ---------------------------------------------------------------- estimators
def ivw(bx, by, sey):
    w = 1 / sey**2; n = len(bx); b = np.sum(w * bx * by) / np.sum(w * bx**2); q = np.sum(w * (by - b * bx)**2)
    se = np.sqrt(1 / np.sum(w * bx**2)) * max(np.sqrt(q / (n - 1)), 1.0); return b, se, 2 * norm.sf(abs(b / se)), chi2.sf(q, n - 1)
def egger(bx, by, sey):
    sg = np.where(bx < 0, -1, 1); x, y = bx * sg, by * sg; w = 1 / sey**2; n = len(bx); X = np.column_stack([np.ones(n), x])
    A = X.T @ (X * w[:, None]); b = np.linalg.solve(A, X.T @ (w * y)); sig = np.sqrt(np.sum(w * (y - X @ b)**2) / (n - 2))
    se = np.sqrt(np.diag(np.linalg.inv(A))) * sig / min(1, sig); p = 2 * tdist.sf(np.abs(b / se), n - 2); return b[1], se[1], p[1], p[0]
def wmed(b, w):
    o = np.argsort(b); bo, wo = b[o], w[o]; cs = (np.cumsum(wo) - 0.5 * wo) / wo.sum(); i = np.max(np.where(cs < 0.5)[0])
    return bo[i] + (bo[i + 1] - bo[i]) * (0.5 - cs[i]) / (cs[i + 1] - cs[i])
def wm(bx, by, sex, sey, nboot=1000):
    w = 1 / (sey**2 / bx**2 + by**2 * sex**2 / bx**4); b = wmed(by / bx, w); rng = np.random.default_rng(SEED)
    se = np.std([wmed(rng.normal(by, sey) / rng.normal(bx, sex), w) for _ in range(nboot)], ddof=1); return b, se, 2 * norm.sf(abs(b / se))

def mv_ivw(BX, by, sey):
    """TwoSampleMR 0.7.13 mv_multiple: weighted regression through the origin; SE = lm SE (fixed-effect SE x residual sigma, also when sigma < 1); normal p.
    (An earlier version divided by min(1, sigma), the rule of the univariable mr_ivw; that is wrong for mv_multiple when sigma < 1 and was
    detected by comparison with the R run of 11 Nov 2026, GlycA + AF to cardioembolic stroke, sigma = 0.931.)"""
    w = 1 / sey**2; K, L = BX.shape; A = BX.T @ (BX * w[:, None]); b = np.linalg.solve(A, BX.T @ (w * by))
    res = by - BX @ b; sig = np.sqrt(np.sum(w * res**2) / (K - L)); se = np.sqrt(np.diag(np.linalg.inv(A))) * sig
    return b, se, 2 * norm.sf(np.abs(b / se)), sig
def cond_F(BX, SEX, j):
    """MVMR::strength_mvmr with gencov = 0: delta from unweighted OLS of exposure j on the other exposures (no intercept);
    sigma2 = se_j^2 + sum_l delta_l^2 se_l^2; F = Q / (K - L + 1)."""
    K, L = BX.shape; oth = [l for l in range(L) if l != j]
    d = np.linalg.lstsq(BX[:, oth], BX[:, j], rcond=None)[0]
    s2 = SEX[:, j]**2 + (SEX[:, oth]**2) @ (d**2)
    return float(np.sum((BX[:, j] - BX[:, oth] @ d)**2 / s2) / (K - L + 1))
def q_a(BX, SEX, by, sey):
    """MVMR::pleiotropy_mvmr (estimator = 'ivw', gencov = 0): Q evaluated at the IVW coefficients; sigma2_A = se_y^2 + sum_j b_j^2 se_j^2."""
    K, L = BX.shape; b = mv_ivw(BX, by, sey)[0]
    q = float(np.sum((by - BX @ b)**2 / (sey**2 + (SEX**2) @ (b**2)))); return q, float(chi2.sf(q, K - L))

def mvmr_dataset(exp_ids, labels, out_fn, name):
    """Union of the exposures' genome-wide significant SNPs, joint clumping, effects on every exposure and the outcome."""
    th = {g: tophits(g) for g in exp_ids}
    allsnp = pd.concat([d[["rsid", "p"]].assign(src=g) for g, d in th.items()]).groupby("rsid").p.min()
    keep = clump(allsnp.index, allsnp.values, name + "_union")
    snps = sorted(keep)
    ex = {g: assoc(g, snps) for g in exp_ids}
    out = out_fn(snps)
    rows = []
    for s in snps:
        if any(s not in ex[g].index for g in exp_ids) or s not in out.index: continue
        r0 = ex[exp_ids[0]].loc[s]; ref = dict(ea=r0.ea, nea=r0.nea, eaf=r0.eaf); bx, sx, ok = [float(r0.beta)], [float(r0.se)], True
        for g in exp_ids[1:]:
            r = ex[g].loc[s]; sg = align(ref, dict(ea=r.ea, nea=r.nea, eaf=r.eaf))
            if sg is None: ok = False; break
            bx.append(sg * float(r.beta)); sx.append(float(r.se))
        if not ok: continue
        ro = out.loc[s]; sg = align(ref, dict(ea=ro.ea, nea=ro.nea, eaf=ro.eaf))
        if sg is None: continue
        rows.append(dict(SNP=s, EA=ref["ea"], NEA=ref["nea"], **{f"beta.{l}": b for l, b in zip(labels, bx)}, **{f"se.{l}": v for l, v in zip(labels, sx)},
                         beta_out=sg * float(ro.beta), se_out=float(ro.se)))
    return pd.DataFrame(rows), len(snps)

def mvmr_summary(d, labels, outcome, analysis):
    BX = d[[f"beta.{l}" for l in labels]].values; SEX = d[[f"se.{l}" for l in labels]].values; by, sey = d.beta_out.values, d.se_out.values
    b, se, p, sig = mv_ivw(BX, by, sey); qa, qap = q_a(BX, SEX, by, sey); rows = []
    for j, l in enumerate(labels):
        rows.append(dict(Analysis=analysis, Outcome=outcome, Exposure=l, N_SNPs=len(d), Beta=b[j], SE=se[j], P=p[j], OR=np.exp(b[j]), OR_lci=np.exp(b[j] - 1.96 * se[j]),
                         OR_uci=np.exp(b[j] + 1.96 * se[j]), Conditional_F=cond_F(BX, SEX, j), Q_A=qa, Q_A_df=len(d) - len(labels), Q_A_P=qap, Residual_sigma=sig))
    return rows

def read_local(fname, snps):
    """Effects of the given rsIDs in a locally downloaded GWAS Catalog file (harmonised format; the older files use hm_* column names)."""
    path = os.path.join(CFG["raw"], fname + ".h.tsv.gz")
    hdr = list(pd.read_csv(path, sep="\t", nrows=0, compression="gzip").columns)
    m = (dict(rsid="rsid", ea="effect_allele", nea="other_allele", eaf="effect_allele_frequency", beta="beta", se="standard_error") if "rsid" in hdr else
         dict(rsid="hm_rsid", ea="hm_effect_allele", nea="hm_other_allele", eaf="hm_effect_allele_frequency", beta="hm_beta", se="standard_error"))
    ids, parts = set(snps), []
    for ch in pd.read_csv(path, sep="\t", usecols=list(m.values()), chunksize=2_000_000, dtype=str, compression="gzip"):
        s = ch[ch[m["rsid"]].isin(ids)]
        if len(s): parts.append(s)
    o = pd.concat(parts).drop_duplicates(m["rsid"])
    o = pd.DataFrame(dict(rsid=o[m["rsid"]], ea=o[m["ea"]].str.upper(), nea=o[m["nea"]].str.upper(), eaf=pd.to_numeric(o[m["eaf"]], errors="coerce"),
                          beta=pd.to_numeric(o[m["beta"]], errors="coerce"), se=pd.to_numeric(o[m["se"]], errors="coerce"))).dropna(subset=["beta", "se"])
    return o.set_index("rsid", drop=False)
