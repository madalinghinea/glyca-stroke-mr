"""Independent re-implementation of the core MR estimators (no R, no TwoSampleMR).

Reads the supplementary tables written by the pipeline (tables/*.csv) and checks that the reported
IVW (random effects), MR-Egger, Cochran Q, Egger intercept and MR-PRESSO observed RSS are reproduced
from the SNP-level harmonised data (Tables S2, S3, S9e).  Usage:  python independent_check.py [project_root]
"""
import sys, os
import numpy as np, pandas as pd

root = sys.argv[1] if len(sys.argv) > 1 else "."
T = lambda f: pd.read_csv(os.path.join(root, "tables", f))

def ivw_random(bx, by, sey):
    w = 1 / sey**2
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    res = by - b * bx
    q = np.sum(w * res**2); df = len(bx) - 1
    se = np.sqrt(1 / np.sum(w * bx**2)) * max(np.sqrt(q / df), 1.0)   # multiplicative random effects
    return b, se, q, df

def egger(bx, by, sey):
    s = np.where(bx < 0, -1.0, 1.0); x, y, w = np.abs(bx), by * s, 1 / sey**2
    X = np.column_stack([np.ones(len(x)), x]); A = X.T @ (w[:, None] * X)
    beta = np.linalg.solve(A, X.T @ (w * y)); res = y - X @ beta
    q = np.sum(w * res**2); df = len(x) - 2
    se = np.sqrt(np.diag(np.linalg.inv(A))) * max(np.sqrt(q / df), 1.0)
    return beta, se, q, df

def presso_rss(bx, by, sey):
    w = 1 / sey**2; rss = 0.0
    for i in range(len(bx)):
        m = np.arange(len(bx)) != i
        b = np.sum(w[m] * bx[m] * by[m]) / np.sum(w[m] * bx[m]**2)
        rss += w[i] * (by[i] - b * bx[i])**2
    return rss

sets = {"MEGASTROKE": T("TableS2_MEGASTROKE_harmonised.csv"), "GIGASTROKE": T("TableS3_GIGASTROKE_harmonised.csv")}
s9e = T("TableS9e_subtype_harmonised_with_steiger.csv")
for k in ("CES", "LAS", "SVS"): sets[k] = s9e[s9e["subtype"] == k]

s4, s5a, s5b, s5c = T("TableS4_MR_results.csv"), T("TableS5a_heterogeneity.csv"), T("TableS5b_pleiotropy_egger_intercept.csv"), T("TableS5c_presso_summary.csv")
s9a, s9b, s9c = T("TableS9a_subtype_MR_results.csv"), T("TableS9b_subtype_heterogeneity.csv"), T("TableS9c_subtype_pleiotropy.csv")

def rep(tab, key, col, **kw):                     # reported value from a table
    t = tab
    for c, v in kw.items(): t = t[t[c] == v]
    return float(t[col].iloc[0])

rows = []
for k, d in sets.items():
    if "mr_keep" in d: d = d[d["mr_keep"].astype(str).str.upper() == "TRUE"]
    bx, by, sey = d["beta.exposure"].values, d["beta.outcome"].values, d["se.outcome"].values
    b, se, q, df = ivw_random(bx, by, sey); be, see, qe, dfe = egger(bx, by, sey); rss = presso_rss(bx, by, sey)
    ov = k in ("MEGASTROKE", "GIGASTROKE")
    t4, th, tp = (s4, s5a, s5b) if ov else (s9a, s9b, s9c)
    kc = "outcome_label" if ov else "subtype"
    ivw_row = dict(method="Inverse variance weighted"); egg_row = dict(method="MR Egger")
    rows += [
        (k, "IVW beta", b, rep(t4, k, "b", **{kc: k}, **ivw_row)), (k, "IVW SE", se, rep(t4, k, "se", **{kc: k}, **ivw_row)),
        (k, "Egger slope", be[1], rep(t4, k, "b", **{kc: k}, **egg_row)), (k, "Egger SE (slope)", see[1], rep(t4, k, "se", **{kc: k}, **egg_row)),
        (k, "Cochran Q (IVW)", q, rep(th, k, "Q", **{kc: k}, **ivw_row)), (k, "Cochran Q (Egger)", qe, rep(th, k, "Q", **{kc: k}, **egg_row)),
        (k, "Egger intercept", be[0], rep(tp, k, "egger_intercept", **{kc: k})),
        (k, "MR-PRESSO observed RSS", rss, rep(s5c, k, "global_rss", outcome_label=k) if k in s5c["outcome_label"].values else np.nan)]
out = pd.DataFrame(rows, columns=["outcome", "quantity", "python", "reported"])
out["abs_diff"] = (out["python"] - out["reported"]).abs()
pd.set_option("display.width", 200); pd.set_option("display.float_format", lambda x: f"{x:.5f}")
print(out.to_string(index=False))
tol = {"IVW beta": 5e-4, "IVW SE": 5e-4, "Egger slope": 5e-4, "Egger SE (slope)": 5e-4, "Cochran Q (IVW)": 0.006, "Cochran Q (Egger)": 0.006,
       "Egger intercept": 5e-5, "MR-PRESSO observed RSS": 0.006}       # tolerance = rounding of the published tables
bad = out[(out["abs_diff"] > out["quantity"].map(tol)) & out["reported"].notna()]
print("\nMismatches beyond table rounding:", len(bad)); print(bad.to_string(index=False) if len(bad) else "none")
out.to_csv(os.path.join(root, "tables", "independent_check_results.csv"), index=False)
sys.exit(1 if len(bad) else 0)
