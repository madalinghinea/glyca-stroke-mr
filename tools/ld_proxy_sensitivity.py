"""LD-proxy sensitivity analysis for the GlycA instruments that are absent from the GIGASTROKE summary statistics.

GIGASTROKE statistics are matched by rsID, so 14 of the 60 GlycA instruments are missing. For each missing instrument this script
 1. converts its GRCh37 position to GRCh38 (UCSC liftOver chain file) and collects GIGASTROKE variants with an rsID within +/- 50 kb,
 2. obtains European LD (1000 Genomes) from the OpenGWAS API (POST /ld/matrix; needs the personal token in the environment variable
    OPENGWAS_JWT, which is never written to disk by this script),
 3. takes the proxy with the highest r2 >= 0.8, aligns the proxy effect to the instrument effect allele by LD phase and
 4. re-estimates IVW, MR-Egger and weighted-median effects with and without the proxies (Table S15, S15b).
Estimators follow TwoSampleMR (random-effects IVW with normal p-value; Egger with t distribution; weighted median with a seeded bootstrap).

Usage: python ld_proxy_sensitivity.py <tables_dir> <raw_dir> <chain_file>      (about 3 minutes; reads four GIGASTROKE files)
"""
import os, sys, time, json
import numpy as np, pandas as pd, requests
from pyliftover import LiftOver
from scipy.stats import norm, t as tdist, chi2

tables, raw, chain = sys.argv[1], sys.argv[2], sys.argv[3]
tok = os.environ.get("OPENGWAS_JWT") or sys.exit("set OPENGWAS_JWT")
HDR = {"Authorization": "Bearer " + tok}
FILES = {"GIGASTROKE": "GCST90104540", "CES": "GCST90104541", "LAS": "GCST90104542", "SVS": "GCST90104543"}
COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"]
MANUAL_GRCH38 = {"rs687621": (9, 133261662)}        # liftOver fails for this variant; position from Ensembl (GRCh38)
WINDOW, R2_MIN, SEED = 50_000, 0.8, 20261105

s1 = pd.read_csv(os.path.join(tables, "TableS1_GlycA_instruments.csv")).set_index("SNP", drop=False)
s3 = pd.read_csv(os.path.join(tables, "TableS3_GIGASTROKE_harmonised.csv")); s3 = s3[s3.mr_keep.astype(str).str.upper() == "TRUE"]
s9e = pd.read_csv(os.path.join(tables, "TableS9e_subtype_harmonised_with_steiger.csv"))
base = {"GIGASTROKE": s3, **{k: s9e[s9e.subtype == k] for k in ("CES", "LAS", "SVS")}}
missing = [s for s in s1.SNP if s not in set(pd.read_csv(os.path.join(tables, "TableS3_GIGASTROKE_harmonised.csv")).SNP)]
print(len(missing), "instruments absent from GIGASTROKE", flush=True)

lo = LiftOver(chain)
pos38 = {}
for s in missing:
    if s in MANUAL_GRCH38: pos38[s] = MANUAL_GRCH38[s]; continue
    c = lo.convert_coordinate(f"chr{int(s1.loc[s, 'chr'])}", int(s1.loc[s, "pos_hg19"]) - 1)
    if c: pos38[s] = (int(s1.loc[s, "chr"]), c[0][1] + 1)

cand = []
for ch in pd.read_csv(os.path.join(raw, FILES["GIGASTROKE"] + ".h.tsv.gz"), sep="\t", usecols=COLS, chunksize=2_000_000, dtype=str, compression="gzip"):
    ch = ch[ch.rsid.fillna("NA").str.startswith("rs")]; c_ = pd.to_numeric(ch.chromosome, errors="coerce"); p_ = pd.to_numeric(ch.base_pair_location, errors="coerce")
    for s, (chrom, pos) in pos38.items():
        sel = ch[(c_ == chrom) & (p_.between(pos - WINDOW, pos + WINDOW))]
        if len(sel): cand.append(sel.assign(target=s, dist=(p_[sel.index] - pos).abs()))
cand = pd.concat(cand)

proxies, unrecovered = {}, []
for s in missing:
    g = cand[(cand.target == s) & (cand.rsid != s)].drop_duplicates("rsid").sort_values("dist").head(450)
    r = requests.post("https://api.opengwas.io/api/ld/matrix", headers=HDR, json={"rsid": [s] + g.rsid.tolist(), "pop": "EUR"}, timeout=180)
    names = r.json().get("snplist", []) if r.status_code == 200 else []
    lab = {n.split("_")[0]: n for n in names}
    if s not in lab: unrecovered.append(dict(Instrument=s, Reason="not in the 1000 Genomes EUR LD reference (indel)")); continue
    M = np.array(r.json()["matrix"], dtype=float); ti = names.index(lab[s])
    best = max(((M[ti, k], n) for k, n in enumerate(names) if k != ti), key=lambda x: x[0] ** 2)
    if best[0] ** 2 < R2_MIN: unrecovered.append(dict(Instrument=s, Reason=f"best proxy {best[1].split('_')[0]} r2 = {best[0]**2:.2f} (< {R2_MIN})")); continue
    proxies[s] = (lab[s], best[1], best[0]); time.sleep(0.5)

pdata = {}
for k, f in FILES.items():
    parts = []
    for ch in pd.read_csv(os.path.join(raw, f + ".h.tsv.gz"), sep="\t", usecols=COLS, chunksize=2_000_000, dtype=str, compression="gzip"):
        sel = ch[ch.rsid.isin({v[1].split("_")[0] for v in proxies.values()})]
        if len(sel): parts.append(sel)
    pdata[k] = pd.concat(parts).drop_duplicates("rsid").set_index("rsid")

used = []
for s, (ilab, plab, r) in proxies.items():
    X1, X2 = ilab.split("_")[1:3]; Y1, Y2 = plab.split("_")[1:3]; ea = s1.loc[s, "effect_allele"]
    partner = (Y1 if r > 0 else Y2) if ea == X1 else (Y2 if r > 0 else Y1)
    rs = plab.split("_")[0]
    for k, d in pdata.items():
        if rs not in d.index: continue
        x = d.loc[rs]; E, O = x.effect_allele.upper(), x.other_allele.upper()
        if {E, O} != {Y1, Y2}: continue
        used.append(dict(outcome=k, Instrument=s, Proxy_SNP=rs, r=r, r2=r * r, Beta_outcome_aligned=(1 if partner == E else -1) * float(x.beta), SE_outcome=float(x.standard_error), p_out=float(x.p_value)))
used = pd.DataFrame(used)

def ivw(bx, by, sey):
    w = 1 / sey**2; n = len(bx); b = np.sum(w * bx * by) / np.sum(w * bx**2); q = np.sum(w * (by - b * bx)**2)
    se = np.sqrt(1 / np.sum(w * bx**2)) * max(np.sqrt(q / (n - 1)), 1.0); return b, se, 2 * norm.sf(abs(b / se)), chi2.sf(q, n - 1)
def egger(bx, by, sey):
    sg = np.where(bx < 0, -1, 1); x, y = bx * sg, by * sg; w = 1 / sey**2; n = len(bx); X = np.column_stack([np.ones(n), x])
    A = X.T @ (X * w[:, None]); b = np.linalg.solve(A, X.T @ (w * y)); sig = np.sqrt(np.sum(w * (y - X @ b)**2) / (n - 2))
    se = np.sqrt(np.diag(np.linalg.inv(A))) * sig / min(1, sig); p = 2 * tdist.sf(np.abs(b / se), n - 2); return b[1], p[1], p[0]
def wmed(b, w):
    o = np.argsort(b); bo, wo = b[o], w[o]; cs = (np.cumsum(wo) - 0.5 * wo) / wo.sum(); i = np.max(np.where(cs < 0.5)[0])
    return bo[i] + (bo[i + 1] - bo[i]) * (0.5 - cs[i]) / (cs[i + 1] - cs[i])
def wm(bx, by, sex, sey, nboot=1000):
    w = 1 / (sey**2 / bx**2 + by**2 * sex**2 / bx**4); b = wmed(by / bx, w); rng = np.random.default_rng(SEED)
    se = np.std([wmed(rng.normal(by, sey) / rng.normal(bx, sex), w) for _ in range(nboot)], ddof=1); return b, 2 * norm.sf(abs(b / se))

rows = []
for k, b in base.items():
    a = used[used.outcome == k]
    ext = pd.concat([b[["SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"]],
                     pd.DataFrame({"SNP": a.Instrument.values, "beta.exposure": s1.loc[a.Instrument, "beta_SD"].values, "se.exposure": s1.loc[a.Instrument, "se"].values,
                                   "beta.outcome": a.Beta_outcome_aligned.values, "se.outcome": a.SE_outcome.values})], ignore_index=True)
    for label, d in (("rsID matching only (primary analysis)", b), ("plus LD proxies (r2 >= 0.8)", ext)):
        bx, by, sex, sey = (d[c].values for c in ("beta.exposure", "beta.outcome", "se.exposure", "se.outcome"))
        bi, si, pi, qp = ivw(bx, by, sey); be, pe, pint = egger(bx, by, sey); bw, pw = wm(bx, by, sex, sey)
        rows.append({"Outcome": "GIGASTROKE (overall)" if k == "GIGASTROKE" else k, "Instrument_set": label, "N_SNPs": len(d), "IVW_OR": np.exp(bi), "IVW_CI_lower": np.exp(bi - 1.96 * si),
                     "IVW_CI_upper": np.exp(bi + 1.96 * si), "IVW_P": pi, "Cochran_Q_P": qp, "Egger_OR": np.exp(be), "Egger_P": pe, "Egger_intercept_P": pint, "Weighted_median_OR": np.exp(bw), "Weighted_median_P": pw})
t15 = pd.DataFrame(rows)
for c in t15.columns[3:]: t15[c] = t15[c].astype(float).map(lambda x: float(f"{x:.4g}"))
t15.to_csv(os.path.join(tables, "TableS15_LD_proxy_sensitivity.csv"), index=False)
out = pd.concat([used.assign(Status="proxy used"), pd.DataFrame(unrecovered).assign(Status="not recovered")], ignore_index=True)
out.to_csv(os.path.join(tables, "TableS15b_LD_proxies_used.csv"), index=False); print(t15[["Outcome", "Instrument_set", "N_SNPs", "IVW_OR", "IVW_P"]].to_string())
