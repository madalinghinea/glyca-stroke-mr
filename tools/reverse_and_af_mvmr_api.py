"""Reverse-direction MR (ischemic stroke -> GlycA) and GlycA + atrial fibrillation multivariable MR for cardioembolic stroke.

Python re-implementation of pipeline step 11 (11_reverse_and_af_mvmr_token.qmd) that talks to the OpenGWAS REST API directly
(POST /tophits, /associations, /ld/clump) through mr_api_lib.py. It needs the personal OpenGWAS token in the environment variable
OPENGWAS_JWT (never written to disk); responses are cached in <cache_dir> (public summary statistics only).
Estimators follow TwoSampleMR (random-effects IVW, normal p; MR-Egger t; seeded weighted median; mv_multiple) and MVMR (strength_mvmr and
pleiotropy_mvmr, gencov = 0). The multivariable code is validated against the saved R dataset of the published GlycA + BMI + smoking analysis
(Tables S6/S7/S7b), which it must reproduce.

Usage: python reverse_and_af_mvmr_api.py <tables_dir> <raw_dir> <cache_dir> <mvdat_446.csv>      (about 20-30 minutes; the API proxy search is slow)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import mr_api_lib as L
from mr_api_lib import *

tables, raw, cache, mvcsv = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
L.init(cache, raw)

# ================================================================ 0a. validation of the multivariable estimators on the saved R dataset (446 SNPs)
print("== validation 1: saved R dataset (GlycA + BMI + smoking -> MEGASTROKE)", flush=True)
mv = pd.read_csv(mvcsv); labs = ["met-d-GlycA", "ieu-b-40", "ieu-b-4877"]; names = ["GlycA", "BMI", "Smoking initiation"]
BX, SX = mv[[f"beta.{l}" for l in labs]].values, mv[[f"se.{l}" for l in labs]].values; by, sey = mv.beta_out.values, mv.se_out.values
b, se, p, sig = mv_ivw(BX, by, sey); qa, qap = q_a(BX, SX, by, sey)
s6 = pd.read_csv(os.path.join(tables, "TableS6_MVMR_results.csv")).set_index("id.exposure"); s7 = pd.read_csv(os.path.join(tables, "TableS7_MVMR_instrument_strength.csv"))
rows = []
for j, (l, nm) in enumerate(zip(labs, names)):
    f_exp = float(s7[s7.exposure.str.contains(l)].conditional_F.iloc[0])
    rows += [dict(Check=f"{nm}: causal estimate (log odds)", Reported_R=float(s6.loc[l, "b"]), Python=float(b[j])), dict(Check=f"{nm}: standard error", Reported_R=float(s6.loc[l, "se"]), Python=float(se[j])),
             dict(Check=f"{nm}: p-value", Reported_R=float(s6.loc[l, "pval"]), Python=float(p[j])), dict(Check=f"{nm}: conditional F", Reported_R=f_exp, Python=cond_F(BX, SX, j))]
rows.append(dict(Check="Q_A", Reported_R=618.54, Python=qa))
val1 = pd.DataFrame(rows); val1["Abs_diff"] = (val1.Python - val1.Reported_R).abs(); val1["Tolerance"] = [0.0001, 0.0001, 0.0001, 0.006] * 3 + [0.006]
val1["Pass"] = val1.Abs_diff <= val1.Tolerance; print(val1.to_string(), flush=True)
if not val1.Pass.all(): sys.exit("validation of the multivariable estimators failed")

# ================================================================ 0b. rebuilding that dataset from the API (checks the data-extraction code)
print("== validation 2: dataset rebuilt from the OpenGWAS API", flush=True)
dv, n_union = mvmr_dataset(labs, ["GlycA", "BMI", "Smoking_initiation"], lambda s: assoc("ebi-a-GCST006908", s), "val3")
val2 = pd.DataFrame(mvmr_summary(dv, ["GlycA", "BMI", "Smoking_initiation"], "MEGASTROKE (ebi-a-GCST006908)", "GlycA + BMI + smoking, rebuilt from the API"))
ov = len(set(dv.SNP) & set(mv.SNP)); print(f"clumped union {n_union}; usable {len(dv)}; shared with the R dataset {ov} of {len(mv)}"); print(val2[["Exposure", "N_SNPs", "Beta", "SE", "P", "Conditional_F", "Q_A"]].to_string(), flush=True)

# ================================================================ A. reverse MR
print("== A. reverse MR", flush=True)
gs, ncase, ncontrol = "GCST90104540", 62100, 1234808
def gws_instruments():
    parts = []
    for ch in pd.read_csv(os.path.join(raw, gs + ".h.tsv.gz"), sep="\t", usecols=COLS, chunksize=2_000_000, dtype=str, compression="gzip"):
        p_ = pd.to_numeric(ch.p_value, errors="coerce"); g = ch[(p_ < 5e-8) & ch.rsid.fillna("").str.startswith("rs")]
        if len(g): parts.append(g.assign(p=p_[g.index]))
    g = pd.concat(parts); g["pos"] = pd.to_numeric(g.base_pair_location); g["chr"] = pd.to_numeric(g.chromosome, errors="coerce"); g = g.sort_values("p", kind="mergesort")
    removed = np.zeros(len(g), bool); lead = []; C, P = g.chr.values, g.pos.values
    for i in range(len(g)):
        if removed[i]: continue
        lead.append(i); removed |= (C == C[i]) & (np.abs(P - P[i]) <= 500000)
    return g.iloc[lead].drop_duplicates("rsid")
g_lead = gws_instruments(); keep = clump(g_lead.rsid, g_lead.p, "gigastroke_rev"); g_inst = g_lead[g_lead.rsid.isin(keep)]
print("GIGASTROKE genome-wide significant leads:", len(g_lead), "after LD clumping:", len(g_inst), flush=True)
mg = tophits("ebi-a-GCST006908"); print("MEGASTROKE OpenGWAS top hits:", len(mg), flush=True)

def reverse(inst, label, n_case, n_ctrl):
    out = assoc("met-d-GlycA", list(inst.SNP)); neff = 4 / (1 / n_case + 1 / n_ctrl); rows = []
    for r in inst.itertuples(index=False):
        if r.SNP not in out.index: continue
        o = out.loc[r.SNP]; sg = align(dict(ea=r.ea, nea=r.nea, eaf=r.eaf), dict(ea=o.ea, nea=o.nea, eaf=o.eaf))
        if sg is None: continue
        z_e, z_o = r.beta / r.se, sg * o.beta / o.se; n_o = float(o.n) if str(o.n).strip() not in ("", "nan", "None") else N_GLYCA
        r2e, r2o = z_e**2 / (z_e**2 + neff), z_o**2 / (z_o**2 + n_o)
        rows.append(dict(Instrument_set=label, SNP=r.SNP, EA=r.ea, NEA=r.nea, EAF=r.eaf, beta_stroke=r.beta, se_stroke=r.se, p_stroke=r.p, beta_GlycA=sg * o.beta, se_GlycA=o.se,
                         p_GlycA=o.p, r2_stroke_approx=r2e, r2_GlycA_approx=r2o, Steiger_direction_correct=bool(r2e > r2o)))
    d = pd.DataFrame(rows); bx, by_, sy_, sx_ = d.beta_stroke.values, d.beta_GlycA.values, d.se_GlycA.values, d.se_stroke.values
    res = []
    for meth, fn in (("Inverse variance weighted", lambda: ivw(bx, by_, sy_)[:3]), ("MR-Egger", lambda: egger(bx, by_, sy_)[:3]), ("Weighted median", lambda: wm(bx, by_, sx_, sy_))):
        b_, se_, p_ = fn(); res.append(dict(Instrument_set=label, Method=meth, N_SNPs=len(d), Beta_SD_GlycA_per_logOR=b_, SE=se_, CI_lower=b_ - 1.96 * se_, CI_upper=b_ + 1.96 * se_, P=p_))
    res[0]["Cochran_Q_P"] = ivw(bx, by_, sy_)[3]; res[1]["Egger_intercept_P"] = egger(bx, by_, sy_)[3]
    res[0]["Steiger_correct_direction"] = f"{int(d.Steiger_direction_correct.sum())}/{len(d)}"
    return d, res
mg_i = pd.DataFrame(dict(SNP=mg.rsid, ea=mg.ea, nea=mg.nea, eaf=mg.eaf, beta=mg.beta, se=mg.se, p=mg.p))
gi_i = pd.DataFrame(dict(SNP=g_inst.rsid, ea=g_inst.effect_allele.str.upper(), nea=g_inst.other_allele.str.upper(), eaf=pd.to_numeric(g_inst.effect_allele_frequency, errors="coerce"),
                         beta=pd.to_numeric(g_inst.beta), se=pd.to_numeric(g_inst.standard_error), p=g_inst.p))
dA1, rA1 = reverse(gi_i, "GIGASTROKE ischemic stroke (GCST90104540)", ncase, ncontrol)
dA2, rA2 = reverse(mg_i, "MEGASTROKE ischemic stroke (ebi-a-GCST006908)", 34217, 406111)
revres = pd.DataFrame(rA1 + rA2); revsnp = pd.concat([dA1, dA2]); print(revres.to_string(), flush=True)

# ================================================================ B. GlycA + AF MVMR
print("== B. GlycA + AF", flush=True)
idsB, labB = ["met-d-GlycA", "ebi-a-GCST006414"], ["GlycA", "AF"]
OUTS = {"GIGASTROKE cardioembolic (GCST90104541)": "GCST90104541", "GIGASTROKE ischemic stroke (GCST90104540)": "GCST90104540", "GIGASTROKE large artery (GCST90104542)": "GCST90104542",
        "GIGASTROKE small vessel (GCST90104543)": "GCST90104543", "MEGASTROKE cardioembolic (GCST006910)": "GCST006910"}
union = pd.concat([tophits(g)[["rsid"]] for g in idsB]).rsid.unique().tolist()
loc = {k: read_local(f, union) for k, f in OUTS.items()}
resB, snpB = [], []
for k, o in loc.items():
    d, n_union = mvmr_dataset(idsB, labB, lambda s, o=o: o[o.index.isin(s)], "glyca_af"); print(k, "union after clumping", n_union, "usable", len(d), flush=True)
    resB += mvmr_summary(d, labB, k, "GlycA + AF"); snpB.append(d.assign(Outcome=k))
mvB = pd.DataFrame(resB); print(mvB[["Outcome", "Exposure", "N_SNPs", "Beta", "SE", "P", "OR", "Conditional_F", "Q_A"]].to_string(), flush=True)

# ---------------------------------------------------------------- write tables
def rnd(df):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f": df[c] = df[c].map(lambda x: float(f"{x:.4g}") if pd.notna(x) else x)
    return df
rnd(revres).to_csv(os.path.join(tables, "TableS16_reverse_MR.csv"), index=False)
rnd(revsnp).to_csv(os.path.join(tables, "TableS16b_reverse_MR_instruments.csv"), index=False)
rnd(mvB).to_csv(os.path.join(tables, "TableS17_MVMR_GlycA_AF.csv"), index=False)
rnd(pd.concat(snpB)).to_csv(os.path.join(tables, "TableS17b_MVMR_GlycA_AF_snp_level.csv"), index=False)
rnd(val1).to_csv(os.path.join(tables, "TableS17c_MVMR_code_validation.csv"), index=False)
rnd(val2).to_csv(os.path.join(tables, "TableS17d_MVMR_dataset_rebuilt_from_API.csv"), index=False)
print("tables written", flush=True)
