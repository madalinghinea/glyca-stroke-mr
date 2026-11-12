"""Independent check of the data-handling steps against the raw GWAS Catalog files (no R, no TwoSampleMR).

 1. Allele harmonisation of the GlycA instruments against GIGASTROKE CES / LAS / SVS (exact rsID matching, TwoSampleMR action = 2
    logic) versus Table S9e.
 2. Atrial-fibrillation lead SNPs (p < 5e-8, greedy pruning by chi-square within 500 kb), AF -> stroke IVW estimates versus
    Table S14, and the cardioembolic IVW estimate after removing AF-related GlycA instruments versus Table S11.

Usage:  python independent_raw_check.py <project_root> <raw_dir>        (about 5 minutes; reads ~2.5 GB of raw files)
Exit status 1 if any comparison fails.
"""
import sys, os, re
import numpy as np, pandas as pd
from scipy.stats import norm

root = sys.argv[1] if len(sys.argv) > 1 else "."
raw = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "data", "raw")
T = lambda f: pd.read_csv(os.path.join(root, "tables", f))
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
OUT_COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"]

def scan(path, key, ids, usecols):
    out = []
    for ch in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=2_000_000, dtype=str, compression="gzip"):
        s = ch[ch[key].isin(ids)]
        if len(s): out.append(s)
    return pd.concat(out)

def harmonise(exp, out, tol=0.08):
    """TwoSampleMR harmonise_data(action = 2): align to the exposure effect allele; palindromic SNPs by allele frequency; drop ambiguous."""
    m = exp.merge(out, on="SNP", how="inner"); rows = []
    for r in m.itertuples(index=False):
        a1e, a2e, a1o, a2o = r.ea_e, r.oa_e, r.ea_o, r.oa_o
        bo = float(r.beta_o); eo = float(r.eaf_o) if r.eaf_o not in (None, "NA", "") and pd.notna(r.eaf_o) else np.nan
        ee = r.eaf_e; pal = (a1e, a2e) in (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")); keep, flip = True, False
        if not pal:
            if (a1e, a2e) == (a1o, a2o) or (COMP.get(a1e), COMP.get(a2e)) == (a1o, a2o): pass
            elif (a1e, a2e) == (a2o, a1o) or (COMP.get(a1e), COMP.get(a2e)) == (a2o, a1o): flip = True
            else: keep = False
        elif np.isnan(eo) or abs(ee - 0.5) < tol or abs(eo - 0.5) < tol: keep = False
        else: flip = (ee - 0.5) * (eo - 0.5) < 0
        rows.append(dict(SNP=r.SNP, keep=keep, beta_exp=r.beta_e, beta_out=-bo if flip else bo, se_out=float(r.se_o)))
    h = pd.DataFrame(rows); return h[h.keep]

def ivw(bx, by, sey):
    w = 1 / sey**2; b = np.sum(w * bx * by) / np.sum(w * bx**2); q = np.sum(w * (by - b * bx)**2)
    se = np.sqrt(1 / np.sum(w * bx**2)) * max(np.sqrt(q / (len(bx) - 1)), 1.0)
    return b, se, 2 * norm.sf(abs(b / se))

res = []                                          # (check, quantity, python, reported, tolerance)
s1 = T("TableS1_GlycA_instruments.csv")
exp = pd.DataFrame(dict(SNP=s1.SNP, ea_e=s1.effect_allele, oa_e=s1.other_allele, eaf_e=s1.eaf, beta_e=s1.beta_SD, chr=s1.chr.astype(int), pos=s1.pos_hg19.astype(float)))

# ---- 1. AF lead SNPs and GlycA lookups in the AF file
af_cols = ["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "hm_effect_allele_frequency", "standard_error", "p_value"]
gws, gly = [], set(exp.SNP)
for ch in pd.read_csv(os.path.join(raw, "GCST006414.h.tsv.gz"), sep="\t", usecols=af_cols, chunksize=2_000_000, dtype=str, compression="gzip"):
    p = pd.to_numeric(ch.p_value, errors="coerce")
    g = ch[(p < 5e-8) & ch.hm_rsid.fillna("").str.startswith("rs") & ch.hm_beta.notna() & ch.standard_error.notna()]
    if len(g): gws.append(g.assign(p=p[g.index], _glyca=False))
    s = ch[ch.hm_rsid.isin(gly)]
    if len(s): gws.append(s.assign(p=pd.to_numeric(s.p_value, errors="coerce"), _glyca=True))
allr = pd.concat(gws); allr["_glyca"] = allr["_glyca"].astype(bool); glyca_af = allr[allr._glyca].drop_duplicates("hm_rsid")
a = allr[~allr._glyca].copy()
a["chr"] = pd.to_numeric(a.hm_chrom, errors="coerce"); a = a[a.chr.notna()].copy(); a["chr"] = a.chr.astype(int)
a["pos"] = a.hm_pos.astype(float); a["beta"] = a.hm_beta.astype(float); a["se"] = a.standard_error.astype(float); a["z2"] = (a.beta / a.se)**2
a = a.sort_values("z2", ascending=False, kind="mergesort").reset_index(drop=True)
removed = np.zeros(len(a), bool); lead = []
for i in range(len(a)):
    if removed[i]: continue
    lead.append(i); removed |= (a.chr.values == a.chr.values[i]) & (np.abs(a.pos.values - a.pos.values[i]) <= 500000)
af = a.iloc[lead].drop_duplicates("hm_rsid")
s14 = T("TableS14_atrial_fibrillation_analyses.csv")
n_rep = int(re.search(r"(\d+) positionally", s14.Exposure.iloc[0]).group(1))
res.append(("AF leads", "number of AF lead SNPs", len(af), n_rep, 0))

# ---- 2. outcome files: GlycA instruments (S9e) and AF leads (S14)
ids = gly | set(af.hm_rsid)
out = {}
for k, f in (("GIGASTROKE", "GCST90104540"), ("CES", "GCST90104541"), ("LAS", "GCST90104542"), ("SVS", "GCST90104543")):
    o = scan(os.path.join(raw, f + ".h.tsv.gz"), "rsid", ids, OUT_COLS).drop_duplicates("rsid")
    out[k] = pd.DataFrame(dict(SNP=o.rsid, ea_o=o.effect_allele, oa_o=o.other_allele, beta_o=o.beta, se_o=o.standard_error, eaf_o=o.effect_allele_frequency))
s9e = T("TableS9e_subtype_harmonised_with_steiger.csv")
s3 = T("TableS3_GIGASTROKE_harmonised.csv"); s3 = s3[s3.mr_keep.astype(str).str.upper() == "TRUE"].assign(subtype="GIGASTROKE")
s9e = pd.concat([s9e, s3[s9e.columns.intersection(s3.columns)]], ignore_index=True)
for k in ("GIGASTROKE", "CES", "LAS", "SVS"):
    h = harmonise(exp, out[k]); sv = s9e[s9e.subtype == k].merge(h, on="SNP", how="outer", indicator=True)
    both = sv[sv._merge == "both"]
    res += [(f"{k} harmonisation", "SNPs retained", len(h), int((s9e.subtype == k).sum()), 0),
            (f"{k} harmonisation", "SNPs in one set only", int((sv._merge != "both").sum()), 0, 0),
            (f"{k} harmonisation", "max |beta.outcome difference|", float((both.beta_out - both["beta.outcome"]).abs().max()), 0.0, 1e-4)]
afx = pd.DataFrame(dict(SNP=af.hm_rsid, ea_e=af.hm_effect_allele, oa_e=af.hm_other_allele, eaf_e=af.hm_effect_allele_frequency.astype(float), beta_e=af.beta))
for k in ("CES", "LAS", "SVS"):
    h = harmonise(afx, out[k]); b, se, p = ivw(h.beta_exp.values, h.beta_out.values, h.se_out.values)
    r = s14[(s14.Analysis.str.startswith("AF")) & (s14.Outcome == k) & (s14.Method == "Inverse variance weighted")].iloc[0]
    res += [(f"AF -> {k}", "SNPs usable", len(h), int(r.N_SNPs), 0), (f"AF -> {k}", "IVW OR", float(np.exp(b)), float(r.OR), 6e-4), (f"AF -> {k}", "IVW p (log10)", float(np.log10(p)), float(np.log10(r.P)), 0.01)]

# ---- 3. cardioembolic estimate after removing AF-related GlycA instruments (Table S11)
lead_xy = af[["chr", "pos"]].values
near = {r.SNP for r in exp.itertuples() if np.any((lead_xy[:, 0] == r.chr) & (np.abs(lead_xy[:, 1] - r.pos) <= 5e5))}
g = glyca_af.set_index("hm_rsid"); low = {s for s in exp.SNP if s in g.index and g.loc[s, "p"] < 1e-3}
ces = s9e[s9e.subtype == "CES"]; s11 = T("TableS11_robustness_analysis_sets.csv")
for name, drop, label in (("AF p < 1e-3", low, "Excluding SNPs with AF p < 1e-3"), ("near AF lead", near, "Excluding SNPs within 500 kb of an AF lead SNP"), ("both", low | near, "Excluding both")):
    d = ces[~ces.SNP.isin(drop)]; b, se, p = ivw(d["beta.exposure"].values, d["beta.outcome"].values, d["se.outcome"].values)
    r = s11[(s11.Outcome == "GIGASTROKE CES (GCST90104541)") & (s11.Analysis_set == label) & (s11.Method == "IVW")].iloc[0]
    res += [(f"CES excl. {name}", "SNPs", len(d), int(r.N_SNPs), 0), (f"CES excl. {name}", "IVW OR", float(np.exp(b)), float(r.OR), 6e-4), (f"CES excl. {name}", "IVW p", float(p), float(r.P), 6e-4)]

df = pd.DataFrame(res, columns=["check", "quantity", "python", "reported", "tolerance"])
df["abs_diff"] = (df.python - df.reported).abs(); df["pass"] = df.abs_diff <= df.tolerance
pd.set_option("display.width", 200); pd.set_option("display.float_format", lambda x: f"{x:.5g}")
print(df.to_string(index=False)); df.to_csv(os.path.join(root, "tables", "independent_raw_check_results.csv"), index=False)
print("\nFAILED:", int((~df["pass"]).sum()), "of", len(df)); sys.exit(1 if (~df["pass"]).any() else 0)
