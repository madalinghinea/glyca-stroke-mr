
"""Recompute Tables S17/S17b/S17c after the mv_multiple SE correction; cardioembolic rows come from the R run (script 11, 11 Nov 2026).
Usage: python recompute_s17.py <folder with mr_api_lib.py> <folder with the earlier Table S17b/S17c> <csv exported from data/processed/mvmr_glyca_af_ces.rds> <output folder>"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, sys.argv[1]); import mr_api_lib as L
tabs_in, r_csv, out = sys.argv[2], sys.argv[3], sys.argv[4]
os.makedirs(out, exist_ok=True)
s17b = pd.read_csv(os.path.join(tabs_in, 'TableS17b_MVMR_GlycA_AF_snp_level.csv')); labels = ['GlycA', 'AF']
CES = 'GIGASTROKE cardioembolic (GCST90104541)'
R = pd.read_csv(r_csv)
Rd = pd.DataFrame({'SNP': R.SNP, 'beta.GlycA': R['b_met-d-GlycA'], 'beta.AF': R['b_ebi-a-GCST006414'], 'se.GlycA': R['se_met-d-GlycA'], 'se.AF': R['se_ebi-a-GCST006414'],
                   'beta_out': R['by'], 'se_out': R['sey'], 'Outcome': CES})
rows, snp = [], []
for k in s17b.Outcome.unique():
    d = Rd if k == CES else s17b[s17b.Outcome == k]
    res = L.mvmr_summary(d, labels, k, 'GlycA + AF')
    for r in res: r['Source'] = 'R: script 11, TwoSampleMR 0.7.13 / MVMR 0.4.9' if k == CES else 'Python (validated against R)'
    rows += res; snp.append(d.assign(Source='R' if k == CES else 'Python'))
t17 = pd.DataFrame(rows)
# S17c: the 13 earlier checks (446-SNP set; sigma > 1) plus the R-versus-Python check on the R dataset of the CES analysis (sigma < 1)
old = pd.read_csv(os.path.join(tabs_in, 'TableS17c_MVMR_code_validation.csv'))
c = t17[t17.Outcome == CES].set_index('Exposure')
rep = {'GlycA': dict(b=-0.1343358, se=0.05855078, p=2.177058e-02, F=32.35035), 'AF': dict(b=0.6537875, se=0.02283289, p=2.566380e-180, F=71.3362)}
new = []
for ex, v in rep.items():
    nm = 'GlycA' if ex == 'GlycA' else 'AF'
    for lab, key, col, tol in (('causal estimate (log odds)', 'b', 'Beta', 1e-6), ('standard error', 'se', 'SE', 1e-6), ('conditional F', 'F', 'Conditional_F', 1e-3)):
        new.append(dict(Check=f'CES, GlycA + AF, sigma < 1: {nm}: {lab}', Reported_R=v[key], Python=c.loc[ex, col], Abs_diff=abs(v[key] - c.loc[ex, col]), Tolerance=tol, Pass=abs(v[key] - c.loc[ex, col]) <= tol))
    pv = c.loc[ex, 'P']; new.append(dict(Check=f'CES, GlycA + AF, sigma < 1: {nm}: p-value (relative difference)', Reported_R=v['p'], Python=pv, Abs_diff=abs(v['p'] - pv) / v['p'], Tolerance=1e-4, Pass=abs(v['p'] - pv) / v['p'] <= 1e-4))
new.append(dict(Check='CES, GlycA + AF: Q_A', Reported_R=108.3767, Python=c.loc['GlycA', 'Q_A'], Abs_diff=abs(108.3767 - c.loc['GlycA', 'Q_A']), Tolerance=1e-3, Pass=abs(108.3767 - c.loc['GlycA', 'Q_A']) <= 1e-3))
new.append(dict(Check='CES, GlycA + AF: Q_A p-value', Reported_R=0.9491208, Python=c.loc['GlycA', 'Q_A_P'], Abs_diff=abs(0.9491208 - c.loc['GlycA', 'Q_A_P']), Tolerance=1e-6, Pass=abs(0.9491208 - c.loc['GlycA', 'Q_A_P']) <= 1e-6))
t17c = pd.concat([old, pd.DataFrame(new)], ignore_index=True)
def rnd(df):
    df = df.copy()
    for col in df.columns:
        if df[col].dtype.kind == 'f': df[col] = df[col].map(lambda x: float('%.4g' % x) if pd.notna(x) else x)
    return df
rnd(t17).to_csv(os.path.join(out, 'TableS17_MVMR_GlycA_AF.csv'), index=False)
rnd(pd.concat(snp, ignore_index=True)[['SNP', 'EA', 'NEA', 'beta.GlycA', 'beta.AF', 'se.GlycA', 'se.AF', 'beta_out', 'se_out', 'Outcome', 'Source']]).to_csv(os.path.join(out, 'TableS17b_MVMR_GlycA_AF_snp_level.csv'), index=False)
t17c.to_csv(os.path.join(out, 'TableS17c_MVMR_code_validation.csv'), index=False)
print(rnd(t17).to_string()); print(t17c.tail(12).to_string())
