#!/usr/bin/env python3
"""Verify the compact CREB5 reproducibility release without refitting DE models."""
from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent
SUMMARY = ROOT / "frozen_results" / "EXTERNAL_VALIDATION_FINAL_RECONSTRUCTED_SUMMARY_v1.json"
SYNTHESIS = ROOT / "frozen_results" / "EXTERNAL_VALIDATION_CREB5_THREE_DATASET_SYNTHESIS_FINAL.csv"
GSE157827 = ROOT / "provenance" / "GSE157827_STEP13B_PRIMARY_CONFIRMATORY_RESULT_LOCK_v1.json"
PROTOCOL = ROOT / "provenance" / "EXTERNAL_VALIDATION_PROTOCOL_LOCK_v1.json"

for p in [SUMMARY, SYNTHESIS, GSE157827, PROTOCOL]:
    if not p.exists():
        raise FileNotFoundError(p)

protocol = json.loads(PROTOCOL.read_text())
assert protocol["scientific_rules_sha256"] == "ea54a1b3e8e2ff361ec1fc6f55ed2330410e347081926cb657f4df8892e8d8b0"
assert protocol["processing_order"] == ["GSE188545", "GSE138852", "GSE174367"]

confirm = json.loads(GSE157827.read_text())
assert confirm["CREB5_CONFIRMATORY_SUCCESS"] is True
assert np.isclose(confirm["CREB5_log2FoldChange"], 0.8566032954734822)
assert np.isclose(confirm["CREB5_nominal_pvalue"], 3.7799716195705296e-05)

x = pd.read_csv(SYNTHESIS)
assert x["accession"].tolist() == ["GSE188545", "GSE138852", "GSE174367"]
assert (x["log2FoldChange"] > 0).all()

p = x["nominal_pvalue"].to_numpy(float)
order = np.argsort(p)
adj_rank = np.empty(len(p))
running = 0.0
for rank, idx in enumerate(order):
    running = max(running, (len(p) - rank) * p[idx])
    adj_rank[idx] = min(1.0, running)
assert np.allclose(adj_rank, x["Holm_adjusted_pvalue"].to_numpy(float), rtol=1e-10, atol=1e-12)

y = x["log2FoldChange"].to_numpy(float)
se = x["SE"].to_numpy(float)
w = 1.0 / se**2
fixed = np.sum(w*y) / np.sum(w)
fixed_se = math.sqrt(1.0/np.sum(w))
fixed_p = 2*norm.sf(abs(fixed/fixed_se))
Q = float(np.sum(w*(y-fixed)**2))
df = len(y)-1
I2 = max(0.0, (Q-df)/Q)*100 if Q > 0 else 0.0
C = np.sum(w) - np.sum(w**2)/np.sum(w)
tau2 = max(0.0, (Q-df)/C)
wr = 1.0/(se**2+tau2)
random = np.sum(wr*y)/np.sum(wr)
random_se = math.sqrt(1.0/np.sum(wr))
random_p = 2*norm.sf(abs(random/random_se))

assert np.isclose(fixed, 0.677998, atol=5e-6)
assert np.isclose(fixed_se, 0.195831, atol=5e-6)
assert np.isclose(fixed_p, 0.000535861203005, rtol=1e-6)
assert np.isclose(Q, 4.826886, atol=5e-6)
assert np.isclose(I2, 58.57, atol=0.02)
assert np.isclose(tau2, 0.169902, atol=5e-6)
assert np.isclose(random, 0.670385, atol=5e-6)
assert np.isclose(random_se, 0.311931, atol=5e-6)
assert np.isclose(random_p, 0.0316230963615, rtol=1e-6)

print("CREB5 release verification: PASS")
print(f"GSE157827 confirmatory log2FC = {confirm['CREB5_log2FoldChange']:.6f}")
print(f"External fixed-effect log2FC = {fixed:.6f}, P = {fixed_p:.12g}")
print(f"External random-effects sensitivity log2FC = {random:.6f}, P = {random_p:.12g}")
print(f"I^2 = {I2:.2f}%")
