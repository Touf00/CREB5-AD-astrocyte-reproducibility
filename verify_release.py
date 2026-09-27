#!/usr/bin/env python3
"""Verify the CREB5 reproducibility release without refitting DE models."""

from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent

# ============================================================
# PRIMARY ANALYSIS FILES
# ============================================================

SUMMARY = (
    ROOT / "frozen_results" /
    "EXTERNAL_VALIDATION_FINAL_RECONSTRUCTED_SUMMARY_v1.json"
)

SYNTHESIS = (
    ROOT / "frozen_results" /
    "EXTERNAL_VALIDATION_CREB5_THREE_DATASET_SYNTHESIS_FINAL.csv"
)

GSE157827 = (
    ROOT / "provenance" /
    "GSE157827_STEP13B_PRIMARY_CONFIRMATORY_RESULT_LOCK_v1.json"
)

PROTOCOL = (
    ROOT / "provenance" /
    "EXTERNAL_VALIDATION_PROTOCOL_LOCK_v1.json"
)

# ============================================================
# SECONDARY CONTEXT FILES
# ============================================================

TABLE5 = (
    ROOT / "frozen_results" /
    "Table5_GSE268599_CREB5_context.csv"
)

LOCAL = (
    ROOT / "frozen_results" /
    "gse268599_creb5_local_pathology.csv"
)

BRAAK = (
    ROOT / "frozen_results" /
    "gse268599_creb5_braak.csv"
)

HPA = (
    ROOT / "frozen_results" /
    "hpa_creb5_context.csv"
)

AGORA_EXPR = (
    ROOT / "frozen_results" /
    "agora_creb5_overall_expression.csv"
)

S4 = (
    ROOT / "supplementary" /
    "Supplementary_Table_S4_GSE268599_CREB5.csv"
)

S5 = (
    ROOT / "supplementary" /
    "Supplementary_Table_S5_HPA_Agora_CREB5.csv"
)

QA = (
    ROOT / "provenance" /
    "QA_gse268599_creb5_context.json"
)

required = [
    SUMMARY,
    SYNTHESIS,
    GSE157827,
    PROTOCOL,
    TABLE5,
    LOCAL,
    BRAAK,
    HPA,
    AGORA_EXPR,
    S4,
    S5,
    QA,
]

for p in required:
    if not p.exists():
        raise FileNotFoundError(p)

# ============================================================
# PRIMARY ANALYSIS VERIFICATION
# ============================================================

protocol = json.loads(PROTOCOL.read_text())

assert protocol["scientific_rules_sha256"] == (
    "ea54a1b3e8e2ff361ec1fc6f55ed2330410e347081926cb657f4df8892e8d8b0"
)

assert protocol["processing_order"] == [
    "GSE188545",
    "GSE138852",
    "GSE174367"
]

confirm = json.loads(GSE157827.read_text())

assert confirm["CREB5_CONFIRMATORY_SUCCESS"] is True

assert np.isclose(
    confirm["CREB5_log2FoldChange"],
    0.8566032954734822
)

assert np.isclose(
    confirm["CREB5_nominal_pvalue"],
    3.7799716195705296e-05
)

x = pd.read_csv(SYNTHESIS)

assert x["accession"].tolist() == [
    "GSE188545",
    "GSE138852",
    "GSE174367"
]

assert (x["log2FoldChange"] > 0).all()

# Holm correction
p = x["nominal_pvalue"].to_numpy(float)

order = np.argsort(p)
adj_rank = np.empty(len(p))
running = 0.0

for rank, idx in enumerate(order):
    running = max(
        running,
        (len(p) - rank) * p[idx]
    )
    adj_rank[idx] = min(1.0, running)

assert np.allclose(
    adj_rank,
    x["Holm_adjusted_pvalue"].to_numpy(float),
    rtol=1e-10,
    atol=1e-12
)

# Fixed-effect synthesis
y = x["log2FoldChange"].to_numpy(float)
se = x["SE"].to_numpy(float)

w = 1.0 / se**2

fixed = np.sum(w * y) / np.sum(w)
fixed_se = math.sqrt(1.0 / np.sum(w))
fixed_p = 2 * norm.sf(abs(fixed / fixed_se))

Q = float(np.sum(w * (y - fixed)**2))
df = len(y) - 1

I2 = (
    max(0.0, (Q - df) / Q) * 100
    if Q > 0 else 0.0
)

C = (
    np.sum(w)
    - np.sum(w**2) / np.sum(w)
)

tau2 = max(
    0.0,
    (Q - df) / C
)

wr = 1.0 / (se**2 + tau2)

random = np.sum(wr * y) / np.sum(wr)
random_se = math.sqrt(1.0 / np.sum(wr))
random_p = 2 * norm.sf(abs(random / random_se))

assert np.isclose(fixed, 0.677998, atol=5e-6)
assert np.isclose(fixed_se, 0.195831, atol=5e-6)
assert np.isclose(
    fixed_p,
    0.000535861203005,
    rtol=1e-6
)

assert np.isclose(Q, 4.826886, atol=5e-6)
assert np.isclose(I2, 58.57, atol=0.02)
assert np.isclose(tau2, 0.169902, atol=5e-6)
assert np.isclose(random, 0.670385, atol=5e-6)
assert np.isclose(random_se, 0.311931, atol=5e-6)

assert np.isclose(
    random_p,
    0.0316230963615,
    rtol=1e-6
)

# ============================================================
# SECONDARY CONTEXT VERIFICATION
# ============================================================

table5 = pd.read_csv(TABLE5)
local = pd.read_csv(LOCAL)
braak = pd.read_csv(BRAAK)
s4 = pd.read_csv(S4)
s5 = pd.read_csv(S5)

expected_regions = [
    "EC",
    "ITG",
    "PFC",
    "V2",
    "V1"
]

assert table5["Region"].tolist() == expected_regions
assert local["region"].tolist() == expected_regions
assert braak["region"].tolist() == expected_regions

assert len(table5) == 5
assert len(s4) == 22
assert len(s5) == 26

# Local pathology: none survives BH
assert (
    table5["Abeta_BH_q"] >= 0.05
).all()

assert (
    table5["pTauTau_BH_q"] >= 0.05
).all()

# Braak B1 vs B3: all five survive additional BH
assert (
    table5[
        "Braak_B1_vs_B3_additional_BH_q"
    ] < 0.05
).all()

assert (
    table5["Braak_direction"] == "B3 > B1"
).all()

expected_braak_q = np.array([
    0.031700,
    0.000478,
    0.0103875,
    0.000478,
    0.000079
])

assert np.allclose(
    table5[
        "Braak_B1_vs_B3_additional_BH_q"
    ].to_numpy(float),
    expected_braak_q,
    rtol=0.01,
    atol=1e-7
)

# HPA contextualization present
hpa = pd.read_csv(HPA)

assert len(hpa) == 1
assert "Oligodendrocyte" in str(
    hpa.loc[0, "Value"]
)

# Agora overall-expression table
agora = pd.read_csv(AGORA_EXPR)

assert len(agora) == 9

assert set(agora["Region"]) == {
    "ACC",
    "CBE",
    "DLPFC",
    "FP",
    "IFG",
    "PCC",
    "PHG",
    "STG",
    "TCX"
}

# Secondary QA report
qa = json.loads(QA.read_text())

assert qa["status"] == "PASS"
assert qa["local_pathology_q_lt_0.05"] == 0
assert qa["braak_B1_vs_B3_q_lt_0.05"] == 5

# ============================================================
# FINAL REPORT
# ============================================================

print("========================================")
print("CREB5 RELEASE VERIFICATION: PASS")
print("========================================")

print()
print("PRIMARY ANALYSIS")
print(
    "GSE157827 confirmatory log2FC = "
    f"{confirm['CREB5_log2FoldChange']:.6f}"
)
print(
    "External fixed-effect log2FC = "
    f"{fixed:.6f}, P = {fixed_p:.12g}"
)
print(
    "External random-effects log2FC = "
    f"{random:.6f}, P = {random_p:.12g}"
)
print(f"I^2 = {I2:.2f}%")

print()
print("SECONDARY CONTEXT")
print("GSE268599 regions = 5")
print("Local pathology surviving BH = 0 / 10")
print("Braak B1-vs-B3 surviving BH = 5 / 5")
print("Supplementary Table S4 rows = 22")
print("Supplementary Table S5 rows = 26")

print()
print("ALL CHECKS PASSED")
