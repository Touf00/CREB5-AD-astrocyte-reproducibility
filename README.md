# CREB5 upregulation is reproducible across human astrocyte datasets in Alzheimer’s disease

This repository contains reproducibility materials for the manuscript:

**CREB5 upregulation is reproducible across human astrocyte datasets in Alzheimer’s disease**

**Author:** Mohamed Tawfik, PhD  
**ORCID:** 0000-0002-5028-7463

## Study design

The primary analysis asks whether a CREB5 astrocyte signal first identified in SEA-AD middle temporal gyrus remains directionally reproducible across independent human Alzheimer’s disease datasets when donor-level inference, frozen target rules, and negative or nonsignificant results are retained.

The primary inferential sequence comprises SEA-AD, GSE160936, GSE157827, GSE188545, GSE138852, and GSE174367.

A separate secondary contextualization layer examines CREB5 in the five-region AD Progression Atlas (GSE268599) and uses the Human Protein Atlas and Agora as orthogonal normal-brain and bulk-AD contextual resources. These secondary resources are not treated as additional primary validation cohorts.

The workflow is organized into five notebooks:

1. `notebooks/01_discovery_replication_clean.ipynb` — SEA-AD discovery/regional replication and GSE160936 validation.
2. `notebooks/02_GSE157827_confirmation_clean.ipynb` — target-insulated GSE157827 confirmation and sensitivity analyses.
3. `notebooks/03_external_validation_clean.ipynb` — prespecified external-validation sequence GSE188545 → GSE138852 → GSE174367 and three-accession synthesis.
4. `notebooks/04_figures_tables_clean.ipynb` — final primary-analysis table and figure assembly from frozen results.
5. `notebooks/05_CREB5_AD_Progression_Context.ipynb` — reproducible secondary contextualization of GSE268599, HPA, and Agora; additional multiplicity calculations; Table 5; Supplementary Tables S4–S5; Supplementary Figure S6; and QA outputs.

## Key inferential rules

- The **human donor** is the biological replicate for the primary differential-expression analyses.
- Raw integer counts are aggregated to donor-level astrocyte pseudobulk before primary differential expression.
- CREB5 is not used to select QC thresholds, clustering features, donor eligibility, or the genome-wide gene universe.
- The external-validation protocol fixed accession order, eligibility, the target-reveal gate, multiplicity control, and meta-analysis before target inspection.
- Nonsignificant prespecified results are retained.
- The external-validation plan was internally frozen in the provenance record; it was **not externally preregistered**.
- GSE268599 is a separate secondary contextualization analysis and is not incorporated into the primary external-validation meta-analysis.
- The five GSE268599 regions belong to one multi-region study and are not treated as five independent cohorts.
- For GSE268599, portal-reported correlation coefficients and P values are distinguished from Benjamini–Hochberg adjustments calculated in the present study.
- HPA and Agora are used as contextual resources rather than independent astrocyte validation datasets.

## Final primary external-validation synthesis

The final three prespecified external-validation accessions are GSE188545, GSE138852, and GSE174367. GSE138852 became evaluable after exact donor reconstruction from author-released genotype-demultiplexed metadata.

The final fixed-effect inverse-variance synthesis is:

- pooled log2FC = **0.677998**
- 95% CI = **0.294169 to 1.061826**
- P = **5.36 × 10⁻4**
- I² = **58.57%**

The DerSimonian–Laird random-effects sensitivity remains positive:

- log2FC = **0.670385**
- 95% CI = **0.058999 to 1.281770**
- P = **0.0316**

These results support reproducibility in direction while retaining substantial between-cohort heterogeneity.

## Secondary AD Progression Atlas contextualization

GSE268599 was examined separately using the AD Progression Atlas sample-level portal.

Across EC, ITG, PFC, V2, and V1, astrocytic CREB5 was higher in the advanced B3 group (Braak V/VI) than in the early B1 group (Braak 0/I/II) in all five displayed regional comparisons.

An additional Benjamini–Hochberg correction across those five portal-reported B1-versus-B3 comparisons retained all five:

- EC: q = **0.0317**
- ITG: q = **4.78 × 10⁻4**
- PFC: q = **0.0104**
- V2: q = **4.78 × 10⁻4**
- V1: q = **7.90 × 10⁻5**

In contrast, none of ten portal-reported sample-level correlations between CREB5 and local Aβ plaque or pTau/total-tau burden survived an additional Benjamini–Hochberg correction across the ten tests.

These Atlas results are treated as secondary disease-stage contextualization, not as five independent replication cohorts and not as evidence that CREB5 is independent of local amyloid or tau pathology.

## HPA and Agora contextualization

Human Protein Atlas normal-brain single-nucleus information indicates that CREB5 is not a constitutive astrocyte-specific transcript; its reported brain cluster specificity is cell-type enhanced in oligodendrocytes and committed oligodendrocyte precursors.

Agora provides complementary bulk-tissue AD context. Regional RNA effects are heterogeneous rather than uniformly positive across the queried brain regions. No CREB5 SRM, TMT, or LFQ protein measurements were displayed in the queried Agora panels. Lack of displayed proteomic data is an availability limitation and is not interpreted as evidence for absence of a protein-level effect.

These resources are used only for orthogonal biological context.

## Provenance

The repository deliberately preserves the chronology of the primary analysis.

An intermediate checkpoint dated 2026-09-15 recorded GSE138852 as structurally unresolved and therefore did not execute Holm correction or the three-study synthesis. Subsequent donor reconstruction resolved the linkage, and the preserved executed Notebook 03 records the later completed three-accession analysis.

Historical intermediate files are retained under `provenance/historical/` and are explicitly labeled as superseded/intermediate. They are **not** the final analysis state.

The original GSE157827 protocol lock predates the confirmatory result lock, and the public Notebook 02 requires the historical protocol rather than generating a new timestamp.

Secondary-context QA outputs and SHA256 manifests are preserved under `provenance/`.

## Repository contents

- `notebooks/` — five cleaned/reproducible analysis notebooks.
- `provenance/` — protocol/result locks, chronology notes, historical records, secondary-context QA, and SHA256 manifests.
- `frozen_results/` — compact primary result tables plus frozen GSE268599, HPA, and Agora contextual source tables.
- `supplementary/` — machine-readable Supplementary Tables S1–S5, Supplementary Figure S6, and external-validation synthesis materials.
- `figures/` — manuscript figure source materials included in the repository release.
- `requirements.txt` and `requirements_notebook*.txt` — software environment specifications.
- `DATASET_MANIFEST.csv` — dataset and contextual-resource manifest.
- `verify_release.py` — compact consistency checks.
- `CITATION.cff` and `LICENSE` — citation metadata and repository license.

Within `supplementary/`, only `Supplementary_Table_S3_External_Validation_AUDITED_FINAL.csv` is labeled as Supplementary Table S3. The aggregate primary external-validation synthesis is stored separately as `External_Validation_Final_Synthesis_AUDITED.csv`.

New secondary-context materials include:

- `Supplementary_Table_S4_GSE268599_CREB5.csv`
- `Supplementary_Table_S4_GSE268599_CREB5.xlsx`
- `Supplementary_Table_S5_HPA_Agora_CREB5.csv`
- `Supplementary_Table_S5_HPA_Agora_CREB5.xlsx`
- `Supplementary_Figure_S6_GSE268599_CREB5.png`
- `Supplementary_Figure_S6_GSE268599_CREB5.pdf`
- `Supplementary_Figure_S6_GSE268599_CREB5.svg`

Large genome-wide result tables and archival materials may be retained in the immutable Zenodo release rather than duplicated unnecessarily in GitHub.

## Data availability

The primary biological datasets are public resources from SEA-AD/CELLxGENE and NCBI GEO accessions GSE160936, GSE157827, GSE188545, GSE138852, and GSE174367.

Secondary contextualization uses GSE268599 through the AD Progression Atlas, together with the Human Protein Atlas and Agora.

Source-dataset terms and licenses remain those of the original repositories/authors; this repository does not relicense third-party source data.

## Reproducibility

The primary compact verification workflow checks the locked GSE157827 CREB5 result and recomputes the final external fixed-effect and random-effects synthesis from accession-level estimates.

Notebook 05 separately freezes the manually transcribed GSE268599 portal outputs and HPA/Agora contextual observations, recomputes the specified Benjamini–Hochberg corrections, generates secondary tables/figures, and runs consistency checks.

The secondary-context package passed its final audit with:

- **19 required files present**
- **22 rows in Supplementary Table S4**
- **26 rows in Supplementary Table S5**
- **0/10 local-pathology correlations significant after BH correction**
- **5/5 displayed Braak B1-versus-B3 comparisons significant after the additional BH correction**

The current archived reproducibility release is published at Zenodo:

**https://doi.org/10.5281/zenodo.23002477**

The previous release remains available in the Zenodo version history.

## Citation

Please cite the associated manuscript/preprint when available.

The current reproducibility release DOI is **10.5281/zenodo.23002477**. Repository citation metadata are provided in `CITATION.cff`.
