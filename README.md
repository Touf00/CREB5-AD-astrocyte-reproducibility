# CREB5 upregulation in human Alzheimer’s disease astrocytes — reproducibility materials

This repository contains the reproducibility materials for the manuscript:

**Sequential cross-cohort validation supports reproducible CREB5 upregulation in human Alzheimer’s disease astrocytes**

**Author:** Mohamed Tawfik, PhD  
**ORCID:** 0000-0002-5028-7463

## Study design

The analysis asks whether a CREB5 astrocyte signal first identified in SEA-AD middle temporal gyrus remains directionally reproducible across independent human Alzheimer’s disease datasets when donor-level inference, frozen target rules, and negative/nonsignificant results are retained.

The workflow is organized into four notebooks:

1. `notebooks/01_discovery_replication_clean.ipynb` — SEA-AD discovery/regional replication and GSE160936 validation.
2. `notebooks/02_GSE157827_confirmation_clean.ipynb` — target-insulated GSE157827 confirmation and sensitivity analyses.
3. `notebooks/03_external_validation_clean.ipynb` — the prespecified external-validation sequence GSE188545 → GSE138852 → GSE174367 and three-accession synthesis.
4. `notebooks/04_figures_tables_clean.ipynb` — final table and figure assembly from frozen results.

## Key inferential rules

- The **human donor** is the biological replicate.
- Raw integer counts are aggregated to donor-level astrocyte pseudobulk before differential expression.
- CREB5 is not used to select QC thresholds, clustering features, donor eligibility, or the genome-wide gene universe.
- The external-validation protocol fixed accession order, eligibility, the target-reveal gate, multiplicity control, and meta-analysis before target inspection.
- Nonsignificant prespecified results are retained.
- The external-validation plan was internally frozen in the provenance record; it was **not externally preregistered**.

## Final external-validation synthesis

The final three prespecified accessions are GSE188545, GSE138852, and GSE174367. GSE138852 became evaluable after exact donor reconstruction from author-released genotype-demultiplexed metadata.

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

## Provenance

The repository deliberately preserves the chronology of the analysis.

An intermediate checkpoint dated 2026-09-15 recorded GSE138852 as structurally unresolved and therefore did not execute Holm correction or the three-study synthesis. Subsequent donor reconstruction resolved the linkage, and the preserved executed Notebook 03 records the later completed three-accession analysis.

Historical intermediate files are retained under `provenance/historical/` and are explicitly labeled as superseded/intermediate. They are **not** the final analysis state.

The original GSE157827 protocol lock predates the confirmatory result lock, and the public Notebook 02 requires the historical protocol rather than generating a new timestamp.

## Repository contents

- `notebooks/` — cleaned public notebooks.
- `provenance/` — protocol/result locks, chronology notes, and historical intermediate records.
- `frozen_results/` — compact result tables and machine-readable final summaries needed to verify the manuscript statistics.
- `supplementary/` — final supplementary tables.
- `figures/` — manuscript figures when included in the archived release.
- `environment.yml` and `requirements.txt` — software environment information.
- `verify_release.py` — compact consistency checks.
- `CHECKSUMS_SHA256.txt` / `RELEASE_MANIFEST.tsv` — file-level release audit.

Large genome-wide result tables and executed archival notebooks are intended for the immutable Zenodo release rather than duplicated unnecessarily in GitHub.

## Data availability

All biological datasets analyzed are public resources. The study uses SEA-AD/CELLxGENE and NCBI GEO accessions GSE160936, GSE157827, GSE188545, GSE138852, and GSE174367. Source-dataset terms and licenses remain those of the original repositories/authors; this repository does not relicense third-party source data.

## Reproducibility

The compact verification script checks the locked GSE157827 CREB5 result and recomputes the final external fixed-effect synthesis from the accession-level estimates.

A frozen archival release is reserved at Zenodo: **https://doi.org/10.5281/zenodo.22884235**.

## Citation

Please cite the associated manuscript/preprint when available. The reproducibility release DOI is **10.5281/zenodo.22884235**. Repository citation metadata are provided in `CITATION.cff`.
