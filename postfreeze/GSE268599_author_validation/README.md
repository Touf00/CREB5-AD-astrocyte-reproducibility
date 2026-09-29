# CREB5 donor-aware secondary analysis in the MGH-AbbVie AD Progression Atlas

## Role
Post-freeze donor-aware secondary validation and neuropathological contextualization of the CREB5 astrocyte signal using author-provided sample-level average-expression data from GSE268599.

This analysis does **not** alter the prespecified three-accession external-validation family or its meta-analysis.

## Inputs
The study authors supplied:
- `ast_sample.level.average.expression_filtered.Rdata`
- associated sample metadata

The original author-provided files are not redistributed here.

## Frozen design
- Target: CREB5
- Biological unit: donor
- 32 donors, 146 donor-region observations
- 8 donors in each pathology group (1-4)
- Primary contrast: Pathology Group 1 vs Groups 3+4
- Primary test: donor-level Welch test after averaging CREB5 across available regions
- Robustness: Mann-Whitney and log1p Welch
- Sensitivity: region-adjusted linear model with donor-clustered HC1 standard errors
- Stage-specific sensitivity: Groups 2, 3 and 4 vs Group 1

## Frozen results
See `creb5_mgh_abbvie_frozen_results.tsv`.

Interpretation guardrail: these data support higher CREB5 expression in donors with greater AD neuropathologic burden; Group 2 was heterogeneous, so the result should not be described as a strictly monotonic four-stage trajectory.
