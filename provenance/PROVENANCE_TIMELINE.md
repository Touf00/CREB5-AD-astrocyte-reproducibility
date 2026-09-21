# Provenance timeline

Historical timestamps below come from preserved original files. Missing historical timestamps are not recreated.

- **2026-09-11 15:40:01 UTC** — `ANALYSIS_PROTOCOL_LOCK_v1.json` frozen for the GSE157827 confirmatory stage.
- **2026-09-11 23:41:26 UTC** — `GSE157827_STEP13B_PRIMARY_CONFIRMATORY_RESULT_LOCK_v1.json` created after donor-level genome-wide analysis; CREB5 was positive with nominal P < 0.05.
- **2026-09-13 17:57:46 UTC** — `EXTERNAL_VALIDATION_PROTOCOL_LOCK_v1.json` frozen. It prespecified GSE188545 → GSE138852 → GSE174367, donor-level pseudobulk, target-reveal gates, Holm correction and a three-accession meta-analysis.
- **2026-09-15 22:21:34 UTC** — an **intermediate** finalization attempt recorded GSE138852 as structurally unresolved and therefore did not execute the planned Holm/meta-analysis. That checkpoint is preserved under `historical/` and is not the final scientific state.
- **Later preserved Notebook 03 execution (cells 53–55; execution counts 49–51)** — donor resolution for GSE138852 was completed from author-released genotype-demultiplexed metadata, the final three-accession Holm correction and synthesis were executed, and a reconstruction manifest was written. The printed reconstruction-manifest SHA-256 was `bfc48bad2ed73dbcff44dfb8fd50c03c45bf8616600b4fac4fc977e9f9119af7`.

The external-validation protocol was an internal prospective analysis record, **not an externally registered preregistration**.

Authoritative release-facing records are the cleaned notebooks, `NOTEBOOK03_FINAL_EXECUTED_EVIDENCE.txt`, `RELEASE_RECONSTRUCTION_NOTE_v1.json`, and the audited final Supplementary Table S3. The files under `historical/` are retained only to make the chronology inspectable.
