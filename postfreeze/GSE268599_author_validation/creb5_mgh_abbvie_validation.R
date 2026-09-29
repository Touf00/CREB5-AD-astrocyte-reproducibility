#!/usr/bin/env Rscript

# CREB5 donor-aware secondary analysis in the MGH-AbbVie AD Progression Atlas
# Inputs supplied by the study authors:
#   ast_sample.level.average.expression_filtered.Rdata
#   AD_progression_meta.xlsx
#
# Primary unit of inference: donor (NOT donor-region rows).
# Primary contrast: Pathology Group 1 vs Groups 3+4 (high AD neuropathologic burden).
# Sensitivity model: region-adjusted linear model with donor-clustered HC1 SEs.

options(stringsAsFactors = FALSE)

rdata_file <- "ast_sample.level.average.expression_filtered.Rdata"
meta_file  <- "AD_progression_meta.xlsx"
out_dir    <- "creb5_mgh_abbvie_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_file)) stop("Missing: ", rdata_file)
if (!file.exists(meta_file)) stop("Missing: ", meta_file)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to read AD_progression_meta.xlsx")
}

env <- new.env(parent = emptyenv())
loaded <- load(rdata_file, envir = env)

obj_name <- NULL
for (nm in loaded) {
  x <- get(nm, envir = env)
  if (is.data.frame(x) && all(c("gene", "ave.exp", "Donor.ID") %in% names(x))) {
    obj_name <- nm
    break
  }
}
if (is.null(obj_name)) {
  stop("Could not locate an expression object with columns gene, ave.exp, Donor.ID. Loaded objects: ", paste(loaded, collapse=", "))
}
expr <- get(obj_name, envir = env)
expr$Donor.ID <- as.character(expr$Donor.ID)

meta <- readxl::read_excel(meta_file)
meta <- as.data.frame(meta)
meta$Donor.ID <- as.character(meta$Donor.ID)
meta$Path..Group. <- as.integer(meta$Path..Group.)

region_map <- c(EC="EC", BA20="ITG", BA46="PFC", V2="V2", V1="V1",
                ITG="ITG", PFC="PFC")
if ("Region" %in% names(expr)) {
  expr$Region_std <- unname(region_map[as.character(expr$Region)])
} else if ("Unified_region" %in% names(expr)) {
  expr$Region_std <- unname(region_map[as.character(expr$Unified_region)])
} else {
  expr$Region_std <- NA_character_
}
meta$Region_std <- unname(region_map[as.character(meta$Region)])

donor_meta <- unique(meta[, c("Donor.ID", "Path..Group.")])
if (anyDuplicated(donor_meta$Donor.ID)) stop("Pathology group is not unique per donor.")
expr <- merge(expr, donor_meta, by="Donor.ID", all.x=TRUE, sort=FALSE)

if (all(is.na(expr$Region_std)) && "donor_region" %in% names(expr)) {
  rr <- sub("^[^_]+_", "", as.character(expr$donor_region))
  expr$Region_std <- unname(region_map[rr])
}

creb <- expr[expr$gene == "CREB5", , drop=FALSE]
if (nrow(creb) == 0) stop("CREB5 not found in expression object.")

checks <- data.frame(
  metric = c("loaded_object", "expression_rows", "expression_genes", "CREB5_rows",
             "CREB5_donors", "metadata_rows", "metadata_donors", "metadata_samples"),
  value = c(obj_name, nrow(expr), length(unique(expr$gene)), nrow(creb),
            length(unique(creb$Donor.ID)), nrow(meta), length(unique(meta$Donor.ID)),
            length(unique(meta$sample_ID)))
)
write.csv(checks, file.path(out_dir, "integrity_checks.csv"), row.names=FALSE)

donor <- aggregate(ave.exp ~ Donor.ID + Path..Group., data=creb, FUN=mean)
donor$contrast <- ifelse(donor$Path..Group. == 1, "Group1_low",
                         ifelse(donor$Path..Group. %in% c(3,4), "Groups3_4_high", "Group2_intermediate"))

primary <- donor[donor$Path..Group. %in% c(1,3,4), ]
primary$high_adnc <- as.integer(primary$Path..Group. %in% c(3,4))
low  <- primary$ave.exp[primary$high_adnc == 0]
high <- primary$ave.exp[primary$high_adnc == 1]

welch <- t.test(high, low, var.equal=FALSE)
wilcx <- wilcox.test(high, low, exact=FALSE)
welch_log <- t.test(log1p(high), log1p(low), var.equal=FALSE)

primary_summary <- data.frame(
  analysis = c("donor_level_raw", "donor_level_raw", "donor_level_raw", "donor_level_log1p"),
  statistic = c("mean_Group1", "mean_Groups3_4", "Welch_p", "Welch_p"),
  value = c(mean(low), mean(high), welch$p.value, welch_log$p.value)
)
primary_summary <- rbind(primary_summary,
  data.frame(analysis="donor_level_raw", statistic="Mann_Whitney_p", value=wilcx$p.value),
  data.frame(analysis="donor_level_raw", statistic="fold_high_over_low", value=mean(high)/mean(low)),
  data.frame(analysis="donor_level_raw", statistic="n_Group1_donors", value=length(low)),
  data.frame(analysis="donor_level_raw", statistic="n_Groups3_4_donors", value=length(high))
)
write.csv(primary_summary, file.path(out_dir, "primary_donor_level_results.csv"), row.names=FALSE)
write.csv(donor, file.path(out_dir, "CREB5_donor_means.csv"), row.names=FALSE)

cluster_vcov_hc1 <- function(model, cluster) {
  X <- model.matrix(model)
  u <- residuals(model)
  cluster <- as.factor(cluster)
  N <- nrow(X); K <- ncol(X); G <- nlevels(cluster)
  bread <- solve(crossprod(X))
  meat <- matrix(0, nrow=K, ncol=K)
  for (g in levels(cluster)) {
    idx <- which(cluster == g)
    Xg <- X[idx, , drop=FALSE]
    ug <- matrix(u[idx], ncol=1)
    sg <- crossprod(Xg, ug)
    meat <- meat + sg %*% t(sg)
  }
  correction <- (G/(G-1)) * ((N-1)/(N-K))
  correction * bread %*% meat %*% bread
}

coef_cluster <- function(model, cluster, term) {
  V <- cluster_vcov_hc1(model, cluster)
  b <- coef(model)[term]
  se <- sqrt(diag(V))[term]
  z <- b/se
  p <- 2*pnorm(abs(z), lower.tail=FALSE)
  ci <- b + c(-1,1)*qnorm(0.975)*se
  data.frame(term=term, estimate=b, SE=se, CI_low=ci[1], CI_high=ci[2], p_value=p)
}

model_dat <- creb[creb$Path..Group. %in% c(1,3,4) & !is.na(creb$Region_std), , drop=FALSE]
model_dat$high_adnc <- as.integer(model_dat$Path..Group. %in% c(3,4))
model_dat$Region_std <- factor(model_dat$Region_std)

m_raw <- lm(ave.exp ~ high_adnc + Region_std, data=model_dat)
r_raw <- coef_cluster(m_raw, model_dat$Donor.ID, "high_adnc")
r_raw$analysis <- "region_adjusted_clustered_raw"

model_dat$log1p_expr <- log1p(model_dat$ave.exp)
m_log <- lm(log1p_expr ~ high_adnc + Region_std, data=model_dat)
r_log <- coef_cluster(m_log, model_dat$Donor.ID, "high_adnc")
r_log$analysis <- "region_adjusted_clustered_log1p"

stage_dat <- creb[creb$Path..Group. %in% c(1,2,3,4) & !is.na(creb$Region_std), , drop=FALSE]
stage_dat$stage <- factor(stage_dat$Path..Group., levels=c(1,2,3,4))
stage_dat$Region_std <- factor(stage_dat$Region_std)
m_stage <- lm(ave.exp ~ stage + Region_std, data=stage_dat)
r_stage <- do.call(rbind, lapply(c("stage2","stage3","stage4"), function(tt) {
  x <- coef_cluster(m_stage, stage_dat$Donor.ID, tt)
  x$analysis <- "region_adjusted_clustered_stage_specific"
  x
}))

clustered <- rbind(r_raw, r_log, r_stage)
write.csv(clustered, file.path(out_dir, "clustered_region_adjusted_results.csv"), row.names=FALSE)

creb$high_group <- ifelse(creb$Path..Group. %in% c(3,4), "Groups3_4_high",
                          ifelse(creb$Path..Group.==1, "Group1_low", "Group2_intermediate"))
region_summary <- aggregate(ave.exp ~ Region_std + high_group, data=creb, FUN=function(x) c(n=length(x), mean=mean(x), sd=sd(x)))
region_out <- data.frame(Region_std=region_summary$Region_std,
                         high_group=region_summary$high_group,
                         n=region_summary$ave.exp[,"n"],
                         mean=region_summary$ave.exp[,"mean"],
                         sd=region_summary$ave.exp[,"sd"])
write.csv(region_out, file.path(out_dir, "region_descriptive_summary.csv"), row.names=FALSE)
