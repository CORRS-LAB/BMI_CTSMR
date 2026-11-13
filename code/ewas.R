# Set the working directory to the project root before running this script.
if (!dir.exists("data")) stop("Set the working directory to the project root containing data.")
dir.create("results", recursive=TRUE, showWarnings=FALSE)
# ---- 1. Discovery BMI model ----
discovery_inputs <- c("data/discovery_model_phenotype.csv",
  "data/discovery_cell_proportions.csv",
  "data/discovery.pre_QC.beta.tsv")
if (!all(file.exists(discovery_inputs))) {
  stop("Download the discovery phenotype and methylation matrix from GEO using the link in the paper, and place them in data/ as described in README.md.")
}
inp <- read.csv("data/discovery_model_phenotype.csv")
prop.est <- read.csv("data/discovery_cell_proportions.csv",
                     check.names = FALSE)
mr.out <- data.table::fread(
  "data/discovery.pre_QC.beta.tsv",
  check.names = FALSE)
# Preserve the pre-exclusion cohort size used by the original zero-value QC.
discovery_qc_reference_n <- nrow(inp)
selected_ids <- inp$sample_id[inp$sample_id != "S6"]
new_ids <- paste0("S", seq_along(selected_ids))
inp <- inp[inp$sample_id != "S6", ]
inp$sample_id <- new_ids
dat <- as.data.frame(mr.out[, c("chr", "start", "end", selected_ids),
                            with = FALSE])
names(dat)[-(1:3)] <- new_ids
data <- dat[rowSums(dat[, -(1:3)] == 0, na.rm = TRUE) < discovery_qc_reference_n * 0.2, ]
prop.est <- prop.est[match(selected_ids, prop.est$sample_id), ]
cells <- names(prop.est)[-1]
prop <- data.matrix(prop.est[, -1])
phy.run <- data.frame(base_bmi = inp$base_bmi)
covs <- inp[, c("base_age", "base_gender")]
data.no.na <- na.omit(data)
mr.no.na <- data.matrix(data.no.na[, -(1:3)])
      alpha.all <- pv.all <- fdr.all <- NULL
      for(j in 1:ncol(phy.run)){
         alpha.out <- pv.out <- NULL     
          y = phy.run[,j]
          na.ids <- which(is.na(y))
          if(length(na.ids) > 0){
            DummyBeta.m <- mr.no.na[,-na.ids]
            pheno.v <- y[-na.ids]
            frac.m <- prop[-na.ids,]
            covs.m <- covs[-na.ids,]
          }else{
            DummyBeta.m <- mr.no.na
            pheno.v <- y
            frac.m <- prop  
            covs.m <- covs
          }
        for(r in 1:nrow(DummyBeta.m)){
           beta.v <- DummyBeta.m[r,]
           design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
           cov.mod <- model.matrix(~ base_age + factor(base_gender), data=covs.m)[,-1]
           design <- cbind(design, cov.mod)
           fit <- summary(lm(beta.v ~ .-1, data = data.frame(design),na.action=na.omit))
           alpha <- as.numeric(fit$coef[14:26,1])
           alpha.out <- rbind(alpha.out,alpha)
           pv <- as.numeric(fit$coef[14:26,4])
           pv.out <- rbind(pv.out, pv)
        }
        fdr.out <- apply(pv.out, 2, function(p) p.adjust(p, method = "BH"))
        alpha.all <- rbind(alpha.all, data.frame(data.no.na[,1:3],alpha.out, phenotype = colnames(phy.run)[j]))
        pv.all <- rbind(pv.all,data.frame(data.no.na[1:3],pv.out, phenotype = colnames(phy.run)[j]))
        fdr.all <- rbind(fdr.all,data.frame(data.no.na[1:3],fdr.out, phenotype = colnames(phy.run)[j]))
      }
   colnames(alpha.all)[c(4:16)] <- paste(cells, ".alpha", sep = "")
   colnames(pv.all)[c(4:16)] <- paste(cells, ".pv", sep = "")
   colnames(fdr.all)[c(4:16)] <- paste(cells, ".fdr", sep = "")
library(plyr)   
out <- join(join(alpha.all[,c(1:16)], pv.all[,c(1:16)], by = c("chr", "start", "end")), 
                       fdr.all[,c(1:16)],by = c("chr", "start", "end"))
BMI.discovery.all <- out
BMI.discovery.sig <- out[rowSums(out[,c(30:42)] < 0.05) >=1,]
# ---- 2. Validation BMI model ----
inp <- read.csv("data/validation_model_phenotype.csv")
prop.est <- read.csv("data/validation_cell_proportions.csv",
                     check.names = FALSE)
prop.est <- prop.est[match(inp$sample_id, prop.est$sample_id), ]
cells <- names(prop.est)[-1]
prop <- data.matrix(prop.est[, -1])
mr.out <- data.table::fread("data/validation_beta.tsv",
                            check.names = FALSE)
dat <- as.data.frame(mr.out[, c("chr", "start", "end", inp$sample_id),
                            with = FALSE])
DMRs <- BMI.discovery.sig
data <- merge(DMRs[, 1:3], dat, by = c("chr", "start", "end"), sort = FALSE)
phy.run <- data.frame(base_bmi = inp$base_bmi)
covs <- inp[, c("base_age", "base_gender")]
data.no.na <- na.omit(data)
mr.no.na <- data.matrix(data.no.na[, -(1:3)])
      alpha.all <- pv.all <- NULL
      for(j in 1:ncol(phy.run)){
         alpha.out <- pv.out <- NULL     
          y = phy.run[,j]
          na.ids <- which(is.na(y))
          if(length(na.ids) > 0){
            DummyBeta.m <- mr.no.na[,-na.ids]
            pheno.v <- y[-na.ids]
            frac.m <- prop[-na.ids,]
            covs.m <- covs[-na.ids,]
          }else{
            DummyBeta.m <- mr.no.na
            pheno.v <- y
            frac.m <- prop  
            covs.m <- covs
          }
        for(r in 1:nrow(DummyBeta.m)){
           beta.v <- DummyBeta.m[r,]
           design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
           cov.mod <- model.matrix(~ base_age + factor(base_gender), data=covs.m)[,-1]
           design <- cbind(design, cov.mod)
           fit <- summary(lm(beta.v ~ .-1, data = data.frame(design),na.action=na.omit))
           alpha <- as.numeric(fit$coef[14:26,1])
           alpha.out <- rbind(alpha.out,alpha)
           pv <- as.numeric(fit$coef[14:26,4])
           pv.out <- rbind(pv.out, pv)
        }
        alpha.all <- rbind(alpha.all, data.frame(data.no.na[,1:3],alpha.out, phenotype = colnames(phy.run)[j]))
        pv.all <- rbind(pv.all,data.frame(data.no.na[1:3],pv.out, phenotype = colnames(phy.run)[j]))
      }
   colnames(alpha.all)[c(4:16)] <- paste(cells, ".alpha", sep = "")
   colnames(pv.all)[c(4:16)] <- paste(cells, ".pv", sep = "")
library(plyr)   
out <- join(alpha.all[,c(1:16)], pv.all[,c(1:16)], by = c("chr", "start", "end"))
BMI.validation.all <- out
# ---- 3. Discovery and validation merge ----
BMI.dis <- BMI.discovery.sig
fdr.dis <- BMI.dis[,c(30:42)]
ctps <- as.numeric(which(colSums(fdr.dis<0.05) >1))
alpha.dis <- BMI.dis[,c(4:16)]
pv.dis <- BMI.dis[,c(17:29)] 
alpha.dis.sig <- alpha.dis[, ctps]
pv.dis.sig <- pv.dis[,ctps]
fdr.dis.sig <- fdr.dis[,ctps]
out <- NULL
for(i in 1:ncol(fdr.dis.sig)){
    ids <- as.numeric(which(fdr.dis.sig[,i]<0.05))
    tmp <- data.frame(BMI.dis[ids,c(1:3)], alpha=alpha.dis.sig[ids,i], pv = pv.dis.sig[ids,i], fdr = fdr.dis.sig[ids,i], cell.type = read.table(text = colnames(alpha.dis.sig)[i], sep = ".")[,1])
    out <- rbind(out, tmp)
}
BMI.vd <- BMI.validation.all
alpha.vd <- BMI.vd[,c(4:16)]
pv.vd <- BMI.vd[,c(17:29)] 
alpha.vd.sig <- alpha.vd[, ctps]
pv.vd.sig <- pv.vd[,ctps]
vd.reshape <- NULL
for(i in 1:ncol(alpha.vd.sig)){
    vd.tmp <- data.frame(BMI.vd[,1:3], BMI.alpha.validate = alpha.vd.sig[,i], BMI.pv.validate = pv.vd.sig[,i], cell.type = read.table(text = colnames(alpha.vd.sig)[i], sep = ".")[,1])
    vd.reshape <- rbind(vd.reshape, vd.tmp)
}
out.dis.vd <- merge(out, vd.reshape, by = c("cell.type", "chr", "start", "end"))
write.csv(out.dis.vd, "results/discovery_validation_associations.csv", row.names = F)
out.dis.vd.sub <- out.dis.vd[out.dis.vd$BMI.pv.validate < 0.05 & (out.dis.vd$BMI.alpha.validate * out.dis.vd$alpha>0),]
write.csv(out.dis.vd.sub, "results/discovery_validation_replicated.csv", row.names = F)
write.csv(BMI.discovery.all, "results/discovery_all.csv", row.names = FALSE)
write.csv(BMI.discovery.sig, "results/discovery_significant.csv", row.names = FALSE)
write.csv(BMI.validation.all, "results/validation_all.csv", row.names = FALSE)
