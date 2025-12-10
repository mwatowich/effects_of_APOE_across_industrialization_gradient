### Differential gene expression analysis of APOE genotypes in Orang Asli, on baseline PBMC expression 

library(dplyr)
library(limma)
library(ggplot2)
library(EMMREML)


### Files 
# Normalized RNAseq counts = exp_pbmc
# K matrix = K_pbmc 
# Metadata = meta_pbmc


### PCA of cell composition 
cell_composition <- meta_pbmc[, c("lymphocytes_perc", "eosinophils_perc", "neutrophils_perc", "monocytes_perc")]
pca <- prcomp(cell_scaled, center = TRUE, scale. = TRUE)
summary(pca)
meta_pbmc <- cbind(meta_pbmc, pca$x[, 1:3]) 


### Model GE ~ APOE_linear ---------

# Design matrix 
X <- model.matrix(~ age + APOE_linear + sex + urb_score + PC1 + PC2 + PC3, data = meta_pbmc)

# Z matrix
Z_pbmc <- diag(nrow(K_pbmc))
rownames(Z_pbmc) <- rownames(K_pbmc)
colnames(Z_pbmc) <- colnames(K_pbmc)

# Run EMMREML for each gene
out_pbmc_apoe <- do.call(rbind, lapply(1:nrow(exp_pbmc), function(i) { 
  print(paste("Starting row:", i))
  y <- as.numeric(exp_pbmc[i, ])
  
  emm <- emmreml(y = y, X = X, Z = Z_pbmc, K = K_pbmc,
                 varbetahat = TRUE, varuhat = TRUE, PEVuhat = TRUE, test = TRUE)
  
  if (!is.null(emm) && !is.null(emm$betahat)) {
    data.frame(
      gene = rownames(exp_pbmc)[i],
      term = colnames(X),
      beta = emm$betahat,
      se = sqrt(emm$varbetahat),
      pval = emm$pvalbeta[, 8]
    )
  } else {
    data.frame(
      gene = rownames(exp_pbmc)[i],
      term = colnames(X),
      beta = NA,
      se = NA,
      pval = NA
    )
  }
}))

# Pivot wider
out_pbmc_apoe <- out_pbmc_apoe |>
  filter(!term=="(Intercept)") %>%
  tidyr::pivot_wider(names_from = term, values_from = c(beta, se, pval)) %>%
  as.data.frame()


### GSEA enrichment analysis ---------------
library(fgsea)

# Load GSEA hallmarks 
hallmarks <- gmtPathways("h.all.v7.5.1.symbols.gmt.txt")
apoe_sig <- out_pbmc_apoe %>% dplyr::select(gene, beta_APOE_linear) %>% arrange(desc(beta_APOE_linear))
apoe_sig_vec <- apoe_sig$beta_APOE_linear
names(apoe_sig_vec) = apoe_sig$gene

fgsea.apoe <- fgsea(pathways = hallmarks, 
                    stats = apoe_sig_vec,
                    minSize = 15,
                    maxSize = 500,
                    eps = 0.0)

fgsea.apoe %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  head(20)

