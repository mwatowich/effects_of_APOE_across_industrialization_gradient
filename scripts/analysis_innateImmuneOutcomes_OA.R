### Linear models of innate immune outcomes in Orang Asli
# Analyses performed similarly for Turkana 

# APOE_linear = continuous variable of APOE genotypes, described in Watowich et al. 

# First, prep data: perform log2() transformations on cell counts, CRP, cytokines

# Specify immune outcomes of interest 
outcomes_immune<-c("neutrophils_perc","lymphocytes_perc","monocytes_perc",
                   "neutrophils2","lymphocytes2","monocytes2",
                   "wbc","NtoL","crp","tnfa","ifn","il10","il1B","il6","il12", 
                   "il2","il4","il8")

###########################################################
# Effects of genotype and urbanicity, genotype*urbanicity #
###########################################################

# Models 
out_immune_oa <- do.call(rbind,lapply(outcomes_immune, function(outcome) {
  df <- na.omit(apoe_oa[, c("age","sex","APOE_linear","urb_score",outcome)])
  
  # Run linear additive model
  model <- lm(scale(df[[outcome]]) ~ scale(age) + sex + scale(urb_score) + scale(APOE_linear), df)
  tidy_res <- tidy(model)[-1,c(1:3,5)]
  # clean up 
  tidy_res$term <- gsub("urb_score","UrbScore",
                        gsub("_linear","",
                             gsub("scale(","",
                                  gsub(")","",tidy_res$term, fixed = T), fixed = T), fixed = T))
  colnames(tidy_res) <- c("covariate","beta","se","pval")
  tidy_res_wide <- tidy_res %>%
    tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval), names_glue = "{.value}_{covariate}") %>% 
    mutate(n = nrow(df), 
           outcome = outcome, 
           type = "additive")
  
  # Run linear interactive model
  model_int <- lm(scale(df[[outcome]]) ~ scale(age) + sex + scale(APOE_linear)*scale(urb_score), df)
  tidy_res_i <- tidy(model_int)[-1,c(1:3,5)]
  # clean up 
  tidy_res_i$term <- gsub(":","x",gsub("urb_score","UrbScore",
                                       gsub("_linear","",
                                            gsub("scale(","",
                                                 gsub(")","",tidy_res_i$term, 
                                                      fixed = T), fixed = T), fixed = T), fixed = T), fixed = T)
  colnames(tidy_res_i) <- c("covariate","beta","se","pval")
  tidy_res_i_wide <- tidy_res_i %>%
    tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval), names_glue = "{.value}_{covariate}") %>% 
    mutate(n = nrow(df), 
           outcome = outcome, 
           type = "interactive")
  
  # Combine 
  out <- dplyr::bind_rows(tidy_res_wide, tidy_res_i_wide)
  out <- out %>% 
    dplyr::select(c("outcome","type","n",everything()))
  return(out)
}))

# Adjust p-values for multiple columns
out_immune_oa <- out_immune_oa %>% 
  group_by(type) %>% 
  mutate(across(starts_with("pval_"), ~ p.adjust(.x, method = "BH"), .names = "padj_{.col}")) %>% 
  as.data.frame()
colnames(out_immune_oa) <- gsub("padj_pval","padj",colnames(out_immune_oa))



###########################################################
#               Effects of urbanicity*lipids              #
###########################################################

# Model 
immuneLipid_out_OA <- do.call(rbind, lapply(outcomes_immune, function(immune) {
  do.call(rbind,lapply(c("hdl","ldl","total_chol"), function(lipid) {
    
    # Subset 
    df <- na.omit(apoe_oa[,c("age","sex","urb_score",lipid,immune)])
    df$lipid <- df[[lipid]]
    
    # Run linear additive model
    model <- lm(scale(df[[immune]]) ~ scale(age) + sex + scale(urb_score)*scale(lipid), df)
    tidy_res <- tidy(model)[-1,c(1:3,5)]
    # clean up 
    tidy_res$term <- gsub(":","x",gsub("urb_score","UrbScore",
                                       gsub("scale(","",
                                            gsub(")","",tidy_res$term, 
                                                 fixed = T), fixed = T), fixed = T), fixed = T)
    colnames(tidy_res) <- c("covariate","beta","se","pval")
    tidy_res_wide <- tidy_res %>%
      tidyr::pivot_wider(names_from = covariate, values_from = c(beta, se, pval), names_glue = "{.value}_{covariate}") %>% 
      mutate(immune_outcome = immune,
             n = nrow(df), 
             lipid = lipid)
    
    # Combine 
    out <- tidy_res_wide %>% 
      dplyr::select(c("immune_outcome","lipid","n",everything()))
    return(out)
  }))
})) %>% as.data.frame()

# Adjust p-values for multiple columns
immuneLipid_out_OA <- immuneLipid_out_OA %>% 
  group_by(immune_outcome, lipid) %>% 
  mutate(across(starts_with("pval_"), ~ p.adjust(.x, method = "BH"), .names = "padj_{.col}")) %>% 
  as.data.frame()
colnames(immuneLipid_out_OA) <- gsub("padj_pval","padj",colnames(immuneLipid_out_OA))
