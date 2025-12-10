### Linear modeling of cholesterol traits. Additive effects of genotype and urbanicity, interactive effects of genotype*urbanicity  

# APOE_linear = continuous variable of APOE genotypes, described in Watowich et al. 

# Combine data for Orang Asli and Turkana 
lipids_df <- rbind(data.frame(
  apoe_turkana[,c("hdl","ldl","total_cholesterol",
                  "APOE_linear",
                  "urb_score","age","sex")], 
  population = "Turkana"),
  data.frame(apoe_oa[,c("hdl","ldl","total_chol",
                        "APOE_linear",
                        "urb_score","age","sex")], 
             population = "Orang Asli") %>% 
    dplyr::rename(total_cholesterol = total_chol)) %>% 
  dplyr::rename(LDL = ldl, 
                HDL = hdl,
                Total_cholesterol = total_cholesterol)

# Run main model of APOE_linear and interactive model
out_chol_both <- as.data.frame(do.call(rbind, lapply(c("Orang Asli","Turkana"), function(pop) {
  do.call(rbind,lapply(c("HDL","LDL","Total_cholesterol"), function(outcome) {
    
    # Subset to complete cases 
    df <- na.omit(lipids_df[which(lipids_df$population == pop),
                            c("age","sex","APOE_linear","urb_score",outcome)]) 
    
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
      mutate(population = pop,
             n = nrow(df), 
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
      mutate(population = pop,
             n = nrow(df), 
             outcome = outcome, 
             type = "interactive")
    
    # Combine 
    out <- dplyr::bind_rows(tidy_res_wide, tidy_res_i_wide)
    out <- out %>% 
      dplyr::select(c("outcome","type","n",everything()))
    return(out)
  }))
})))

# Adjust p-values for multiple columns
out_chol_both <- out_chol_both %>% 
  group_by(type, population) %>% 
  mutate(across(starts_with("pval_"), ~ p.adjust(.x, method = "BH"), .names = "padj_{.col}")) %>% 
  as.data.frame()
colnames(out_chol_both) <- gsub("padj_pval","padj",colnames(out_chol_both))
