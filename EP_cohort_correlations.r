### Creates data frames for the analysis ###
corr_frame <- data.frame(Subject = as.numeric(TMStable$ID), 
                         P3 = UPDRStable$P3, 
                         P3_UL = UPDRStable$P3UL, 
                         P3_LL = UPDRStable$P3LL, 
                         P3_MAS = UPDRStable$P3_MAS, 
                         P3_LAS = UPDRStable$P3_LAS, 
                         AI = UPDRStable$AI, 
                         mtMA = TMStable$mtMA, 
                         mtLA = TMStable$mtLA, 
                         silMA = TMStable$silMA, 
                         silLA = TMStable$silLA,
                         aucMA = TMStable$aucMA, 
                         aucLA = TMStable$aucLA,
                         mtAI = TMStable$mtAI, 
                         silAI = TMStable$silAI, 
                         silRat = TMStable$silRat, 
                         aucAI = TMStable$aucAI)
corr_frame_ePD <- data.frame(Subject = as.numeric(TMStable_ePD$ID), 
                         P3 = UPDRStable_ePD$P3, 
                         P3_UL = UPDRStable_ePD$P3UL, 
                         P3_LL = UPDRStable_ePD$P3LL, 
                         P3_MAS = UPDRStable_ePD$P3_MAS, 
                         P3_LAS = UPDRStable_ePD$P3_LAS, 
                         AI = UPDRStable_ePD$AI, 
                         mtMA = TMStable_ePD$mtMA, 
                         mtLA = TMStable_ePD$mtLA, 
                         silMA = TMStable_ePD$silMA, 
                         silLA = TMStable_ePD$silLA,
                         aucMA = TMStable_ePD$aucMA, 
                         aucLA = TMStable_ePD$aucLA,
                         mtAI = TMStable_ePD$mtAI, 
                         silAI = TMStable_ePD$silAI, 
                         silRat = TMStable_ePD$silRat, 
                         aucAI = TMStable_ePD$aucAI)
corr_frame_aPD <- data.frame(Subject = as.numeric(TMStable_aPD$ID), 
                         P3 = UPDRStable_aPD$P3, 
                         P3_UL = UPDRStable_aPD$P3UL, 
                         P3_LL = UPDRStable_aPD$P3LL, 
                         P3_MAS = UPDRStable_aPD$P3_MAS, 
                         P3_LAS = UPDRStable_aPD$P3_LAS, 
                         AI = UPDRStable_aPD$AI, 
                         mtMA = TMStable_aPD$mtMA, 
                         mtLA = TMStable_aPD$mtLA, 
                         silMA = TMStable_aPD$silMA, 
                         silLA = TMStable_aPD$silLA,
                         aucMA = TMStable_aPD$aucMA, 
                         aucLA = TMStable_aPD$aucLA,
                         mtAI = TMStable_aPD$mtAI, 
                         silAI = TMStable_aPD$silAI, 
                         silRat = TMStable_aPD$silRat, 
                         aucAI = TMStable_aPD$aucAI)

###################################
### Within-subjects correlation ###
###################################
corr_mask <- c(NA, NA, NA, NA, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA, NA, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
cvars <- c("P3", "P3_UL", "P3_LL", "P3_MAS","P3_LAS","AI","mtMA","mtLA","silMA","silLA","aucMA","aucLA","mtAI","silAI","aucAI")

### All PD subjects ###
my.corr_mat <- rmcorr_mat(participant = Subject, 
                          variables = cvars, 
                          dataset = corr_frame)
my.corr_mat$summary$p.vals <- my.corr_mat$summary$p.vals * corr_mask
my.corr_mat$summary$p.vals.adj <- p.adjust(my.corr_mat$summary$p.vals, method = "fdr")
corr_tab <- my.corr_mat$summary %>% add_significance("p.vals", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
print(subset(corr_tab, p.vals < 0.05))

### Early PD subjects ###
my.corr_mat <- rmcorr_mat(participant = Subject, 
                          variables = cvars, 
                          dataset = corr_frame_ePD)
my.corr_mat$summary$p.vals <- my.corr_mat$summary$p.vals * corr_mask
my.corr_mat$summary$p.vals.adj <- p.adjust(my.corr_mat$summary$p.vals, method = "fdr")
corr_tab <- my.corr_mat$summary %>% add_significance("p.vals", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
print(subset(corr_tab, p.vals < 0.05))

### Advanced PD subjects ###
my.corr_mat <- rmcorr_mat(participant = Subject, 
                          variables = cvars, 
                          dataset = corr_frame_aPD)
my.corr_mat$summary$p.vals <- my.corr_mat$summary$p.vals * corr_mask
my.corr_mat$summary$p.vals.adj <- p.adjust(my.corr_mat$summary$p.vals, method = "fdr")
corr_tab <- my.corr_mat$summary %>% add_significance("p.vals", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
print(subset(corr_tab, p.vals < 0.05))


####################################
### Between-subjects correlation ###
####################################

# Between-subjects correlation function
rmcorr_btw <- function(dataset, participant, measure1, measure2){
  subject_means <- aggregate(cbind(dataset[[measure1]], dataset[[measure2]]) ~ dataset[[participant]], data = dataset, FUN = mean)
  times <- as.numeric(table(as.numeric(dataset[[participant]][(!is.na(dataset[[measure1]]))&(!is.na(dataset[[measure2]]))])))
  corcoef <- cor.test(rep(subject_means$V1, times = times),rep(subject_means$V2, times = times))
  n <- length(times)
  r_value <- corcoef$estimate
  t_value <- (r_value * sqrt(n - 2)) / sqrt(1 - r_value^2)
  p_value <- 2 * (1 - pt(abs(t_value), n - 2))
  return(data.frame(r = r_value, ci.low = corcoef$conf.int[1], ci.up = corcoef$conf.int[2], p = p_value))
}

# Between-subjects correlation table function
rmcorr_btw_tab <- function(dataset, participant){
  clinVars = c("P3", "P3_UL", "P3_LL", "P3_MAS","P3_LAS","AI")
  tmsVars = c("mtMA","mtLA","silMA","silLA","aucMA","aucLA","mtAI","silAI","aucAI")
  nr = length(clinVars)*length(tmsVars)
  bstable <- data.frame(measure1 = numeric(nr),
                        measure2 = numeric(nr),
                        r = numeric(nr),
                        lowerCI = numeric(nr),
                        upperCI = numeric(nr),
                        p.vals = numeric(nr),
                        p.vals.adj = numeric(nr))
  k <- 1
  for(i in clinVars){
    for(j in tmsVars){
      tmp <- rmcorr_btw(dataset,participant,i,j)
      bstable$measure1[k] <- i
      bstable$measure2[k] <- j
      bstable$r[k] <- tmp$r
      bstable$lowerCI[k] <- tmp$ci.low
      bstable$upperCI[k] <- tmp$ci.up
      bstable$p.vals[k] <- tmp$p
      k <- k + 1
    }
  }
  bstable$p.vals.adj <- p.adjust(bstable$p.vals, method = "fdr")
  bstable <- add_significance(bstable, "p.vals", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  return(bstable)
}

### All PD subjects ###
corr_tab <- rmcorr_btw_tab(participant = "Subject", dataset = corr_frame)
print(subset(corr_tab, p.vals < 0.05))

### Early PD subjects ###
corr_tab <- rmcorr_btw_tab(participant = "Subject", dataset = corr_frame_ePD)
print(subset(corr_tab, p.vals < 0.05))

### Advanced PD subjects ###
corr_tab <- rmcorr_btw_tab(participant = "Subject", dataset = corr_frame_aPD)
print(subset(corr_tab, p.vals < 0.05))