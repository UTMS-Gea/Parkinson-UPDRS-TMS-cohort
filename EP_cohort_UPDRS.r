if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if (!requireNamespace("rstatix", quietly = TRUE)) install.packages("rstatix")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggprism", quietly = TRUE)) install.packages("ggprism")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("robustlmm", quietly = TRUE)) install.packages("robustlmm")
if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")

library(ggpubr)
library(rstatix)
library(tidyr)
library(ggprism)
library(cowplot)
library(car)
library(lme4)
library(robustlmm)
library(emmeans)
library(lmerTest)

### Load table ###
fpath <- "~/R/data_UDPRS_TMS/"
UPDRStable <- read.table(file = paste(fpath,"UPDRS.csv",sep = ""), sep = ",", header = TRUE)
UPDRStable$Period <- factor(UPDRStable$Period)
UPDRStable$ID <- factor(UPDRStable$ID , levels=as.character(1:22))
color1 <- "#00AFBB"
color2 <- "#5594c9"
color3 <- "#f8766d"

### Create 'Stage' column ###
UPDRStable$Stage <- ifelse(UPDRStable$HY<2.5, "Early", "Advanced")

### Create MAS & LAS columns ###
for(i in 1:nrow(UPDRStable)){
  if(UPDRStable$MAS[i]=="R"){
    UPDRStable$P3_MAS[i] <- UPDRStable$P3right[i]
    UPDRStable$P3_LAS[i] <- UPDRStable$P3left[i]
  }else{
    UPDRStable$P3_MAS[i] <- UPDRStable$P3left[i]
    UPDRStable$P3_LAS[i] <- UPDRStable$P3right[i]
  }
}

### Asymetry index ###
UPDRStable$AI <- (UPDRStable$P3_MAS-UPDRStable$P3_LAS)/(UPDRStable$P3_MAS+UPDRStable$P3_LAS)

### Levene Test ###
leveneTest(P1 ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(P2 ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(P3 ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(P3UL ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(P3LL ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(P3_MAS ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(P3_LAS ~ Period, UPDRStable)$`Pr(>F)`[1]
leveneTest(AI ~ Period, UPDRStable)$`Pr(>F)`[1]

### Plot Function ###
plot_UPDRS <- function(dataset, metric, model, strTitle, pcolor, ylims){
  n <- nrow(dataset)/4
  stat.test <- emmeans(model, pairwise ~ Period)$contrasts %>% summary() %>% add_significance("p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  stat.test$group1 = c('T1','T1','T1','T2','T2','T4')
  stat.test$group2 = c('T2','T4','T5','T4','T5','T5')
  tmp <- stat.test$p.value < 0.05
  tmp[tmp==FALSE] <- NA
  tmp <- tmp + ifelse(!is.na(tmp), cumsum(!is.na(tmp)) * 0.16, 0)
  stat.test$y.position = tmp * max(dataset[[deparse(substitute(metric))]], na.rm = TRUE)
  stat.test$xmin = c(1,1,1,2,2,4)
  stat.test$xmax = c(2,4,5,4,5,5)
  if(missing(ylims)){
    ylims <- c(ceiling(max(stat.test$y.position, na.rm = TRUE)),10,5)
  }
  summary_stats <- dataset %>%
    group_by(PeriodNum) %>%
    dplyr::summarise(
      mean = mean({{metric}}, na.rm = TRUE),
      ci_lower = mean({{metric}}, na.rm = TRUE) - 1.96 * (sd({{metric}}, na.rm = TRUE) / sqrt(n)),
      ci_upper = mean({{metric}}, na.rm = TRUE) + 1.96 * (sd({{metric}}, na.rm = TRUE) / sqrt(n))
    )
  p1 <- ggplot(dataset,aes(x = PeriodNum, y = {{metric}}, group = ID)) +
    geom_line(linewidth = 0.5, colour = pcolor, alpha = 0.3) +
    geom_point(size = 4, colour = pcolor, alpha = 0.2) +
    scale_y_continuous(breaks = seq(0,ylims[1],ylims[2]), minor_breaks = seq(0,ylims[1],ylims[3]), guide = guide_prism_minor()) +
    coord_cartesian(ylim = c(0, ylims[1])) +
    scale_x_continuous(breaks = c(1, 2, 4, 5), labels = c("2018", "2019", "2021", "2022")) +
    geom_vline(xintercept = 3, linetype = "dotted", linewidth = 1, alpha = 0.3) +
    xlab("") +
    ylab("Score") +
    stat_pvalue_manual(stat.test,
                       label = "p.value.signif",
                       label.size = 8,
                       tip.length = 0.03,
                       hide.ns = TRUE) +
    geom_point(inherit.aes = FALSE, data = summary_stats, aes(x = PeriodNum, y = mean), size = 3, color = "black") +
    geom_line(data = summary_stats, aes(x = PeriodNum, y = mean, group = 1), linewidth = 0.2, color = "black") +
    geom_errorbar(inherit.aes = FALSE, data = summary_stats, aes(x = PeriodNum, ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "black") +
    theme_prism() +
    theme(axis.text.x = element_text(angle = 45, size = 11))
  title <- ggdraw() + draw_label(strTitle, fontface = 'bold', size = 16)
  p4 <- plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1))
  return(p4)
}

### p-value aproximation function ###
rlmer_aprox_pval <- function(rlmer_model, predictor) {
  robust_coef <- fixef(rlmer_model)
  robust_vcov <- vcov(rlmer_model)
  model_formula <- formula(rlmer_model)
  model_data <- model.frame(rlmer_model)
  lmer_model <- lmer(model_formula, data = model_data)
  X <- model.matrix(terms(rlmer_model), data = model_data)
  predictor_cols <- grep(predictor, colnames(X))
  X_predictor <- X[, predictor_cols, drop = FALSE]
  beta_predictor <- robust_coef[predictor_cols]
  vcov_predictor <- robust_vcov[predictor_cols, predictor_cols]
  Wald_stat <- as.numeric(t(beta_predictor) %*% solve(vcov_predictor) %*% beta_predictor)
  df1 <- length(beta_predictor)
  anova_lmer <- anova(lmer_model, type = "III")
  df2 <- anova_lmer[predictor, "DenDF"]
  F_stat <- Wald_stat / df1
  p_value <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  partial_eta_squared <- (F_stat * df1) / (F_stat * df1 + df2)  
  return(data.frame(F_stat = F_stat, df1 = df1, df2 = df2, p = p_value, p_eta_sqrt = partial_eta_squared))
}

#######################
### All PD subjects ###
#######################

### Models ###
mod_P1 <- rlmer(P1 ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_P2 <- rlmer(P2 ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3 <- rlmer(P3 ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3UL <- rlmer(P3UL ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3LL <- rlmer(P3LL ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3MAS <- rlmer(P3_MAS ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3LAS <- rlmer(P3_LAS ~ Period + (1|ID), data = UPDRStable, REML = TRUE)
mod_AI <- rlmer(AI ~ Period + (1|ID), data = UPDRStable, REML = TRUE)

round(rlmer_aprox_pval(mod_P1 ,"Period"),2)
round(rlmer_aprox_pval(mod_P2 ,"Period"),2)
round(rlmer_aprox_pval(mod_P3 ,"Period"),2)
round(rlmer_aprox_pval(mod_P3UL ,"Period"),2)
round(rlmer_aprox_pval(mod_P3LL ,"Period"),2)
round(rlmer_aprox_pval(mod_P3MAS ,"Period"),2)
round(rlmer_aprox_pval(mod_P3LAS ,"Period"),2)
round(rlmer_aprox_pval(mod_AI ,"Period"),2)

### Plots ###

pP1 <- plot_UPDRS(dataset = UPDRStable, metric = P1, model = mod_P1, strTitle = "", color1, c(50,10,5))# c(40,10,5)
pP2 <-  plot_UPDRS(dataset = UPDRStable, metric = P2, model = mod_P2, strTitle = "", color1, c(80,10,5))# c(60,10,5)
pP3 <- plot_UPDRS(dataset = UPDRStable, metric = P3, model = mod_P3, strTitle = "", color1, c(140,20,10))# c(130,20,10)
pP3UL <- plot_UPDRS(dataset = UPDRStable, metric = P3UL, model = mod_P3UL, strTitle = "", color1, c(70,10,5))# c(70,10,5)
pP3LL <- plot_UPDRS(dataset = UPDRStable, metric = P3LL, model = mod_P3LL, strTitle = "", color1, c(50,10,5))# c(50,10,5)
pP3_MAS <- plot_UPDRS(dataset = UPDRStable, metric = P3_MAS, model = mod_P3MAS, strTitle = "", color1, c(50,10,5))# c(50,10,5)
pP3_LAS <- plot_UPDRS(dataset = UPDRStable, metric = P3_LAS, model = mod_P3LAS, strTitle = "", color1, c(50,10,5))# c(50,10,5)
pAI <- plot_UPDRS(dataset = UPDRStable, metric = AI, model = mod_AI, strTitle = "", color1, c(1,0.2,0.1))


#########################
### Early PD subgroup ###
#########################

### Filter early PD subjects (H&Y 1-2)
subSample <- c(1,2,3,5,6,9,10,12,13,14,15,16,17,18,20)
UPDRStable_ePD <- UPDRStable[UPDRStable$ID %in% subSample,]

### Models ###
mod_P1e <- rlmer(P1 ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_P2e <- rlmer(P2 ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_P3e <- rlmer(P3 ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_P3ULe <- rlmer(P3UL ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_P3LLe <- rlmer(P3LL ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_P3MASe <- rlmer(P3_MAS ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_P3LASe <- rlmer(P3_LAS ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)
mod_AIe <- rlmer(AI ~ Period + (1|ID), data = UPDRStable_ePD, REML = TRUE)

round(rlmer_aprox_pval(mod_P1e ,"Period"),2)
round(rlmer_aprox_pval(mod_P2e ,"Period"),2)
round(rlmer_aprox_pval(mod_P3e ,"Period"),2)
round(rlmer_aprox_pval(mod_P3ULe ,"Period"),2)
round(rlmer_aprox_pval(mod_P3LLe ,"Period"),2)
round(rlmer_aprox_pval(mod_P3MASe ,"Period"),2)
round(rlmer_aprox_pval(mod_P3LASe ,"Period"),2)
round(rlmer_aprox_pval(mod_AIe ,"Period"),2)

### Plots ###

pP1e <- plot_UPDRS(dataset = UPDRStable_ePD, metric = P1, model = mod_P1e, strTitle = "", color2, c(50,10,5))# c(40,10,5)
pP2e <-  plot_UPDRS(dataset = UPDRStable_ePD, metric = P2, model = mod_P2e, strTitle = "", color2, c(80,10,5))# c(60,10,5)
pP3e <- plot_UPDRS(dataset = UPDRStable_ePD, metric = P3, model = mod_P3e, strTitle = "", color2, c(140,20,10))# c(130,20,10)
pP3ULe <- plot_UPDRS(dataset = UPDRStable_ePD, metric = P3UL, model = mod_P3ULe, strTitle = "", color2, c(70,10,5))# c(70,10,5)
pP3LLe <- plot_UPDRS(dataset = UPDRStable_ePD, metric = P3LL, model = mod_P3LLe, strTitle = "", color2, c(50,10,5))# c(50,10,5)
pP3_MASe <- plot_UPDRS(dataset = UPDRStable_ePD, metric = P3_MAS, model = mod_P3MASe, strTitle = "", color2, c(50,10,5))# c(50,10,5)
pP3_LASe <- plot_UPDRS(dataset = UPDRStable_ePD, metric = P3_LAS, model = mod_P3LASe, strTitle = "", color2, c(50,10,5))# c(50,10,5)
pAIe <- plot_UPDRS(dataset = UPDRStable_ePD, metric = AI, model = mod_AIe, strTitle = "", color2, c(1,0.2,0.1))


############################
### Advanced PD subgroup ###
############################

### Filter advanced PD subjects (H&Y 2.5-4)
subSample <- c(4,7,8,11,19,21,22)
UPDRStable_aPD <- UPDRStable[UPDRStable$ID %in% subSample,]

### Models ###
mod_P1a <- rlmer(P1 ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_P2a <- rlmer(P2 ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_P3a <- rlmer(P3 ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_P3ULa <- rlmer(P3UL ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_P3LLa <- rlmer(P3LL ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_P3MASa <- rlmer(P3_MAS ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_P3LASa <- rlmer(P3_LAS ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)
mod_AIa <- rlmer(AI ~ Period + (1|ID), data = UPDRStable_aPD, REML = TRUE)

round(rlmer_aprox_pval(mod_P1a ,"Period"),2)
round(rlmer_aprox_pval(mod_P2a ,"Period"),2)
round(rlmer_aprox_pval(mod_P3a ,"Period"),2)
round(rlmer_aprox_pval(mod_P3ULa ,"Period"),2)
round(rlmer_aprox_pval(mod_P3LLa ,"Period"),2)
round(rlmer_aprox_pval(mod_P3MASa ,"Period"),2)
round(rlmer_aprox_pval(mod_P3LASa ,"Period"),2)
round(rlmer_aprox_pval(mod_AIa ,"Period"),2)

### Plots ###

pP1a <- plot_UPDRS(dataset = UPDRStable_aPD, metric = P1, model = mod_P1a, strTitle = "", color3, c(50,10,5))# c(40,10,5)
pP2a <-  plot_UPDRS(dataset = UPDRStable_aPD, metric = P2, model = mod_P2a, strTitle = "", color3, c(80,10,5))# c(60,10,5)
pP3a <- plot_UPDRS(dataset = UPDRStable_aPD, metric = P3, model = mod_P3a, strTitle = "", color3, c(140,20,10))# c(130,20,10)
pP3ULa <- plot_UPDRS(dataset = UPDRStable_aPD, metric = P3UL, model = mod_P3ULa, strTitle = "", color3, c(70,10,5))# c(70,10,5)
pP3LLa <- plot_UPDRS(dataset = UPDRStable_aPD, metric = P3LL, model = mod_P3LLa, strTitle = "", color3, c(50,10,5))# c(50,10,5)
pP3_MASa <- plot_UPDRS(dataset = UPDRStable_aPD, metric = P3_MAS, model = mod_P3MASa, strTitle = "", color3, c(50,10,5))# c(50,10,5)
pP3_LASa <- plot_UPDRS(dataset = UPDRStable_aPD, metric = P3_LAS, model = mod_P3LASa, strTitle = "", color3, c(50,10,5))# c(50,10,5)
pAIa <- plot_UPDRS(dataset = UPDRStable_aPD, metric = AI, model = mod_AIa, strTitle = "", color3, c(1,0.2,0.1))

# Fig. 2: 1500 x 1200
t1 <- ggdraw() + draw_label("UPDRS - Part I", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("UPDRS - Part II", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t3 <- ggdraw() + draw_label("UPDRS - Part III", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("UPDRS - Part III\n(asymmetry index)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color1)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color2)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color3)
plot_grid(NULL, t1, t2, t3, t4, 
          L1, pP1, pP2, pP3, pAI, 
          L2, pP1e, pP2e, pP3e, pAIe, 
          L3, pP1a, pP2a, pP3a, pAIa, 
          ncol = 5, rel_heights = c(0.15,1,1,1), rel_widths = c(0.1,1,1,1,1))

# Sup. Fig. 1: 1500 x 1200
t1 <- ggdraw() + draw_label("UPDRS - Part III\n(upper limbs)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("UPDRS - Part III\n(lower limbs)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t3 <- ggdraw() + draw_label("UPDRS - Part III\n(most affected side)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("UPDRS - Part III\n(least affected side)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color1)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color2)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color3)
plot_grid(NULL, t1, t2, t3, t4, 
          L1, pP3UL, pP3LL, pP3_MAS, pP3_LAS, 
          L2, pP3ULe, pP3LLe, pP3_MASe, pP3_LASe, 
          L3, pP3ULa, pP3LLa, pP3_MASa, pP3_LASa, 
          ncol = 5, rel_heights = c(0.15,1,1,1), rel_widths = c(0.1,1,1,1,1))

### Analysis by Sex ###
mod_P1_sex <- rlmer(P1 ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_P2_sex <- rlmer(P2 ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3_sex <- rlmer(P3 ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3UL_sex <- rlmer(P3UL ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3LL_sex <- rlmer(P3LL ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3MAS_sex <- rlmer(P3_MAS ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3LAS_sex <- rlmer(P3_LAS ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)
mod_AI_sex <- rlmer(AI ~ Period + Sex + (1|ID), data = UPDRStable, REML = TRUE)

rlmer_aprox_pval(mod_P1_sex ,"Sex")
rlmer_aprox_pval(mod_P2_sex ,"Sex")
rlmer_aprox_pval(mod_P3_sex ,"Sex")
rlmer_aprox_pval(mod_P3UL_sex ,"Sex")
rlmer_aprox_pval(mod_P3LL_sex ,"Sex")
rlmer_aprox_pval(mod_P3MAS_sex ,"Sex")
rlmer_aprox_pval(mod_P3LAS_sex ,"Sex")
rlmer_aprox_pval(mod_AI_sex ,"Sex")


### Analysis by stage ###
mod_P1_stage <- rlmer(P1 ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_P2_stage <- rlmer(P2 ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3_stage <- rlmer(P3 ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3UL_stage <- rlmer(P3UL ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3LL_stage <- rlmer(P3LL ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3MAS_stage <- rlmer(P3_MAS ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_P3LAS_stage <- rlmer(P3_LAS ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)
mod_AI_stage <- rlmer(AI ~ Period + Stage + (1|ID), data = UPDRStable, REML = TRUE)

rlmer_aprox_pval(mod_P1_stage ,"Stage")
rlmer_aprox_pval(mod_P2_stage ,"Stage")
rlmer_aprox_pval(mod_P3_stage ,"Stage")
rlmer_aprox_pval(mod_P3UL_stage ,"Stage")
rlmer_aprox_pval(mod_P3LL_stage ,"Stage")
rlmer_aprox_pval(mod_P3MAS_stage ,"Stage")
rlmer_aprox_pval(mod_P3LAS_stage ,"Stage")
rlmer_aprox_pval(mod_AI_stage ,"Stage")
