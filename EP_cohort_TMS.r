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

fpath <- "~/R/data_UDPRS_TMS/"
TMStable <- read.table(file = paste(fpath,"TMS.csv",sep = ""), sep = ",", header = TRUE)
TMStable$Period <- factor(TMStable$Period)
TMStable$ID <- factor(TMStable$ID , levels=as.character(1:22))
color1 <- "#00AFBB"
color2 <- "#5594c9"
color3 <- "#f8766d"

### Create 'Stage' column ###
TMStable$Stage <- ifelse(TMStable$HY<2.5, "Early", "Advanced")

## S50 Amplitude (130%UMR) ###
TMStable$ampMA <- RCtable[RCtable$Intensity==130,]$ampMA
TMStable$ampLA <- RCtable[RCtable$Intensity==130,]$ampLA

### Asymetry indices ###
TMStable$mtAI <- (TMStable$mtMA-TMStable$mtLA)/(TMStable$mtMA+TMStable$mtLA)
TMStable$aucAI <- (TMStable$aucMA-TMStable$aucLA)/(TMStable$aucMA+TMStable$aucLA)
TMStable$silAI <- (TMStable$silMA-TMStable$silLA)/(TMStable$silMA+TMStable$silLA)
TMStable$ampAI <- (TMStable$ampMA-TMStable$ampLA)/(TMStable$ampMA+TMStable$ampLA)

TMStable$silRat <- TMStable$silLA/TMStable$silMA
tmpDF <- data.frame(T1 = TMStable$silRat[TMStable$Period == "T1"],
                    T2 = TMStable$silRat[TMStable$Period == "T2"],
                    T4 = TMStable$silRat[TMStable$Period == "T4"],
                    T5 = TMStable$silRat[TMStable$Period == "T5"])
silRatAve <- rep(rowMeans(tmpDF, na.rm = TRUE),4)
TMStable$silRat <- (TMStable$silRat-silRatAve)/silRatAve

### Levene Test ###
leveneTest(mtMA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(mtLA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(latMA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(latLA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(durMA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(durLA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(silMA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(silLA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(aucMA ~ Period, TMStable)$`Pr(>F)`[1]
leveneTest(aucLA ~ Period, TMStable)$`Pr(>F)`[1]

### Plot Function ###
plot_TMS <- function(dataset, metric, model, strTitle, yLabel, pcolor, ylims){
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
    ylims <- c(0,ceiling(max(stat.test$y.position, na.rm = TRUE)),10,5)
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
    scale_y_continuous(breaks = seq(ylims[1],ylims[2],ylims[3]), minor_breaks = seq(ylims[1],ylims[2],ylims[4]), guide = guide_prism_minor()) +
    coord_cartesian(ylim = c(ylims[1], ylims[2])) +
    scale_x_continuous(breaks = c(1, 2, 4, 5), labels = c("2018", "2019", "2021", "2022")) +
    geom_vline(xintercept = 3, linetype = "dotted", linewidth = 1, alpha = 0.3) +
    xlab("") +
    ylab(yLabel) +
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
  p4 <- plot_grid(title, p1, ncol = 1, rel_heights = c(0.15, 1))
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
mod_mtMA <- rlmer(mtMA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_mtLA <- rlmer(mtLA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_latMA <- rlmer(latMA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_latLA <- rlmer(latLA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_durMA <- rlmer(durMA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_durLA <- rlmer(durLA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_silMA <- rlmer(silMA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_silLA <- rlmer(silLA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_ampMA <- rlmer(ampMA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_ampLA <- rlmer(ampLA ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_aucMA <- rlmer(aucMA ~ Period + mtMA + (1|ID), data = TMStable, REML = TRUE)
mod_aucLA <- rlmer(aucLA ~ Period + mtLA + (1|ID), data = TMStable, REML = TRUE)
mod_mtAI <- rlmer(mtAI ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_aucAI <- rlmer(aucAI ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_silAI <- rlmer(silAI ~ Period + (1|ID), data = TMStable, REML = TRUE)
mod_ampAI <- rlmer(ampAI ~ Period + (1|ID), data = TMStable, REML = TRUE)

round(rlmer_aprox_pval(mod_mtMA,"Period"),3)
round(rlmer_aprox_pval(mod_mtLA,"Period"),3)
round(rlmer_aprox_pval(mod_latMA,"Period"),3)
round(rlmer_aprox_pval(mod_latLA,"Period"),3)
round(rlmer_aprox_pval(mod_durMA,"Period"),3)
round(rlmer_aprox_pval(mod_durLA,"Period"),3)
round(rlmer_aprox_pval(mod_silMA,"Period"),3)
round(rlmer_aprox_pval(mod_silLA,"Period"),3)
round(rlmer_aprox_pval(mod_ampMA,"Period"),3)
round(rlmer_aprox_pval(mod_ampLA,"Period"),3)
round(rlmer_aprox_pval(mod_aucMA,"Period"),3)
round(rlmer_aprox_pval(mod_aucLA,"Period"),3)
round(rlmer_aprox_pval(mod_mtAI,"Period"),3)
round(rlmer_aprox_pval(mod_aucAI,"Period"),3)
round(rlmer_aprox_pval(mod_silAI,"Period"),3)
round(rlmer_aprox_pval(mod_ampAI,"Period"),3)

### Plots ###
pmtMA <- plot_TMS(dataset = TMStable, metric = mtMA, model = mod_mtMA, strTitle = "", yLabel = "%MSA", color1, c(0,150,20,10))# c(0,140,20,10)
pmtLA <- plot_TMS(dataset = TMStable, metric = mtLA, model = mod_mtLA, strTitle = "", yLabel = "%MSA", color1, c(0,150,20,10))# c(0,140,20,10)
platMA <- plot_TMS(dataset = TMStable, metric = latMA, model = mod_latMA, strTitle = "", yLabel = "Time (ms)", color1, c(0,30,10,5))
platLA <- plot_TMS(dataset = TMStable, metric = latLA, model = mod_latLA, strTitle = "", yLabel = "Time (ms)", color1, c(0,30,10,5))
pdurMA <- plot_TMS(dataset = TMStable, metric = durMA, model = mod_durMA, strTitle = "", yLabel = "Time (ms)", color1, c(0,20,10,5))
pdurLA <- plot_TMS(dataset = TMStable, metric = durLA, model = mod_durLA, strTitle = "", yLabel = "Time (ms)", color1, c(0,20,10,5))
psilMA <- plot_TMS(dataset = TMStable, metric = silMA, model = mod_silMA, strTitle = "", yLabel = "Time (ms)", color1, c(0,180,20,10))# c(0,140,20,10)
psilLA <- plot_TMS(dataset = TMStable, metric = silLA, model = mod_silLA, strTitle = "", yLabel = "Time (ms)", color1, c(0,180,20,10))# c(0,140,20,10)
pampMA <- plot_TMS(dataset = TMStable, metric = ampMA, model = mod_ampMA, strTitle = "", yLabel = "Time (ms)", color1, c(0,12,3,1))
pampLA <- plot_TMS(dataset = TMStable, metric = ampLA, model = mod_ampLA, strTitle = "", yLabel = "Time (ms)", color1, c(0,12,3,1))
paucMA <- plot_TMS(dataset = TMStable, metric = aucMA, model = mod_aucMA, strTitle = "", yLabel = "AUC (mV∙ms)", color1, c(0,70,20,10))# c(0,70,20,10))
paucLA <- plot_TMS(dataset = TMStable, metric = aucLA, model = mod_aucLA, strTitle = "", yLabel = "AUC (mV∙ms)", color1, c(0,70,20,10))# c(0,70,20,10))
pmtAI <- plot_TMS(dataset = TMStable, metric = mtAI, model = mod_mtAI, strTitle = "", yLabel = "Index", color1, c(-1,1,0.5,0.25))
paucAI <- plot_TMS(dataset = TMStable, metric = aucAI, model = mod_aucAI, strTitle = "", yLabel = "Index", color1, c(-1,1,0.5,0.25))
psilAI <- plot_TMS(dataset = TMStable, metric = silAI, model = mod_silAI, strTitle = "", yLabel = "Index", color1, c(-1,1,0.5,0.25))
pampAI <- plot_TMS(dataset = TMStable, metric = ampAI, model = mod_ampAI, strTitle = "", yLabel = "Index", color1, c(-1,1,0.5,0.25))


#########################
### Early PD subgroup ###
#########################

### Filter early PD subjects (H&Y 1-2)
subSample <- c(1,2,3,5,6,9,10,12,13,14,15,16,17,18,20)
TMStable_ePD <- TMStable[TMStable$ID %in% subSample,]

### Models ###
mod_mtMAe <- rlmer(mtMA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_mtLAe <- rlmer(mtLA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_latMAe <- rlmer(latMA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_latLAe <- rlmer(latLA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_durMAe <- rlmer(durMA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_durLAe <- rlmer(durLA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_silMAe <- rlmer(silMA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_silLAe <- rlmer(silLA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_ampMAe <- rlmer(ampMA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_ampLAe <- rlmer(ampLA ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_aucMAe <- rlmer(aucMA ~ Period + mtMA + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_aucLAe <- rlmer(aucLA ~ Period + mtLA + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_mtAIe <- rlmer(mtAI ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_aucAIe <- rlmer(aucAI ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_silAIe <- rlmer(silAI ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)
mod_ampAIe <- rlmer(ampAI ~ Period + (1|ID), data = TMStable_ePD, REML = TRUE)

round(rlmer_aprox_pval(mod_mtMAe,"Period"),3)
round(rlmer_aprox_pval(mod_mtLAe,"Period"),3)
round(rlmer_aprox_pval(mod_latMAe,"Period"),3)
round(rlmer_aprox_pval(mod_latLAe,"Period"),3)
round(rlmer_aprox_pval(mod_durMAe,"Period"),3)
round(rlmer_aprox_pval(mod_durLAe,"Period"),3)
round(rlmer_aprox_pval(mod_silMAe,"Period"),3)
round(rlmer_aprox_pval(mod_silLAe,"Period"),3)
round(rlmer_aprox_pval(mod_ampMAe,"Period"),3)
round(rlmer_aprox_pval(mod_ampLAe,"Period"),3)
round(rlmer_aprox_pval(mod_aucMAe,"Period"),3)
round(rlmer_aprox_pval(mod_aucLAe,"Period"),3)
round(rlmer_aprox_pval(mod_mtAIe,"Period"),3)
round(rlmer_aprox_pval(mod_aucAIe,"Period"),3)
round(rlmer_aprox_pval(mod_silAIe,"Period"),3)
round(rlmer_aprox_pval(mod_ampAIe,"Period"),3)

### Plots ###
pmtMAe <- plot_TMS(dataset = TMStable_ePD, metric = mtMA, model = mod_mtMAe, strTitle = "", yLabel = "%MSA", color2, c(0,150,20,10))# c(0,140,20,10)
pmtLAe <- plot_TMS(dataset = TMStable_ePD, metric = mtLA, model = mod_mtLAe, strTitle = "", yLabel = "%MSA", color2, c(0,150,20,10))# c(0,140,20,10)
platMAe <- plot_TMS(dataset = TMStable_ePD, metric = latMA, model = mod_latMAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,30,10,5))
platLAe <- plot_TMS(dataset = TMStable_ePD, metric = latLA, model = mod_latLAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,30,10,5))
pdurMAe <- plot_TMS(dataset = TMStable_ePD, metric = durMA, model = mod_durMAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,20,10,5))
pdurLAe <- plot_TMS(dataset = TMStable_ePD, metric = durLA, model = mod_durLAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,20,10,5))
psilMAe <- plot_TMS(dataset = TMStable_ePD, metric = silMA, model = mod_silMAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,180,20,10))# c(0,140,20,10)
psilLAe <- plot_TMS(dataset = TMStable_ePD, metric = silLA, model = mod_silLAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,180,20,10))# c(0,140,20,10)
pampMAe <- plot_TMS(dataset = TMStable_ePD, metric = ampMA, model = mod_ampMAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,12,3,1))
pampLAe <- plot_TMS(dataset = TMStable_ePD, metric = ampLA, model = mod_ampLAe, strTitle = "", yLabel = "Time (ms)", color2, c(0,12,3,1))
paucMAe <- plot_TMS(dataset = TMStable_ePD, metric = aucMA, model = mod_aucMAe, strTitle = "", yLabel = "AUC (mV∙ms)", color2, c(0,70,20,10))# c(0,70,20,10))
paucLAe <- plot_TMS(dataset = TMStable_ePD, metric = aucLA, model = mod_aucLAe, strTitle = "", yLabel = "AUC (mV∙ms)", color2, c(0,70,20,10))# c(0,70,20,10))
pmtAIe <- plot_TMS(dataset = TMStable_ePD, metric = mtAI, model = mod_mtAIe, strTitle = "", yLabel = "Index", color2, c(-1,1,0.5,0.25))
paucAIe <- plot_TMS(dataset = TMStable_ePD, metric = aucAI, model = mod_aucAIe, strTitle = "", yLabel = "Index", color2, c(-1,1,0.5,0.25))
psilAIe <- plot_TMS(dataset = TMStable_ePD, metric = silAI, model = mod_silAIe, strTitle = "", yLabel = "Index", color2, c(-1,1,0.5,0.25))
pampAIe <- plot_TMS(dataset = TMStable_ePD, metric = ampAI, model = mod_ampAIe, strTitle = "", yLabel = "Index", color2, c(-1,1,0.5,0.25))

############################
### Advanced PD subgroup ###
############################

### Filter advanced PD subjects (H&Y 2.5-4)
subSample <- c(4,7,8,11,19,21,22)
TMStable_aPD <- TMStable[TMStable$ID %in% subSample,]

### Models ###
mod_mtMAa <- rlmer(mtMA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_mtLAa <- rlmer(mtLA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_latMAa <- rlmer(latMA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_latLAa <- rlmer(latLA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_durMAa <- rlmer(durMA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_durLAa <- rlmer(durLA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_silMAa <- rlmer(silMA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_silLAa <- rlmer(silLA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_ampMAa <- rlmer(ampMA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_ampLAa <- rlmer(ampLA ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_aucMAa <- rlmer(aucMA ~ Period + mtMA + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_aucLAa <- rlmer(aucLA ~ Period + mtLA + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_mtAIa <- rlmer(mtAI ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_aucAIa <- rlmer(aucAI ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_silAIa <- rlmer(silAI ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)
mod_ampAIa <- rlmer(ampAI ~ Period + (1|ID), data = TMStable_aPD, REML = TRUE)

round(rlmer_aprox_pval(mod_mtMAa,"Period"),3)
round(rlmer_aprox_pval(mod_mtLAa,"Period"),3)
round(rlmer_aprox_pval(mod_latMAa,"Period"),3)
round(rlmer_aprox_pval(mod_latLAa,"Period"),3)
round(rlmer_aprox_pval(mod_durMAa,"Period"),3)
round(rlmer_aprox_pval(mod_durLAa,"Period"),3)
round(rlmer_aprox_pval(mod_silMAa,"Period"),3)
round(rlmer_aprox_pval(mod_silLAa,"Period"),3)
round(rlmer_aprox_pval(mod_ampMAa,"Period"),3)
round(rlmer_aprox_pval(mod_ampLAa,"Period"),3)
round(rlmer_aprox_pval(mod_aucMAa,"Period"),3)
round(rlmer_aprox_pval(mod_aucLAa,"Period"),3)
round(rlmer_aprox_pval(mod_mtAIa,"Period"),3)
round(rlmer_aprox_pval(mod_aucAIa,"Period"),3)
round(rlmer_aprox_pval(mod_silAIa,"Period"),3)
round(rlmer_aprox_pval(mod_ampAIa,"Period"),3)

### Plots ###
pmtMAa <- plot_TMS(dataset = TMStable_aPD, metric = mtMA, model = mod_mtMAa, strTitle = "", yLabel = "%MSA", color3, c(0,150,20,10))# c(0,140,20,10)
pmtLAa <- plot_TMS(dataset = TMStable_aPD, metric = mtLA, model = mod_mtLAa, strTitle = "", yLabel = "%MSA", color3, c(0,150,20,10))# c(0,140,20,10)
platMAa <- plot_TMS(dataset = TMStable_aPD, metric = latMA, model = mod_latMAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,30,10,5))
platLAa <- plot_TMS(dataset = TMStable_aPD, metric = latLA, model = mod_latLAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,30,10,5))
pdurMAa <- plot_TMS(dataset = TMStable_aPD, metric = durMA, model = mod_durMAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,20,10,5))
pdurLAa <- plot_TMS(dataset = TMStable_aPD, metric = durLA, model = mod_durLAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,20,10,5))
psilMAa <- plot_TMS(dataset = TMStable_aPD, metric = silMA, model = mod_silMAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,180,20,10))# c(0,140,20,10)
psilLAa <- plot_TMS(dataset = TMStable_aPD, metric = silLA, model = mod_silLAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,180,20,10))# c(0,140,20,10)
pampMAa <- plot_TMS(dataset = TMStable_aPD, metric = ampMA, model = mod_ampMAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,12,3,1))
pampLAa <- plot_TMS(dataset = TMStable_aPD, metric = ampLA, model = mod_ampLAa, strTitle = "", yLabel = "Time (ms)", color3, c(0,12,3,1))
paucMAa <- plot_TMS(dataset = TMStable_aPD, metric = aucMA, model = mod_aucMAa, strTitle = "", yLabel = "AUC (mV∙ms)", color3, c(0,70,20,10))# c(0,70,20,10))
paucLAa <- plot_TMS(dataset = TMStable_aPD, metric = aucLA, model = mod_aucLAa, strTitle = "", yLabel = "AUC (mV∙ms)", color3, c(0,70,20,10))# c(0,70,20,10))
pmtAIa <- plot_TMS(dataset = TMStable_aPD, metric = mtAI, model = mod_mtAIa, strTitle = "", yLabel = "Index", color3, c(-1,1,0.5,0.25))
paucAIa <- plot_TMS(dataset = TMStable_aPD, metric = aucAI, model = mod_aucAIa, strTitle = "", yLabel = "Index", color3, c(-1,1,0.5,0.25))
psilAIa <- plot_TMS(dataset = TMStable_aPD, metric = silAI, model = mod_silAIa, strTitle = "", yLabel = "Index", color3, c(-1,1,0.5,0.25))
pampAIa <- plot_TMS(dataset = TMStable_aPD, metric = ampAI, model = mod_ampAIa, strTitle = "", yLabel = "Index", color3, c(-1,1,0.5,0.25))


# Fig. 3: 1500 x 1200
t1 <- ggdraw() + draw_label("Resting Motor Threshold\n(Most affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("Resting Motor Threshold\n(Least affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t3 <- ggdraw() + draw_label("Cortical Silent Period\n(Most affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("Cortical Silent Period\n(Least affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color1)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color2)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color3)
plot_grid(NULL, t1, t2, t3, t4, 
          L1, pmtMA, pmtLA, psilMA, psilLA, 
          L2, pmtMAe, pmtLAe, psilMAe, psilLAe, 
          L3, pmtMAa, pmtLAa, psilMAa, psilLAa, 
          ncol = 5, rel_heights = c(0.15,1,1,1), rel_widths = c(0.1,1,1,1,1))

# Sup. Fig. 2: 1500 x 1200
t1 <- ggdraw() + draw_label("Latency\n(Most affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("Latency\n(Least affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t3 <- ggdraw() + draw_label("Duration\n(Most affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("Duration\n(Least affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color1)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color2)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color3)
plot_grid(NULL, t1, t2, t3, t4, 
          L1, platMA, platLA, pdurMA, pdurLA, 
          L2, platMAe, platLAe, pdurMAe, pdurLAe, 
          L3, platMAa, platLAa, pdurMAa, pdurLAa,
          ncol = 5, rel_heights = c(0.15,1,1,1), rel_widths = c(0.1,1,1,1,1))

# Sup. Fig. 3: 1500 x 1200
t1 <- ggdraw() + draw_label("MEP Amplitude @130%rMT\n(Most affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("MEP Amplitude @130%rMT\n(Least affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t3 <- ggdraw() + draw_label("Area Under the Curve\n(Most affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("Area Under the Curve\n(Least affected hemisphere)", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color1)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color2)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color3)
plot_grid(NULL, t1, t2, t3, t4, 
          L1, pampMA, pampLA, paucMA, paucLA, 
          L2, pampMAe, pampLAe, paucMAe, paucLAe, 
          L3, pampMAa, pampLAa, paucMAa, paucLAa, 
          ncol = 5, rel_heights = c(0.15,1,1,1), rel_widths = c(0.1,1,1,1,1))

# Sup. Fig. 4: 1500 x 1200
t1 <- ggdraw() + draw_label("Resting Motor Threshold\nAsymmetry index", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("Area Under the Curve\nAsymmetry index", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t3 <- ggdraw() + draw_label("Cortical Silent Period\nAsymmetry index", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("MEP Amplitude @130%rMT\nAsymmetry index", fontface = 'bold', size = 18, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color1)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color2)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = color3)
plot_grid(NULL, t1, t2, t3, t4, 
          L1, pmtAI, paucAI, psilAI, pampAI, 
          L2, pmtAIe, paucAIe, psilAIe, pampAIe, 
          L3, pmtAIa, paucAIa, psilAIa, pampAIa, 
          ncol = 5, rel_heights = c(0.15,1,1,1), rel_widths = c(0.1,1,1,1,1))


### Analysis by Sex ###

### Models ###
mod_mtMA <- rlmer(mtMA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_mtLA <- rlmer(mtLA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_latMA <- rlmer(latMA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_latLA <- rlmer(latLA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_durMA <- rlmer(durMA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_durLA <- rlmer(durLA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_silMA <- rlmer(silMA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_silLA <- rlmer(silLA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_ampMA <- rlmer(ampMA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_ampLA <- rlmer(ampLA ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_aucMA <- rlmer(aucMA ~ Period + Sex + mtMA + (1|ID), data = TMStable, REML = TRUE)
mod_aucLA <- rlmer(aucLA ~ Period + Sex + mtLA + (1|ID), data = TMStable, REML = TRUE)
mod_mtAI <- rlmer(mtAI ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_aucAI <- rlmer(aucAI ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_silAI <- rlmer(silAI ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)
mod_ampAI <- rlmer(ampAI ~ Period + Sex + (1|ID), data = TMStable, REML = TRUE)

rlmer_aprox_pval(mod_mtMA ,"Sex")
rlmer_aprox_pval(mod_mtLA ,"Sex")
rlmer_aprox_pval(mod_latMA ,"Sex")
rlmer_aprox_pval(mod_latLA ,"Sex")
rlmer_aprox_pval(mod_durMA ,"Sex")
rlmer_aprox_pval(mod_durLA ,"Sex")
rlmer_aprox_pval(mod_silMA ,"Sex")
rlmer_aprox_pval(mod_silLA ,"Sex")
rlmer_aprox_pval(mod_ampMA ,"Sex")
rlmer_aprox_pval(mod_ampLA ,"Sex")
rlmer_aprox_pval(mod_aucMA ,"Sex")
rlmer_aprox_pval(mod_aucLA ,"Sex")
rlmer_aprox_pval(mod_mtAI ,"Sex")
rlmer_aprox_pval(mod_aucAI ,"Sex")
rlmer_aprox_pval(mod_silAI ,"Sex")
rlmer_aprox_pval(mod_ampAI ,"Sex")


### Analysis by Stage ###

### Models ###
mod_mtMA <- rlmer(mtMA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_mtLA <- rlmer(mtLA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_latMA <- rlmer(latMA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_latLA <- rlmer(latLA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_durMA <- rlmer(durMA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_durLA <- rlmer(durLA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_silMA <- rlmer(silMA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_silLA <- rlmer(silLA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_ampMA <- rlmer(ampMA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_ampLA <- rlmer(ampLA ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_aucMA <- rlmer(aucMA ~ Period + Stage + mtMA + (1|ID), data = TMStable, REML = TRUE)
mod_aucLA <- rlmer(aucLA ~ Period + Stage + mtLA + (1|ID), data = TMStable, REML = TRUE)
mod_mtAI <- rlmer(mtAI ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_aucAI <- rlmer(aucAI ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_silAI <- rlmer(silAI ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)
mod_ampAI <- rlmer(ampAI ~ Period + Stage + (1|ID), data = TMStable, REML = TRUE)

rlmer_aprox_pval(mod_mtMA ,"Stage")
rlmer_aprox_pval(mod_mtMA ,"Stage")
rlmer_aprox_pval(mod_mtLA ,"Stage")
rlmer_aprox_pval(mod_latMA ,"Stage")
rlmer_aprox_pval(mod_latLA ,"Stage")
rlmer_aprox_pval(mod_durMA ,"Stage")
rlmer_aprox_pval(mod_durLA ,"Stage")
rlmer_aprox_pval(mod_silMA ,"Stage")
rlmer_aprox_pval(mod_silLA ,"Stage")
rlmer_aprox_pval(mod_ampMA ,"Stage")
rlmer_aprox_pval(mod_ampLA ,"Stage")
rlmer_aprox_pval(mod_aucMA ,"Stage")
rlmer_aprox_pval(mod_aucLA ,"Stage")
rlmer_aprox_pval(mod_mtAI ,"Stage")
rlmer_aprox_pval(mod_aucAI ,"Stage")
rlmer_aprox_pval(mod_silAI ,"Stage")
