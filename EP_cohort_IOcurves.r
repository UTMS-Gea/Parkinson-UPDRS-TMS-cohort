### Loads the data and creates data frames for the analysis ###
fpath <- paste(getwd(),"/Parkinson-UPDRS-TMS-cohort-main/data_UDPRS_TMS/",sep = "")
RCtable <- read.table(file = paste(fpath,"IOCurves.csv",sep = ""), sep = ",", header = TRUE)
RCtable$Period <- factor(RCtable$Period)
RCtable$ID <- factor(RCtable$ID)
RCtable$Intensity <- factor(RCtable$Intensity)
RCtable2 <- data.frame(ID = rep(RCtable$ID,2),
                       Period = rep(RCtable$Period,2),
                       Intensity = rep(RCtable$Intensity,2),
                       Side = c(rep("MA",nrow(RCtable)), rep("LA",nrow(RCtable))),
                       Amplitude = c(RCtable$ampMA, RCtable$ampLA))
colorA <- "#00AFBB"
colorB <- "#5594c9"
colorC <- "#f8766d"
colorA1 <- "#00D0DE"
colorA2 <- "#00929C"
colorB1 <- "#62ABE8"
colorB2 <- "#4980AD"
colorC1 <- "#FF9393"
colorC2 <- "#D6665E"
ic50 <- data.frame(LA = numeric(12), MA = numeric(12))

### T1, All ###
RCtableT1 <- subset(RCtable2, Period=="T1")

x <- as.numeric(as.character(unlist(subset(RCtableT1, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT1, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[1] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT1, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT1, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[1] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T1_All <- ggerrorplot(RCtableT1, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorA1, colorA2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorA1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorA2, 
            size = 0.5)

### T2, All ###
RCtableT2 <- subset(RCtable2, Period=="T2")

x <- as.numeric(as.character(unlist(subset(RCtableT2, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT2, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[2] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT2, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT2, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[2] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T2_All <- ggerrorplot(RCtableT2, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorA1, colorA2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorA1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorA2, 
            size = 0.5)

### T4, All ###
RCtableT4 <- subset(RCtable2, Period=="T4")

x <- as.numeric(as.character(unlist(subset(RCtableT4, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT4, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[3] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT4, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT4, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[3] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T4_All <- ggerrorplot(RCtableT4, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorA1, colorA2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorA1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorA2, 
            size = 0.5)

### T5, All ###
RCtableT5 <- subset(RCtable2, Period=="T5")

x <- as.numeric(as.character(unlist(subset(RCtableT5, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT5, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[4] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT5, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT5, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[4] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T5_All <- ggerrorplot(RCtableT5, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorA1, colorA2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorA1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorA2, 
            size = 0.5)

### Filter early PD subjects (H&Y 1-2)
subMuestra <- c(1,2,3,5,6,9,10,12,13,14,15,16,17,18,20)
RCtable_ePD <- RCtable2[RCtable2$ID %in% subMuestra,]

### T1, early ###
RCtableT1 <- subset(RCtable_ePD, Period=="T1")

x <- as.numeric(as.character(unlist(subset(RCtableT1, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT1, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[5] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT1, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT1, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[5] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T1_e <- ggerrorplot(RCtableT1, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorB1, colorB2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorB1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorB2, 
            size = 0.5)

### T2, early ###
RCtableT2 <- subset(RCtable_ePD, Period=="T2")

x <- as.numeric(as.character(unlist(subset(RCtableT2, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT2, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[6] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT2, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT2, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[6] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T2_e <- ggerrorplot(RCtableT2, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorB1, colorB2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorB1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorB2, 
            size = 0.5)

### T4, early ###
RCtableT4 <- subset(RCtable_ePD, Period=="T4")

x <- as.numeric(as.character(unlist(subset(RCtableT4, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT4, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[7] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT4, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT4, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[7] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T4_e <- ggerrorplot(RCtableT4, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorB1, colorB2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorB1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorB2, 
            size = 0.5)

### T5, early ###
RCtableT5 <- subset(RCtable_ePD, Period=="T5")

x <- as.numeric(as.character(unlist(subset(RCtableT5, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT5, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[8] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT5, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT5, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[8] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T5_e <- ggerrorplot(RCtableT5, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorB1, colorB2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorB1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorB2, 
            size = 0.5)

### Filter advanced PD subjects (H&Y 1-2)
subMuestra <- c(4,7,8,11,19,21,22)
RCtable_aPD <- RCtable2[RCtable2$ID %in% subMuestra,]

### T1, advanced ###
RCtableT1 <- subset(RCtable_aPD, Period=="T1")

x <- as.numeric(as.character(unlist(subset(RCtableT1, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT1, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[9] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT1, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT1, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[9] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T1_a <- ggerrorplot(RCtableT1, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorC1, colorC2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorC1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorC2, 
            size = 0.5)

### T2, advanced ###
RCtableT2 <- subset(RCtable_aPD, Period=="T2")

x <- as.numeric(as.character(unlist(subset(RCtableT2, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT2, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[10] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT2, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT2, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[10] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T2_a <- ggerrorplot(RCtableT2, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorC1, colorC2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorC1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorC2, 
            size = 0.5)

### T4, advanced ###
RCtableT4 <- subset(RCtable_aPD, Period=="T4")

x <- as.numeric(as.character(unlist(subset(RCtableT4, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT4, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[11] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT4, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT4, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[11] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T4_a <- ggerrorplot(RCtableT4, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorC1, colorC2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorC1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorC2, 
            size = 0.5)

### T5, advanced ###
RCtableT5 <- subset(RCtable_aPD, Period=="T5")

x <- as.numeric(as.character(unlist(subset(RCtableT5, Side == "LA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT5, Side == "LA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$LA[12] <- coef(fit)["c"]
LAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

x <- as.numeric(as.character(unlist(subset(RCtableT5, Side == "MA", select = "Intensity"))))
y <- as.numeric(unlist(subset(RCtableT5, Side == "MA", select = "Amplitude")))
fit <- nlsLM(y ~ b/(1+exp((c-x)/d)), 
             start = list(b = 0, c = 140, d = 20),
             lower = c(0, 0, 0),
             upper = c(15, 180, 1000))
ic50$MA[12] <- coef(fit)["c"]
MAy <- predict(fit, newdata = list(x = seq(100, 180, by = 1)))

T5_a <- ggerrorplot(RCtableT5, x="Intensity", y="Amplitude", desc_stat = "mean_ci", color = "Side", error.plot = "errorbar", add = "mean", position = position_dodge(width = 0)) +
  xlab("Intensity (%rMT)") +
  ylab("Amplitude (mV)") +
  scale_y_continuous(breaks = seq(0,9,2), minor_breaks = seq(0,9,1)) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_x_discrete(breaks = seq(100, 180, by = 20), labels = seq(100, 180, by = 20)) +
  scale_color_manual(name = "Hemisphere", labels = c("Less affected","More affected"), values = c(colorC1, colorC2)) +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 45, size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(1,1,1,0), "lines"), legend.position = "none") +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = LAy), 
            aes(x = x, y = y), 
            color = colorC1, 
            size = 0.5) +
  geom_line(data = data.frame(x = seq(1, 9, by = 0.1), y = MAy), 
            aes(x = x, y = y), 
            color = colorC2, 
            size = 0.5)

ic50
t1 <- ggdraw() + draw_label("2018", fontface = 'bold', size = 16, hjust = 0.5, vjust = 0.5)
t2 <- ggdraw() + draw_label("2019", fontface = 'bold', size = 16, hjust = 0.5, vjust = 0.5)
t4 <- ggdraw() + draw_label("2021", fontface = 'bold', size = 16, hjust = 0.5, vjust = 0.5)
t5 <- ggdraw() + draw_label("2022", fontface = 'bold', size = 16, hjust = 0.5, vjust = 0.5)
L1 <- ggdraw() + draw_label("All PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = colorA)
L2 <- ggdraw() + draw_label("Early PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = colorB)
L3 <- ggdraw() + draw_label("Advanced PD", fontface = 'bold', size = 18, angle = 90, hjust = 0.5, vjust = 0.5, color = colorC)

# Fig. 3: 1000x1000
plot_grid(NULL, t1, t2, t4, t5,
          L1, T1_All, T2_All, T4_All, T5_All, 
          L2, T1_e, T2_e, T4_e, T5_e, 
          L3, T1_a, T2_a, T4_a, T5_a, 
          ncol = 5, rel_heights = c(0.1,1,1,1), rel_widths = c(0.1,1,1,1,1))
