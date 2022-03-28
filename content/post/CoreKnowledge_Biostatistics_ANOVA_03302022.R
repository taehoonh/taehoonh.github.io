# This R script is for Core Knowledge Seminar on 03/30/2022 about Biostatistics-ANOVA by Taehoon Ha.

# R packages we need
library(dplyr)
library(readxl)
library(multcomp)
library(ggplot2)
library(lme4)
library(lmerTest)

# For those who do not have packages above, please uncomment the following codes and install packages.
# install.packages(c("dplyr", "readxl", "multcomp"))

# Today's dataset - Diet and Ovx
dat.work <- read_xlsx("CoreKnowledge.xlsx", sheet = 1)
# Factorize the categorical variable as preprocessing
dat.work$GroupName <- factor(dat.work$GroupName, levels = c("LowFat", "LowFat.OVX", "HighFat", "HighFat.OVX"))

#########################################################
# One-way ANOVA using cell means model
cellmeans_model <- lm(MouseWt ~ GroupName -1, data = dat.work)
summary(cellmeans_model)

# As anova p-value is significant, we should move on to post-hoc analysis to identify which group mean(s) are significant.
anova(cellmeans_model)

# Post-hoc analysis
aov.cellmeans <- aov(MouseWt ~ GroupName -1, data = dat.work)
summary(glht(aov.cellmeans, linfct=mcp(GroupName="Tukey")), test = adjusted("none"))

# In the contrast matrix K,
# - 1st column: LowFat.shamOVX group
# - 2nd column: LowFat.OVX group
# - 3rd column: HighFat.shamOVX group
# - 4th column: HighFat.OVX group
K <- rbind("OVX effect in LF"     = c(-1, 1, 0, 0),
           "OVX effect in HF"     = c(0, 0, -1, 1),
           "HF effect in shamOVX" = c(-1, 0, 1, 0),
           "HF effect in OVX"     = c(0, -1, 0, 1),
           "OVX effect"           = c(-1, 1, -1, 1),
           "HF effect"            = c(-1, -1, 1, 1), 
           "OVX HF Interation"    = c(1, -1, -1, 1))
summary(glht(aov.cellmeans, linfct=mcp(GroupName=K)), test = adjusted(type="none"))

# Visualizing the group differences using boxplot
boxplot(MouseWt ~ GroupName, data = dat.work,
        xlab = "Group", ylab = "Mouse Weight",
        main = "Average Mouse Weight by Diet and OVX")

# Weighted linear regression could make the data better fit: inverse variance weighting
wts <- 1/rep(tapply(dat.work$MouseWt, dat.work$GroupName, var),each=10)
cellmeans_model <- lm(MouseWt~ GroupName-1, data=dat.work, weights=wts)

# Model diagnostics: assumption check (1) - normality
qqnorm(rstudent(cellmeans_model))
qqline(rstudent(cellmeans_model))

# Model diagnostics: assumption check (2) - homegeneity of variances
plot(1:nrow(dat.work), rstudent(cellmeans_model), pch = 3,
     xlab="Index", ylab="Studentized Residual")

#########################################################

# Two-way ANOVA using effects model
dat.work$Diet <- factor(dat.work$Diet, levels = c("LowFat", "HighFat"))
dat.work$OVX <- factor(dat.work$OVX, levels = c("shamOVX", "OVX"))
effects_model <- lm(MouseWt ~ Diet * OVX, data = dat.work)
summary(effects_model)
anova(effects_model)
K2 <- rbind("LF.OVX - LF" = c(0,0,1,0),
            "HF.OVX - HF" = c(0,0,1,1),
            "HF - LF" = c(0,1,0,0),
            "HF.OVX - LF.OVX" = c(0,1,0,1),
            "LF.OVX + HF.OVX - LF - HF" = c(0,0,2,1),
            "HF + HF.OVX - LF - LF.OVX" = c(0,2,0,1), 
            "HF.OVX - HF - LF.OVX + LF"= c(0, 0, 0, 1))
summary(glht(effects_model, linfct=K2), test=adjusted(type="none"))

# Interaction plot between Diet and Ovx
interaction.plot(dat.work$Diet, dat.work$OVX, dat.work$MouseWt, 
                 xlab="Diet Type", ylab="Average Mouse Weight",
                 legend=F, lty=2:1)
legend("topleft", legend=levels(dat.work$OVX), lty=2:1, bty="n")

#########################################################

# Partial F-test
reduced_model <- lm(MouseWt ~ Diet, data = dat.work)
full_model <- lm(MouseWt ~ Diet + OVX + Diet * OVX, data = dat.work)
anova(reduced_model, full_model)

#########################################################

# Mixed-effects model

# Today's dataset - longitudinal blood pressure between drug types
blood_data <- read_xlsx("CoreKnowledge.xlsx", sheet = 2) %>%
    mutate(id = factor(id),
           trt = factor(trt))

# First 12 patients' bp change over time
ggplot(blood_data %>% filter(id %in% 1:12), aes(x=time, y=bp)) +
    geom_point(size=0.5) +
    stat_smooth(method = "lm",se=F,size=0.5)+
    facet_wrap("id", labeller = label_both)+ 
    theme_bw() +
    ggtitle("First 12 Patients' BP") +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

# Interaction plot between Drug type and Time
blood_data$time <- factor(blood_data$time)
with(blood_data, interaction.plot(x.factor = time, trace.factor = trt, response = d,
                                  xlab = "Day", ylab = "Difference of BP", 
                                  main = "Interaction plot"))

# Mixed effects modeling
blood_fit <- lmer(d ~ trt + time + trt * time + (1|id), data = blood_data)
summary(blood_fit)