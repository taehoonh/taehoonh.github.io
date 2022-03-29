# This R script is for Core Knowledge Seminar on 03/30/2022 about Biostatistics-ANOVA by Taehoon Ha.

# Load R packages we need today
library(dplyr)
library(readxl)
library(multcomp)
library(ggplot2)
library(lme4)
library(lmerTest)

# For those who do not have packages above, please uncomment the following line and install packages.
# install.packages(c("dplyr", "readxl", "multcomp", "ggplot2", "lme4", "lmerTest"))

# Today's dataset - Diet (low fat vs. high fat) and Ovx (shamOVX vs. OVX)
# Place the dataset and this R script file in the same folder and click [Session]-[Set working directory]-[To source file location] on the menu.
dat.work <- read_xlsx("CoreKnowledge.xlsx", sheet = 1)

# Preprocessing Factorize the categorical variable as preprocessing
dat.work$GroupName <- factor(dat.work$GroupName, levels = c("LowFat", "LowFat.OVX", "HighFat", "HighFat.OVX"))


#########################################################
# [1] One-way ANOVA

# 1. One-way ANOVA using cell means model
cellmeans_model <- lm(MouseWt ~ GroupName -1, data = dat.work)
summary(cellmeans_model)



# 2. As anova p-value is significant, we should move on to post-hoc analysis to identify which group mean(s) are significant.
anova(cellmeans_model)



# 3. Post-hoc analysis
aov.cellmeans <- aov(MouseWt ~ GroupName -1, data = dat.work)



# 3-1. Example 1: Tukey's method - compare all possible pairs
summary(glht(aov.cellmeans, linfct=mcp(GroupName="Tukey")), test = adjusted("none"))

# 3-2. Example 2: Dunnett's method - compare control to each treatment group. 
# Here, LowFat group is the reference, so compare LowFat to each of other groups.
summary(glht(aov.cellmeans, linfct=mcp(GroupName="Dunnett")), test = adjusted("none"))

# 3-3. Example 3: Contrast Matrix K
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



# 4. Visualizing the group differences using boxplot
boxplot(MouseWt ~ GroupName, data = dat.work,
        xlab = "Group", ylab = "Mouse Weight",
        main = "Average Mouse Weight by Diet and OVX")



# 5. To see if we can improve our model - using weighted linear regression could improve the model fitting: inverse variance weighting
wts <- 1/rep(tapply(dat.work$MouseWt, dat.work$GroupName, var), each = 10)
cellmeans_model <- lm(MouseWt ~ GroupName - 1, data = dat.work, weights = wts)



# 6. Model diagnostics
# 6-1. assumption check (1) - Normality
qqnorm(rstudent(cellmeans_model))
qqline(rstudent(cellmeans_model))

# 6-2. assumption check (2) - Homogeneity of variances
plot(1:nrow(dat.work), rstudent(cellmeans_model), pch = 3,
     xlab="Index", ylab="Studentized Residual")

# > Both normality and homogeneity of variances assumptions were met!

#########################################################

# [2] Two-way ANOVA

# 1. Preprocessing the dataset
dat.work$Diet <- factor(dat.work$Diet, levels = c("LowFat", "HighFat"))
dat.work$OVX <- factor(dat.work$OVX, levels = c("shamOVX", "OVX"))

# 2. Two-way ANOVA using effects model
effects_model <- lm(MouseWt ~ Diet * OVX, data = dat.work)
summary(effects_model)

# 3. According to the ANOVA result, both Diet, Ovx, and the interaction of Diet and Ovx are all significant.
anova(effects_model)

# 4. Post-hoc analysis
# - 1st column: LowFat group mean (reference group)
# - 2nd column: Difference in mean between LowFat and HighFat
# - 3rd column: Difference in mean between LowFat.OVX and LowFat
# - 4th column: (HF.OVX - HF) - (LF.OVX - LF)
# > As cell means model and effects model are different in desinging contrast matrix, we need to precisely know what it means.
K2 <- rbind("LF.OVX - LF"               = c(0, 0, 1, 0),
            "HF.OVX - HF"               = c(0, 0, 1, 1),
            "HF - LF"                   = c(0, 1, 0, 0),
            "HF.OVX - LF.OVX"           = c(0, 1, 0, 1),
            "LF.OVX + HF.OVX - LF - HF" = c(0, 0, 2, 1),
            "HF + HF.OVX - LF - LF.OVX" = c(0, 2, 0, 1), 
            "HF.OVX - HF - LF.OVX + LF" = c(0, 0, 0, 1))
summary(glht(effects_model, linfct=K2), test=adjusted(type="none"))



# 5. Visualizing the interaction using line graph
# > Depending on the diet type and OVX, the slopes are different with each other. This is because Diet is significantly associated with Ovx.
interaction.plot(dat.work$Diet, dat.work$OVX, dat.work$MouseWt, 
                 xlab="Diet Type", ylab="Average Mouse Weight",
                 legend=F, lty=2:1)
legend("topleft", legend=levels(dat.work$OVX), lty=2:1, bty="n")

#########################################################

# [3] Partial F-test
# - Goal: To compare a full (complex) model to a reduced (simpler) model.
# - Caution: A full (complex) model should include all the predictors that a reduced (simpler) model has. ["nested model"]
# - To compare between non-nested models, use AIC, BIC, or Vuong's test.

reduced_model <- lm(MouseWt ~ Diet, data = dat.work)
full_model <- lm(MouseWt ~ Diet + OVX + Diet * OVX, data = dat.work)
anova(reduced_model, full_model)
# > As p<0.05, there is enough evidence that the full (complex) model is better than the reduced (simpler) model.

#########################################################

# [4] Mixed-effects model
# - Limitation of repeated measures ANOVA: if there is one sinlge NA value in the dataset, that sample is not available.
# - Solution: imputation.
# - Alternatively, we can try mixed effects model instead.

# 1. Today's dataset2 - longitudinal blood pressure between drug types
blood_data <- read_xlsx("CoreKnowledge.xlsx", sheet = 2) %>%
    mutate(id = factor(id),
           trt = factor(trt))



# 2. First 12 patients' bp change over time
# - Current issue: blood pressure change pattern differs from patient to patient.
# - Can we better measure the treatment effect and minimize the patients' effect (variations)?
ggplot(blood_data %>% filter(id %in% 1:12), aes(x=time, y=bp)) +
    geom_point(size=0.5) +
    stat_smooth(method = "lm",se=F,size=0.5)+
    facet_wrap("id", labeller = label_both)+ 
    theme_bw() +
    ggtitle("First 12 Patients' BP") +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

# 3. Interaction plot between Drug type and Time: Drug is associated with Time.
blood_data$time <- factor(blood_data$time)
with(blood_data, interaction.plot(x.factor = time, trace.factor = trt, response = d,
                                  xlab = "Day", ylab = "Difference of BP", 
                                  main = "Interaction plot"))

# 4. Mixed effects modeling
blood_fit <- lmer(d ~ trt + time + trt * time + (1|id), data = blood_data)
summary(blood_fit)
# > Treatment effect was significant. However, control drug decreased more blood pressure.
# > Time effect was not significant.
# > Interaction between treatment and time was not significant. 
# > Overall, we can conclude that not enough evidence was shown that the new (treatment) drug significantly decreased the blood pressure than the control (current) drug.
