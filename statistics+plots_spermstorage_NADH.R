#############################################################################
## statistical analysis of NAD(P)H FLIM data ################################ 
## author: Cornelia Wetzker #################################################
## license: CC-BY 4.0 International #########################################
#############################################################################
#
library(lme4)
library(nlme)
library(car)
library(ggplot2)
library(Rmisc)
#
rm(list=ls()) 
#
### DATA IMPORT
flimN<-read.csv("D:/path_to_data/NADH_LT_data_table.csv")
str(flimN)
#
### DATA FORMATTING
#
#factorization of fluorophore, sex, cohort, mating state, magnification being zoom factor, limited photon counts, exclusion of sample based on quality thresholds
#
flimN$fluor<-as.factor(flimN$fluor)
flimN$sex<-as.factor(flimN$sex)
flimN$coh<-as.factor(flimN$coh)
levels(flimN$coh) <- c("mated females", "mated males", "unmated males")
flimN$mat<-as.factor(flimN$mat)
flimN$magn<-as.factor(flimN$magn)
flimN$lim_phot_number<-as.factor(flimN$lim_phot_number)
flimN$excl<-as.factor(flimN$excl)
str(flimN)

#Exclusion of samples based on quality measures (data quality and photon counts)
#
flimxN <- subset(flimN, lim_phot_number=="no")
flimxN <- subset(flimxN, excl=="n")
str(flimxN)
#
# Subselection of high magnification images
flimxN_high <- subset(flimxN, magn=="high")   
#
# Exclusion of unmated males, comparison of mated males and females
flimxN_high_mated <- subset(flimxN_high, coh != "unmated males")    
#
### STATISTICAL ANALYSES
#
# Test for normal distribution of tm (mated only)
hist(flimxN_high_mated$NADH_tm)     
shapiro.test(flimxN_high_mated$NADH_tm)   
#
# Test for normal distribution of a1% (mated only)
hist((flimxN_high_mated$NADH_a1))      
shapiro.test(flimxN_high_mated$NADH_a1)
#
# Test for homoscedasticity using Levene's test
mod_tm <- lmer(NADH_tm~sex*dpm+(1|sp_id), data=flimxN_high_mated)
leveneTest(residuals(mod_tm) ~ flimxN_high_mated$sex)
boxplot(residuals(mod_tm) ~ flimxN_high_mated$sex)
#
mod_a1 <- lmer(NADH_a1~sex*dpm+(1|sp_id), data=flimxN_high_mated)
leveneTest(residuals(mod_a1) ~ flimxN_high_mated$sex)
boxplot(residuals(mod_a1) ~ flimxN_high_mated$sex)
#
# Summary of statistics
#
summ_high_NADH_a1 <- summarySE(flimxN_high_mated, measurevar = "NADH_a1", groupvars = c("coh", "dpm"))   
summ_high_NADH_a1
summ_high_NADH_tm <- summarySE(flimxN_high_mated, measurevar = "NADH_tm", groupvars = c("coh", "dpm"))   
summ_high_NADH_tm
#
# Linear mixed model for NADH a1%
lmer_mated_a1_high <- with(flimxN_high_mated,lmer(NADH_a1~dpm*sex+(1|sp_id)))
Anova(lmer_mated_a1_high,type="III")
plot(lmer_mated_a1_high)
#
# Linear mixed model for NADH tm
lmer_mated_tm_high <- with(flimxN_high_mated,lmer(NADH_tm~dpm*sex+(1|sp_id)))
Anova(lmer_mated_tm_high,type="III")
plot(lmer_mated_tm_high)
#
### VISUALIZATION OF DATA 
#
# Plotting of a1% of mated males and females across days post mating
summ_high_mated <- summarySE(flimxN_high_mated, measurevar = "NADH_a1", groupvars = c("coh", "dpm"))   # for mean and SE for dataset
summ_high_mated
flimxN_high_fma <- subset(flimxN_high_mated, coh == "mated females")  
lm_high_fma <- with(flimxN_high_fma,lm(NADH_a1~dpm))    
lm_high_fma    
flimxN_high_mma <- subset(flimxN_high_mated, coh == "mated males")  
lm_high_mma <- with(flimxN_high_mma,lm(NADH_a1~dpm))    
lm_high_mma    
p <- ggplot(summ_high_mated, aes(x=dpm, y=NADH_a1, group=coh, color=coh, shape=coh)) + 
  geom_point(size=3)+
  scale_color_manual(values=c("#E69F00", "deepskyblue3")) +       
  scale_shape_manual(values=c(15, 17))+
  geom_errorbar(aes(ymin=NADH_a1-se, ymax=NADH_a1+se), width=.1,
                position=position_dodge(0.05)) +
  scale_x_continuous(name="Days post mating") +                                    
  scale_y_continuous(name="% short NADH LT fraction") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15)) +
  labs(title="NAD(P)H a1%: Mated males and females, high magn.") +
  geom_abline(intercept=80.7270, slope=-0.3689, color="#E69F00", linetype=2) +       
  geom_abline(intercept=81.07153, slope=-0.06159, color="deepskyblue3", linetype=2)  
print(p)
#
# Plotting of tm of mated males and females across days post mating
#
summ_high_mated <- summarySE(flimxN_high_mated, measurevar = "NADH_tm", groupvars = c("coh", "dpm"))   # generates mean and SE for dataset
summ_high
flimxN_high_fma <- subset(flimxN_high_mated, coh == "mated females")  
lm_high_fma <- with(flimxN_high_fma,lm(NADH_tm~dpm))    
lm_high_fma    
flimxN_high_mma <- subset(flimxN_high_mated, coh == "mated males")  
lm_high_mma <- with(flimxN_high_mma,lm(NADH_tm~dpm))    
lm_high_mma    
p <- ggplot(summ_high_mated, aes(x=dpm, y=NADH_tm, group=coh, color=coh, shape=coh)) + 
  geom_point(size=3)+
  scale_color_manual(values=c("#E69F00", "deepskyblue3")) +       
  scale_shape_manual(values=c(15, 17))+
  geom_errorbar(aes(ymin=NADH_tm-se, ymax=NADH_tm+se), width=.1,
                position=position_dodge(0.05)) +
  scale_x_continuous(name="Days post mating") +                                     
  scale_y_continuous(name="Mean NAD(P)H LT") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15)) +
  labs(title="Mean NAD(P)H LT: Mated males and females, high magn.") +
  geom_abline(intercept=796.695, slope=5.414, color="#E69F00", linetype=2) +      
  geom_abline(intercept=782.852, slope=1.535, color="deepskyblue3", linetype=2)   
print(p)
#
