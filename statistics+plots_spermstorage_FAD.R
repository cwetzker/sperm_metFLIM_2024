#############################################################################
## statistical analysis of FAD FLIM data #################################### 
## author: Cornelia Wetzker #################################################
## license: BSD-3-Clause license ############################################
#############################################################################
#
# packages used
library(lme4)
library(nlme)
library(car)
library(ggplot2)
library(Rmisc)
#
rm(list=ls()) 
#
### DATA IMPORT
flimF<-read.csv("E:/path_to_data/FAD_LT_data_table.csv")
str(flimF)
#
flimF$fluor<-as.factor(flimF$fluor)
flimF$sex<-as.factor(flimF$sex)
flimF$coh<-as.factor(flimF$coh)
levels(flimF$coh) <- c("mated females", "mated males", "unmated males")
flimF$mat<-as.factor(flimF$mat)
flimF$magn<-as.factor(flimF$magn)
flimF$lim_phot_number<-as.factor(flimF$lim_phot_number)
str(flimF)
#
#Exclusion of samples based on quality measures (data quality and photon counts)
#
# minimum of 1000 photons per region with binning 3
flimxF <- subset(flimF, ph_px_bin3 > 1000)    
#
# exclusion of images with CHI2 > 1.4
flimxF <- subset(flimxF, Chi2 < 1.4)
str(flimxF)
summary(flimxF)
#
flimxF_high <- flimxF
#
# Exclusion of unmated males, comparison of mated males and females
flimxF_high_mated <- subset(flimxF_high, coh != "unmated males")   # selection mated males + females
#
### STATISTICAL ANALYSES
#
# test for normal distribution of tm (mated only)
hist(flimxF_high_mated$FAD_tm)     
shapiro.test(flimxF_high_mated$FAD_tm) 
#
# test for normal distribution of a1% (mated only)
hist((flimxF_high_mated$FAD_a1))      
shapiro.test(flimxF_high_mated$FAD_a1)
#
### Test for homoscedasticity using Levene's test ####
#
mod_tm <- lmer(FAD_tm~sex*dpm+(1|sp_id), data=flimxF_high_mated)
leveneTest(residuals(mod_tm) ~ flimxF_high_mated$sex)
boxplot(residuals(mod_tm) ~ flimxF_high_mated$sex)
#
mod_a1 <- lmer(FAD_a1~sex*dpm+(1|sp_id), data=flimxF_high_mated)
leveneTest(residuals(mod_a1) ~ flimxF_high_mated$sex)
boxplot(residuals(mod_a1) ~ flimxF_high_mated$sex)
#
# Summary of statistics
#
summ_high_FAD_a1 <- summarySE(flimxF_high, measurevar = "FAD_a1", groupvars = c("coh", "dpm"))   # generates mean and SE for a1 dataset
summ_high_FAD_a1
summ_high_FAD_tm <- summarySE(flimxF_high, measurevar = "FAD_tm", groupvars = c("coh", "dpm"))   # generates mean and SE for a1 dataset
summ_high_FAD_tm
#
# Linear mixed model for FAD a1%
lmer_mated_a1_high <- with(flimxF_high_mated,lmer(FAD_a1~dpm*sex+(1|sp_id)))
Anova(lmer_mated_a1_high,type="III")
plot(lmer_mated_a1_high)
#
# Linear mixed model for FAD tm
lmer_mated_tm_high <- with(flimxF_high_mated,lmer(FAD_tm~dpm*sex+(1|sp_id)))
Anova(lmer_mated_tm_high,type="III")
plot(lmer_mated_tm_high)
#
### VISUALIZATION OF DATA 
#
# Plotting of a1% of mated males and females across days post mating
summ_high_mated <- summarySE(flimxF_high_mated, measurevar = "FAD_a1", groupvars = c("coh", "dpm"))   # for mean and SE for dataset
summ_high_mated
summary(flimxF_high_mated$ph_px_bin3)
flimxF_high_fma <- subset(flimxF_high_mated, coh == "mated females")  
lm_high_fma <- with(flimxF_high_fma,lm(FAD_a1~dpm))    
lm_high_fma    
flimxF_high_mma <- subset(flimxF_high_mated, coh == "mated males")  
lm_high_mma <- with(flimxF_high_mma,lm(FAD_a1~dpm))    
lm_high_mma    
p<- ggplot(summ_high_mated, aes(x=dpm, y=FAD_a1, group=coh, color=coh, shape=coh)) + 
  geom_point(size=3)+
  scale_color_manual(values=c("#E69F00", "deepskyblue3")) +       
  scale_shape_manual(values=c(15, 17))+
  geom_errorbar(aes(ymin=FAD_a1-se, ymax=FAD_a1+se), width=.1,
                position=position_dodge(0.05)) +
  scale_x_continuous(name="Days post mating") +                                     
  scale_y_continuous(name="% short FAD LT fraction") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15)) +
  labs(title="FAD a1%: Mated males and females, high magn.") +
  geom_abline(intercept=66.5466, slope=-0.5128, color="#E69F00", linetype=2) +          
  geom_abline(intercept=75.7283, slope=-0.2031, color="deepskyblue3", linetype=2)    
print(p)
#
# Plotting of tm of mated males and females across days post mating
#
summ_high_mated <- summarySE(flimxF_high_mated, measurevar = "FAD_tm", groupvars = c("coh", "dpm"))   # generates mean and SE for dataset
summ_high_mated
flimxF_high_fma <- subset(flimxF_high_mated, coh == "mated females")  
lm_high_fma <- with(flimxF_high_fma,lm(FAD_tm~dpm))    
lm_high_fma    
flimxF_high_mma <- subset(flimxF_high_mated, coh == "mated males")  
lm_high_mma <- with(flimxF_high_mma,lm(FAD_tm~dpm))    
lm_high_mma    
#
p<- ggplot(summ_high_mated, aes(x=dpm, y=FAD_tm, group=coh, color=coh, shape=coh)) + 
  geom_point(size=3)+
  scale_color_manual(values=c("#E69F00", "deepskyblue3")) +       
  scale_shape_manual(values=c(15, 17))+
  geom_errorbar(aes(ymin=FAD_tm-se, ymax=FAD_tm+se), width=.1,
                position=position_dodge(0.05)) +
  scale_x_continuous(name="Days post mating") +                                     
  scale_y_continuous(name="Mean FAD LT") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15)) +
  labs(title="Mean FAD LT: Mated males and females, high magn.") +
  geom_abline(intercept=865.501, slope=3.645, color="#E69F00", linetype=2) +      
  geom_abline(intercept=167.05, slope=17.71, color="deepskyblue3", linetype=2)   
print(p)
#
