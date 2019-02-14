#How has drosophila morphology changed based on the two different selection regimes? 

library(tidyverse)
library(lme4)
library(MuMIn)
library(emmeans)
library(psych)
library(lmerTest)

options(na.action = "na.fail")

# The IV population is included in some of the datasets, in addition to the A-
# and C-selected populations. IV is the base population from which all selected
# populations were ultimately derived. It is assumed to have been at (or close
# to) equilibrium when the new populations were founded because of a long
# history of maintenance on a 14d cycle. Note that IV is intermediate to A and C
# in schedule and also often in phenotype as well, suggesting that selection has
# operated in both directions (e.g., towards larger and smaller size).
#
# For the morphology (Pitnick) data, our aim will be to characterize the
# differences between A- and C-selected males for their length, body mass,
# testis mass, and GSI and assess the directionality of response using the IV
# data point.
#
# For each character measured, we will want to see if there is a difference
# between the A and C selection treatments, consisting of populations A1-5 vs
# C1-5 with replicate number considered a factor due to shared history of
# same-numbered populations in each treatment (e.g. A1 and C1).
#
# Measurement error has produced some outlying values for which we discussed
# windsorizing to clip extreme high and low weights. I think the GSI (as
# calculated on Scott’s spreadsheet) probably captures all of the wierdnesses
# and could be used as the basis for removing outliers. I assume we will remove
# the individual altogether if anything is squirrelly (i.e.,
#
# Notes on characters to analyze:
#
# 1) Thorax length should be a straightforward two-factor model (selection and
# replicate population number). There is no information in the ID# of the male.
#
# 2) “Condition” was calculated in an unfamiliar and non-intuitive way (length /
# soma mass, which seems backwards) and it may be more appropriate to use a
# condition index based upon residuals, such as the SMI employed by Dakin and
# Mont 2014. Even then, is the best index of condition based upon soma weight
# alone or total weight (testes included?). Note too that these data tend to
# constitute two clusters (A = small and light; C = big and heavy) and some of
# the overall regression slope is dictated by the pull of the A-blob and C-blob,
# while the relationship within each blob may differ.
#
# 3) Soma weight like thorax should be quite simple to analyze.
#
# 4) Testis weight like thorax should be quite simple to analyze.
#
# 5) GSI was calculated as gonad / (gonad + body) in the current spreadsheet
# calculations but could be approached in the same way condition is.

#Read in Pitnick's body size and morphology data
morph <- read.csv("file:///C:/Users/11arc/OneDrive/Documents/Montogomerie Work/Drosophila Sperm/PitnickData_ACO-CO.GSI..csv", as.is = T)
#use "Line" as a unique ID for each population (population is not unique)


morph$Selection <- factor(morph$Selection, levels=c("ACO", "IV", "CO"))
morph$Line <- factor(morph$Line)
morph$Testes[ morph$Testes< (-0.015)] <- NA #Remove that extreme outlier-- clearly not right and unclear what it should be. 

#Decided to 
morph$Testes2 <- ifelse(morph$Testes>0, morph$Testes, min(morph$Testes[morph$Testes>0], na.rm=T))

which(morph$Testes<0)
##############################################################################################
#Has thorax length changed depending on the selection regime?
names(morph)[names(morph)=="Thorax..mm."] <- "Thorax"
hist(morph$Thorax) #looks pretty OK-- no weight outliers. 

ggplot(morph, aes(x=Line, y=Thorax, fill=Selection))+
  geom_boxplot(show.legend = F)+
  geom_jitter(width = 0.1)+
  labs(x="Selection regime", y="Thorax length (mm)")+
  theme_classic(base_family = "serif", base_size = 14)
#clear trend. This really is super easy. 



mod_thorax <- lmer(Thorax ~ Selection + (1|Line), data=morph, REML=F)
plot(mod_thorax)
hist(resid(mod_thorax)) #looks pretty normal
plot(resid(mod_thorax)~factor(morph$Selection))
plot(resid(mod_thorax)~factor(morph$Line))
#this looks great. Fits really good


summary(mod_thorax) #don't really need the random effect. 
dredge(mod_thorax)
anova(mod_thorax)
#There are in fact real differences between thorax length in the different populations

emmeans(mod_thorax, list(pairwise ~ Selection), adjust = "tukey")
#All groups differ


ThoraxLables <- morph %>% group_by(Selection) %>% summarize(y =max(Thorax))

ThoraxLables$Label <- c("A", "B", "C")

ggplot(morph, aes(x=Selection, y=Thorax, fill=Selection))+
  geom_boxplot(show.legend = F, alpha=0.8)+
  labs(x="Selection regime", y="Thorax length (mm)")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_text(data=ThoraxLables,  aes(x=Selection, y=y+0.01, label = Label))+
  scale_fill_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Thorax length.jpeg", units="in", height=3, width=4, dev="jpeg")


###############################################################################################################
#Has testes weight changed depending on the selection regime?

hist(morph$Testes)

morph$Testes2 [morph$Selection=="CO" & morph$Testes2 >0.009] <- NA #TDiff is high-- probably should be removed

testes <- morph %>% filter(!is.na(Testes2))

ggplot(testes, aes(x=Selection, y=Testes2, fill=Selection, group=Line))+
  geom_boxplot(show.legend = F)+
  labs(x="Selection regime", y="Testes weight (g)")+
  theme_classic(base_family = "serif", base_size = 14)


ggplot(morph, aes(x=Testes2)) +
  geom_histogram(binwidth = 0.001)+
  facet_grid(Selection~.)


#"Outlier" in IV is not an outlier, just larger than the others (TDiff is normal)



mod_testes <- lmer(Testes2 ~ Selection + (1|Line), data=testes, REML=F)
plot(mod_testes)
hist(resid(mod_testes)) #Doesn't a bit left skewed....
shapiro.test(resid(mod_testes))
plot(resid(mod_testes)~factor(testes$Selection)) # population IV has slightly less variation but overall pretty OK
plot(resid(mod_testes)~factor(testes$Line)) #looks all good
#this looks OK. Fits decently, but could be better (residuals slightly left skewed) 


summary(mod_testes) #don't really need the random effect. 
dredge(mod_testes)
anova(mod_testes)
#There are in fact real differences between stoma mass in the different populations


emmeans(mod_testes, list(pairwise ~ Selection), adjust = "tukey")
# ACO is lower than CO and IV-- those two are the same. 


TestesLables <- testes %>% group_by(Selection) %>% summarize(y =max(Testes2))

TestesLables$Label <- c("A", "B", "B")


ggplot(testes, aes(x=Selection, y=Testes2, fill=Selection))+
  geom_boxplot(show.legend = F, alpha=0.8)+
  labs(x="Selection regime", y="Testes mass (g)")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_text(data=TestesLables,  aes(x=Selection, y=y+0.001, label = Label))+
  scale_fill_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Testes Mass.jpeg", units="in", height=3, width=4, dev="jpeg")

###############################################################################################################################
#Has soma weight changed depending on the selection regime?

morph$Soma2 <- morph$Soma
morph$Soma2[morph$Soma>0.25 & morph$Selection=="ACO"] <- NA #SDiff very high-- remove this point


soma <- morph %>% filter(!is.na(Soma2))


ggplot(soma, aes(x=Selection, y=Soma, fill=Selection))+
  geom_boxplot(show.legend = F)+
  labs(x="Selection regime", y="Soma weight (g)")+
  theme_classic(base_family = "serif", base_size = 14)
#clear trend. 

ggplot(morph, aes(x=Soma)) +
  geom_histogram()+
  facet_grid(Selection~.)


mod_soma <- lmer(Soma ~ Selection + (1|Line), data=soma, REML=F)
plot(mod_soma)
hist(resid(mod_soma)) #looks great
shapiro.test(resid(mod_soma))
plot(resid(mod_soma)~factor(soma$Selection)) # population IV has slightly less variation but overall pretty OK
plot(resid(mod_soma)~factor(soma$Line)) #looks all good


summary(mod_soma) #don't really need the random effect. 
dredge(mod_soma)
anova(mod_soma)
#There are in fact real differences between stoma mass in the different populations


emmeans(mod_soma, list(pairwise ~ Selection), adjust = "tukey")
#all selection regimes differ


SomaLabels <- soma %>% group_by(Selection) %>% summarize(y =max(Soma2))

SomaLabels$Label <- c("A", "B", "C")


ggplot(soma, aes(x=Selection, y=Soma2, fill=Selection))+
  geom_boxplot(show.legend = F, alpha=0.8)+
  labs(x="Selection regime", y="Soma mass (g)")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_text(data=SomaLabels,  aes(x=Selection, y=y+0.01, label = Label))+
  scale_fill_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Soma Mass.jpeg", units="in", height=3, width=4, dev="jpeg")


##################################################################################################################################
#Has relative investment in testes to soma changed depending on selection regime?

gonad <- morph %>% filter(!is.na(Testes2) & !is.na(Soma2))


ggplot(gonad, aes(y=Testes2, x=Soma2, color=Selection))+
  geom_point()+
  geom_smooth( method="lm")+
  labs(x="Soma", y="Testes")+
  theme_classic(base_family = "serif", base_size = 14)



mod_residgonad <- lmer(Testes2 ~ Soma2*Selection + (1|Line), data=gonad, REML=F) 
# could maybe consider adding a random slope for each lineage but I suspect that's overkill
plot(mod_residgonad)
hist(resid(mod_residgonad)) #Looking mediocre, little of a tail
shapiro.test(resid(mod_residgonad))
plot(resid(mod_residgonad)~factor(gonad$Selection)) #looks great
plot(resid(mod_residgonad)~factor(gonad$Line)) #looks all good
plot(resid(mod_residgonad)~gonad$Soma2) #looks all good
#overall looks pretty OK

emmeans(mod_residgonad)


summary(mod_residgonad) #don't really need the random effect. 
dredge(mod_residgonad)
anova(mod_residgonad)
#There are in fact real differences between relative gonad weight in the different populations

lm_residgonad <- lm(Testes2 ~Soma2+Selection , data=gonad, REML=F)


mod_residgonad <- lmer(Testes2 ~ Soma2+Selection + (1|Line), data=gonad, REML=F) 
summary(mod_residgonad)

#no indication the interaction term does anything. However, I THINK (?) this
#data is showing that CO puts more investment into the testes relative to the
#Soma, regardless of testes or soma size, than CO or IV


newdata<- data.frame(Selection=factor(c(rep("ACO", 20), rep("IV", 20), rep("CO", 20))), 
                     Soma2=c(seq(min(gonad$Soma2[gonad$Selection=="ACO"]), max(gonad$Soma2[gonad$Selection=="ACO"]), length.out=20),
                             seq(min(gonad$Soma2[gonad$Selection=="IV"]), max(gonad$Soma2[gonad$Selection=="IV"]), length.out=20),
                             seq(min(gonad$Soma2[gonad$Selection=="CO"]), max(gonad$Soma2[gonad$Selection=="CO"]), length.out=20)), 
                     Predicted=NA, 
                     lcl=NA, 
                     ucl=NA)

newdata$Predicted <- predict(mam_residgonad, newdata,re.form=~0) #level=0 tells us to ignore the random effect and just pick the mean nest!

ggplot()+
  geom_point(data=gonad, aes(x=Soma2, y=Testes2, color=Selection, shape=Selection), size=2)+
  geom_line(data=newdata, aes(x=Soma2, y=Predicted, color=Selection))+
  labs(x="Soma mass (g)", y="Testes mass (g)", color="Selection \nregime", shape="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  scale_color_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Testes mass by Soma Mass.jpeg", units="in", height=3, width=6, dev="jpeg")




######################################################################################################################################
#Has the relative investment in reproductive organs to non-reproductive organs changed based on selection regimes?

gonad$Mass <- gonad$Soma2 + gonad$Testes2

ggplot(gonad, aes(x=Thorax, y=Mass, color=Selection))+
  geom_point()+
  labs(x="Thorax", y="Mass")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_smooth(method="lm")






mod_mass <- lmer(Mass ~  Selection* Thorax  + (1|Line), data=gonad, REML=F)
plot(mod_mass)
hist(resid(mod_mass)) #Looking very good
shapiro.test(resid(mod_mass))
plot(resid(mod_mass)~factor(gonad$Selection)) #looks great
plot(resid(mod_mass)~factor(gonad$Line)) #looks all good
plot(resid(mod_mass)~gonad$Thorax) #variance might increase a bit with thorax-- might consider accounting for that. 




summary(mod_mass) #don't really need the random effect. 
dredge(mod_mass)
anova(mod_mass)
# no interaction, but both main effects are real


mam_mass <- lmer(Mass ~  Thorax + Selection + (1|Line), data=gonad, REML=F)

anova(mam_mass)




newdata<- data.frame(Selection=factor(c(rep("ACO", 20), rep("IV", 20), rep("CO", 20))), 
                     Thorax=c(seq(min(gonad$Thorax[gonad$Selection=="ACO"]), max(gonad$Thorax[gonad$Selection=="ACO"]), length.out=20),
                             seq(min(gonad$Thorax[gonad$Selection=="IV"]), max(gonad$Thorax[gonad$Selection=="IV"]), length.out=20),
                             seq(min(gonad$Thorax[gonad$Selection=="CO"]), max(gonad$Thorax[gonad$Selection=="CO"]), length.out=20)), 
                     Predicted=NA, 
                     lcl=NA, 
                     ucl=NA)

newdata$Predicted <- predict(mam_mass, newdata,re.form=~0) #level=0 tells us to ignore the random effect and just pick the mean nest!


ggplot()+
  geom_point(data=gonad, aes(x=Thorax, y=Mass, color=Selection, shape=Selection), size=2)+
  geom_line(data=newdata, aes(x=Thorax, y=Predicted, color=Selection))+
  labs(x="Thorax length (mm))", y="Body mass (g)", color="Selection \nregime", shape="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  scale_color_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Body mass by thorax length.jpeg", units="in", height=3, width=6, dev="jpeg")





############################################################################################################
#Calculate scaled mass index based on soma weight and thorax. 

ggplot(gonad, aes(x=log(Thorax), y=log(Soma2), color=Selection))+
  geom_point()+
  geom_abline(slope=3.7, intercept = -1.1)


#Calculate SMI either based off of entire population, or only ACO selected population
##For the entire set of flies. 
smiMod <- lmodel2::lmodel2(log(Soma2) ~ log(Thorax), data=gonad)
BSMA <- smiMod$regression.results[3,]

#Based off of only ACO
smiMod_ACO <- lmodel2::lmodel2(log(Soma2) ~ log(Thorax), data=gonad %>% filter(Selection=="ACO"))
BSMA_Aco <- smiMod_ACO$regression.results[3,]

#Calculate SMI-- set L0 as the mean body mass of group IV
L0 <- mean(gonad$Thorax[gonad$Selection=="IV"])
gonad$SMI <- gonad$Soma2 * ((L0/gonad$Thorax)^BSMA[[3]])
gonad$SMI_ACO <- gonad$Soma2 * ((L0/gonad$Thorax)^BSMA_Aco[[3]])

################################################################################################################
#Does selection regime affect SMI?

mod_condition <- lmer(SMI_ACO ~ Selection + (1|Line), data=gonad, REML=F)
plot(mod_condition)
hist(resid(mod_condition)) #very good
shapiro.test(resid(mod_condition))
plot(resid(mod_condition)~factor(gonad$Selection)) #looks great
plot(resid(mod_condition)~factor(gonad$Line)) #looks all good
#Resids not quite normal, but not horrible. Probably OK


summary(mod_condition) #don't really need the random effect. 
dredge(mod_condition)
anova(mod_condition)
#SMI not influenced by selection regime


GonadLabels <- gonad %>% group_by(Selection) %>% summarize(y =max(SMI_ACO))
GonadLabels$Label <- c("A", "B", "C")


ggplot(gonad, aes(x=Selection, y=SMI_ACO, fill=Selection))+
  geom_boxplot(show.legend = F, alpha=0.8)+
  labs(x="Selection regime", y="Scaled mass index")+
  theme_classic(base_family = "serif", base_size = 14)+
  #geom_text(data=SomaLabels,  aes(x=Selection, y=y+0.01, label = Label))+
  scale_fill_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Soma Mass.jpeg", units="in", height=3, width=4, dev="jpeg")



#################################################################################################################
#############Does selection regime affect the relationship between SMI and testes weight? 
#SMI based off of everyone
ggplot(gonad, aes(x=SMI, y=Testes2, color=Selection))+
  geom_point()+
  geom_smooth(method="lm")

mod_SMI <- lmer(Testes2 ~ SMI * Selection + (1|Line), data=gonad, REML=F)

plot(mod_SMI)
hist(resid(mod_SMI)) #not great
shapiro.test(resid(mod_SMI))
plot(resid(mod_SMI)~factor(gonad$Selection)) #looks great
plot(resid(mod_SMI)~factor(gonad$Line)) #looks all good
plot(resid(mod_SMI)~gonad$SMI) #Looks fantastic
#Resids not quite normal, but not horrible. Probably OK


summary(mod_SMI) #don't really need the random effect. 
dredge(mod_SMI)
anova(mod_SMI)
# no interaction and SMI doesn't really look like it does anything much. 


ggplot(gonad, aes(x=SMI, y=Testes2, color=Selection, shape=Selection))+
  geom_point(size=2)+
  #geom_smooth(method="lm", se=F)+
  labs(x="Scaled Mass Index", y="Testes mass (g)", color="Selection \nregime", shape="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  scale_color_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Testes Mass by SMI.jpeg", units="in", height=3, width=6, dev="jpeg")

#SMI based off of ACO only

ggplot(gonad, aes(x=SMI_ACO, y=Testes2, color=Selection))+
  geom_point()+
  geom_smooth(method="lm")

mod_SMI_ACO <- lmer(Testes2 ~ SMI_ACO * Selection + (1|Line), data=gonad, REML=F)

plot(mod_SMI_ACO)
hist(resid(mod_SMI_ACO)) #fine
shapiro.test(resid(mod_SMI_ACO))
plot(resid(mod_SMI_ACO)~factor(gonad$Selection)) #looks great
plot(resid(mod_SMI_ACO)~factor(gonad$Line)) #looks all good
plot(resid(mod_SMI_ACO)~gonad$SMI_ACO) #Looks fantastic
#Looks oK

summary(mod_SMI_ACO) #don't really need the random effect. 
dredge(mod_SMI_ACO)
anova(mod_SMI_ACO)
#Looks like we probably don't need the interaction but thats it


mam_SMI_ACO <- lmer(Testes2 ~ SMI_ACO + Selection + (1|Line), data=gonad, REML=F)

newdata<- data.frame(Selection=factor(c(rep("ACO", 20), rep("IV", 20), rep("CO", 20))), 
                     SMI_ACO=c(seq(min(gonad$SMI_ACO[gonad$Selection=="ACO"]), max(gonad$SMI_ACO[gonad$Selection=="ACO"]), length.out=20),
                              seq(min(gonad$SMI_ACO[gonad$Selection=="IV"]), max(gonad$SMI_ACO[gonad$Selection=="IV"]), length.out=20),
                              seq(min(gonad$SMI_ACO[gonad$Selection=="CO"]), max(gonad$SMI_ACO[gonad$Selection=="CO"]), length.out=20)), 
                     Predicted=NA, 
                     lcl=NA, 
                     ucl=NA)

newdata$Predicted <- predict(mam_SMI_ACO, newdata,re.form=~0) #level=0 tells us to ignore the random effect and just pick the mean nest!



ggplot(gonad, aes(x=SMI_ACO, y=Testes2, color=Selection, shape=Selection))+
  geom_point(size=2)+
  geom_line(data=newdata, aes(x=SMI_ACO, y=Predicted, color=Selection))+
  labs(x="Scaled Mass Index", y="Testes mass (g)", color="Selection \nregime", shape="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  scale_color_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Testes Mass by SMI (from ACO pop).jpeg", units="in", height=3, width=6, dev="jpeg")

