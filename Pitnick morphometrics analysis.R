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

# #could winsorize testes mass to get rid of high and low values
# #better to winsorize by group-- once you winsorize by group you've no longer caused major issues. 
# morph <- morph %>% group_by(Selection) %>% mutate(WinGSoma= winsor(Soma, trim=0.03), 
#                                                   WinGTestes = winsor(Testes, trim=0.03))
# #quantile (morph$Testes, c(0.05, 0.95))
# hist(morph$WinGTestes)
# summary(morph$WinGTestes)
#this rather messes up the distribution

morph$Testes2 [morph$Selection=="CO" & morph$Testes2 >0.009] <- NA #TDiff is high-- probably should be removed

testes <- morph %>% filter(!is.na(Testes2))

ggplot(testes, aes(x=Selection, y=Testes2, fill=Selection))+
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


ggplot(gonad, aes(x=Testes2, y=Soma2, color=Selection))+
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




summary(mod_residgonad) #don't really need the random effect. 
dredge(mod_residgonad)
anova(mod_residgonad)
#There are in fact real differences between relative gonad weight in the different populations

mam_residgonad <- lmer(Testes2 ~Soma2+Selection + (1|Line), data=gonad, REML=F)

#no indication the interaction term does anything. However, I THINK (?) this
#data is showing that CO puts more investment into the testes relative to the
#Soma, regardless of testes or soma size, than CO or IV


ggplot(gonad, aes(x=Testes2, y=Soma2, color=Selection, shape=Selection))+
  geom_point(size=2)+
  geom_smooth( method="lm", se=F)+
  labs(x="Soma mass (g)", y="Testes mass (g)", color="Selection \nregime", shape="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  scale_color_grey()



######################################################################################################################################
#Has the relative investment in reproductive organs to non-reproductive organs changed based on selection regimes?

gonad$Reproductive <- gonad$Soma2 + gonad$Testes2

ggplot(gonad, aes(x=Thorax, y=Reproductive, color=Selection))+
  geom_point()+
  labs(x="Thorax", y="Reproductive tissue")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_smooth(method="lm")






mod_rep <- lmer(Reproductive ~  Thorax * Selection + (1|Line), data=gonad, REML=F)
plot(mod_rep)
hist(resid(mod_rep)) #Looking very good
shapiro.test(resid(mod_rep))
plot(resid(mod_rep)~factor(gonad$Selection)) #looks great
plot(resid(mod_rep)~factor(gonad$Line)) #looks all good
plot(resid(mod_rep)~gonad$Thorax) #variance might increase a bit with thorax-- might consider accounting for that. 




summary(mod_rep) #don't really need the random effect. 
dredge(mod_rep)
anova(mod_rep)
# no interaction, but both main effects are real


mam_rep <- lmer(Reproductive ~  Thorax + Selection + (1|Line), data=gonad, REML=F)








