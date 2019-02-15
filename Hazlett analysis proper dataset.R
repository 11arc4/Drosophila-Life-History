#Analysis of sperm 


library(tidyverse)
library(MASS)
library(countreg)


# The questions are related to life history tradeoffs and the cost of
# gametogenesis: 

#. Are males in populations selected to develop rapidly and
# breed early limited in fitness by the rate or amount spermatogenesis? Evidence
# would consist of low sperm numbers available within the first 2 days of adult
# life, reduced numbers of sperm bundles etc. 

#. What have males had to give up
# to achieve super-rapid generation time? [The answer is: Almost everything, as
# the size / GSI data will show. 

#. By extension: are females in these
# populations likely to encounter sperm-depleted males that should be avoided?
# Evidence would consist of males in the mixed-sex treatment often having very
# few sperm relative to those in the male-only treatment.
#
# To develop the story, we will be adding data on mating rate, body size, GSI,
# sperm length and sperm competition (P2) outcomes. 


#C population is long lifecyle
#A population is short lifecycle


sperm <- read.csv("~/Montogomerie Work/Drosophila Sperm/Hazlett Sperm Count Data.csv", na.strings = "", as.is=T)
levels(as.factor(sperm$Treatment))
names(sperm) <- c("ID", "Population", "Treatment", "Hours",  "SpermBundles", "MatureSperm")

sperm <- sperm %>% mutate(Treatment2 = ifelse(Treatment %in% c("Male and female", "Male and Females"), "Mixed", ifelse(!is.na(Treatment), "Males only", NA) ), 
                          Time=factor(Hours, levels=c(2,24,48)))



#What is "ID"? Doesn't seem to be fly ID so I'm not sure what that is. Is it a
#pseudoreplicate? Should that be included as a random effect?


summary(sperm)
sperm$Population <- factor(sperm$Population, levels=c("A", "C"))
sperm$Treatment2 <- factor(sperm$Treatment2, levels=c("Males only", "Mixed"))
malesperm <- sperm %>% filter(Treatment2=="Males only")
sperm2 <- sperm %>% filter(Time !=2)



###################################################################################
#Does the amount of mature sperm found in a fly vary dependings on the
#population it is from (only 2 options), or how long the fly has been alive?
ggplot(data=sperm, aes(x=Time, y=MatureSperm, fill=Population))+
  geom_boxplot()+
  #stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))+
  facet_grid(~Treatment2)

hist(sperm2$MatureSperm)



mod_mature <- glm(MatureSperm~Time*Population*Treatment2, data=sperm, family="poisson")

hist(sperm2$MatureSperm) # lot of tails. 

plot(mod_mature)
hist(resid(mod_mature)) #Looks OK
plot(resid(mod_mature)~sperm$Time) #Might want to account for increasing variance with time. 
plot(resid(mod_mature)~sperm$Population) #good
plot(resid(mod_mature)~sperm$Treatment2) #good
#Shoot that's huge. Maybe we should be using a quassi poisson or neg bin
AER::dispersiontest(mod_mature) #80.5
summary(mod_mature)


##Try a quassi poisson-- problematic because now we have to use quassi likelihood and it's just not comparable to pitnick analyses. 
# mod_mature_QP <- glm(MatureSperm~Time*Treatment2*Population, data=sperm2, family="quasipoisson")
# summary(mod_mature_QP)
# anova(mod_mature_QP)
# plot(mod_mature_QP)
# hist(resid(mod_mature)) 
# dredge(mod_mature_QP)
# car::Anova(mod_mature_QP)
# 
# emmeans(mod_mature_QP, list(pairwise ~ Time*Treatment2*Population), adjust = "tukey")
# 
# #R makes it hellish to get QAIC (which is what you need for a Quassi poisson)
# #DO not think this is owrking
# dfun <- function(object){
# with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
#   }
# 
# dredge(mod_mature_QP,rank="QAIC", chat=dfun(mod_mature))

#Try neg. binomial analysis--better because now we can use te normal likelihood based metrics. 
mod_mature_nb <- glm.nb(MatureSperm~Time*Treatment2*Population, data=sperm)
plot(mod_mature_nb) #Scale location plot slight trend down but I don't think you'd hardly notice it if it wasn't for the red line. 
hist(resid(mod_mature_nb)) #Looks OK
plot(resid(mod_mature_nb)~sperm$Time) #Might want to account for changing variance with time: probably a lot of this is due to larger sample size in 24 and 48 though....
plot(resid(mod_mature_nb)~sperm$Population) #good
plot(resid(mod_mature_nb)~sperm$Treatment2) #good

summary(mod_mature_nb)
AICc(mod_mature_nb, mod_mature) #Holyyyy this is way way way better. 

car::Anova(mod_mature_nb)
dredge(mod_mature_nb)
#appears we need to keep Population and the interaction between time and
#treatment but don't need the 3 way interaction or other two way interactions

mam_mature_nb <- glm.nb(MatureSperm ~ Time * Treatment2 + Population, data=sperm)
summary(mam_mature_nb)
emmeans(mod_mature_nb, list(pairwise ~ Time*Treatment2+Population), adjust = "tukey")


MSLables <- sperm %>% group_by(Treatment2, Time, Population ) %>% summarize(y =max(MatureSperm))

MSLables$Label <- c("A", "B", "CD", "EF", "FG", "FG", "H", "CH", "DE", "FG")



ggplot(data=sperm, aes(x=Time, y=MatureSperm, fill=Population))+
  #geom_violin(position=position_dodge(1))+
  geom_boxplot( alpha=0.8, position=position_dodge(1))+
  labs(x="Hours", y="Mature sperm", fill="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_text(data=MSLables,  aes(x=Time, group=Population, y=y+120, label = Label), position=position_dodge(1))+
  scale_fill_grey(labels=c("ACO", "CO"))+
  facet_grid(~Treatment2)

ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Mature sperm.jpeg", units="in", height=3, width=7, dev="jpeg")


############################################################################################
#Does the amount of immature sperm found in a fly vary dependings on the
#population it is from (only 2 options), or how long the fly has been alive?
#Will do males ONLY

ggplot(data=sperm, aes(x=Time, y=SpermBundles, fill=Population))+
  #geom_violin(position=position_dodge(1))+
  geom_boxplot()+
  #stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))+
  facet_grid(~Treatment2)

hist(sperm$SpermBundles)


mod_bundle <- glm(SpermBundles ~ Time*Population*Treatment2, data=sperm, family="poisson")
plot(mod_bundle) #looks pretty good
hist(resid(mod_bundle)) #Looks OK
plot(resid(mod_bundle)~sperm$Time) #good
plot(resid(mod_bundle)~sperm$Population) #good
plot(resid(mod_bundle)~sperm$Treatment2) #good
AER::dispersiontest(mod_bundle) #A bit overdispersed, but not nearly as bad as the mature sperm (5.7). Still worth trying a neg. bin though



mod_bundle_nb <- glm.nb(SpermBundles ~ Time*Population*Treatment2, data=sperm)
plot(mod_bundle_nb) #looks pretty good
hist(resid(mod_bundle_nb)) #Looks OK. One low outlier
plot(resid(mod_bundle_nb)~sperm$Time) #good
plot(resid(mod_bundle_nb)~sperm$Population) #good
plot(resid(mod_bundle_nb)~sperm$Treatment2) #good

#row 59 (male only, 24,  A, only 1 sperm bundle) is an outlier. 

AICc(mod_bundle, mod_bundle_nb)
#Yup that neg bin is much better. Accounts for overdispersion nicely. 

car::Anova(mod_bundle_nb)
anova(mod_bundle_nb)
dredge(mod_bundle_nb)
#looks like there are populaiton and time main effect and a maybe a hint of evidence there might be a time by population interaction
mam_bundle_nb <- glm.nb(SpermBundles ~ Time+Population, data=sperm)
mam_bundle_nb2<- glm.nb(SpermBundles ~ Population, data=sperm)
mam_bundle_nb3 <- glm.nb(SpermBundles ~ Time*Population, data=sperm)
mam_bundle_nb4 <- glm.nb(SpermBundles ~ Time, data=sperm)
mam_null <- glm.nb(SpermBundles ~ 1, data=sperm)

anova(mam_bundle_nb, mam_bundle_nb2, mam_bundle_nb3, mam_bundle_nb4 , mam_null)
AICc(mam_bundle_nb, mam_bundle_nb2, mam_bundle_nb3, mam_null)
car::Anova(mam_bundle_nb)
summary(mam_bundle_nb)
emmeans(mam_bundle_nb, list(pairwise ~ Time*Population), adjust = "tukey")

SBLables <- sperm %>% group_by(Time, Population) %>% summarize(y =max(SpermBundles))

SBLables$Label <- c("A", "B", "A", "B","A", "B")

ggplot(data=sperm, aes(x=Time, y=SpermBundles, fill=Population))+
  #geom_violin(position=position_dodge(1))+
  geom_boxplot( alpha=0.8, position=position_dodge(1))+
  labs(x="Hours", y="Sperm bundles", fill="Selection \nregime")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_text(data=SBLables,  aes(x=Time, group=Population, y=y+5, label = Label), position=position_dodge(1))+
  scale_fill_grey(labels=c("ACO", "CO"))


ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Sperm Bundles.jpeg", units="in", height=3, width=4, dev="jpeg")

