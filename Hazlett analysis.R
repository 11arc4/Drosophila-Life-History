#Analysis of sperm 


library(tidyverse)


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


sperm <- read.csv("~/Montogomerie Work/Drosophila Sperm/Hazlett Sperm Count Data Sheet 2.csv", na.strings = "", as.is=T)

sperm <- sperm %>% mutate(Treatment2 = ifelse(Treatment %in% c("Male and female", "Male and Females"), "Mixed", ifelse(!is.na(Treatment), "Males only", NA) ))

sperm$Time <- as.numeric(gsub( "\\D", "", sperm$Hours..Since.Eclosion.))            

sperm2 <- sperm %>% filter(!is.na(Time))
#There are about 500 NAs-- who are those? WHAT are those. I guess they are ones
#where we don't know the treatment or day old so I've excluded them. 

#What is "ID"? Doesn't seem to be fly ID so I'm not sure what that is. Is it a
#pseudoreplicate? Should that be included as a random effect?

names(sperm2) <- c("ID", "Population", "Treatment", "Hours", "Teste", "Section", "ImatureSperm", "MatureSperm", "Treatment2", "Time")





###################################################################################
#Does the amount of immature sperm found in a fly vary dependings on the
#population it is from (only 2 options), or how long the fly has been alive?
ggplot(data=sperm2, aes(x=factor(Time), y=MatureSperm, fill=Population))+
  geom_violin(position=position_dodge(1))+
  stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))+
  facet_grid(Treatment2~.)

############################################################################################
#Does the amount of mature sperm found in a fly vary dependings on the
#population it is from (only 2 options), or how long the fly has been alive?

ggplot(data=sperm2, aes(x=factor(Time), y=ImatureSperm, fill=Population))+
  geom_violin(position=position_dodge(1))+
  stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))+
  facet_grid(Treatment2~.)




hist(sperm2$MatureSperm)
hist(sperm2$ImatureSperm)

#Both look pretty well zero-inflated. 



mod <- glm(ImatureSperm ~ Time*Population*Treatment2, data=sperm2, family="poisson")
plot(mod)
summary(mod)
#calculate dispersion parameter. 

#GOF test shows that we are NOT fitting the data well- probably because it's zero inflated. 
1 - pchisq(summary(mod)$deviance, 
           summary(mod)$df.residual
) #should be >0.05 if we fit well


#Lets try a negative binomial from MASS
library(MASS)
mod.nb <- glm.nb(ImatureSperm ~ Time*Population*Treatment2, data=sperm2)
plot(mod.nb)
summary(mod.nb)


1 - pchisq(summary(mod.nb)$deviance, 
           summary(mod.nb)$df.residual
)
#still no good. Probably still because we aren't accounting for the zero inflation


#Lets try zero inflated from pscl
library(countreg)
#1 means that the probability of a 0 is the same in all groups-- probably won't fit
mod.zip <- zeroinfl(ImatureSperm ~ Time*Population*Treatment2|1, data = sperm2)
summary(mod.zip)
rootogram(mod.zip, main = "ZIP", ylim = c(-5, 15), max = 50) #ooo year really not fitting that well. We are still super overdispersed. 
qqrplot(mod.zip, main = "ZIP")
#We aren't fitting that super well yet. 




mod.zip2 <- zeroinfl(ImatureSperm ~ Time*Population*Treatment2|Time*Population*Treatment2, data = sperm2)
rootogram(mod.zip2, main = "ZIP", ylim = c(-5, 15), max = 50) #ooo year really not fitting that well. We are still super overdispersed. 
qqrplot(mod.zip2, main = "ZIP")
#This is not any better. Changes almost NOTHING


