#Analysis of sperm 


library(tidyverse)
library(MASS)
library(countreg)

calculateDispersionParameter <- function(mod, data){
  E1 <- resid(mod, type="pearson")
  return(sum(E1^2)/(nrow(data)/length(coef(mod))))
}

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

sperm <- sperm %>% mutate(Treatment2 = ifelse(Treatment %in% c("Male and female", "Male and Females"), "Mixed", ifelse(!is.na(Treatment), "Males only", NA) ), 
                          Time=factor(Hours))



#What is "ID"? Doesn't seem to be fly ID so I'm not sure what that is. Is it a
#pseudoreplicate? Should that be included as a random effect?

names(sperm) <- c("ID", "Population", "Treatment", "Hours",  "SpermBundles", "MatureSperm", "Treatment2", "Time")





###################################################################################
#Does the amount of mature sperm found in a fly vary dependings on the
#population it is from (only 2 options), or how long the fly has been alive?
ggplot(data=sperm, aes(x=Time, y=MatureSperm, fill=Population))+
  geom_violin(position=position_dodge(1))+
  stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))+
  facet_grid(Treatment2~.)

hist(sperm2$MatureSperm)








############################################################################################
#Does the amount of immature sperm found in a fly vary dependings on the
#population it is from (only 2 options), or how long the fly has been alive?
#Will do males ONLY

ggplot(data=sperm, aes(x=Time, y=SpermBundles, fill=Population))+
  geom_violin(position=position_dodge(1))+
  stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))+
  facet_grid(Treatment2~.)

hist(sperm2$ImatureSperm)
#looks pretty well zero-inflated.

malesperm <- sperm2 %>% filter(Treatment2=="Males only")

