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

ggplot(data=malesperm, aes(x=Time, y=ImatureSperm, fill=Population))+
  geom_violin(position=position_dodge(1))+
  stat_summary(fun.y=median, geom="point", size=2, aes(group=Population), position=position_dodge(1))


ggplot(data=malesperm, aes(x=Time, y=ImatureSperm, fill=Population))+
  geom_boxplot(position=position_dodge(1))



nd <- malesperm %>% group_by(Population, Time) %>% summarise(median(ImatureSperm),
                                                             mean(ImatureSperm)
                                                             )



mod <- glm(ImatureSperm ~ Time*Population, data=malesperm, family="poisson")
plot(mod)
summary(mod)
#calculate dispersion parameter. 
calculateDispersionParameter(mod,malesperm)
#37.0 is pretty over dispered. Not good at all. 
pois<- cbind(nd, 
           Mean = predict(mod, newdata=nd, type="response"), 
           SE = predict(mod, newdata=nd, type="response", se.fit=T)$se.fit
)


#Lets try a negative binomial from MASS
mod.nb <- glm.nb(ImatureSperm ~ Time*Population, data=malesperm)
plot(mod.nb)
summary(mod.nb)
rootogram(mod.nb) #That's not horrible. It's not great either but it's not awful
calculateDispersionParameter(mod.nb, malesperm)
#much better but still overdispersed (dp=4.6)

1 - pchisq(summary(mod.nb)$deviance,
           summary(mod.nb)$df.residual
)
#Doesn't fit properly

#numbers match almost too well? Appears that their mean matches exactly and all changing type of model does is change SE-- makes sense. 
nb<- cbind(nd, 
      Mean = predict(mod.nb, newdata=nd, type="response"), 
      SE = predict(mod.nb, newdata=nd, type="response", se.fit=T)$se.fit
)


#Lets try zero inflated from pscl

mod.zip <- zeroinfl(ImatureSperm ~ Time*Population|Time*Population, data = malesperm)
rootogram(mod.zip, main = "ZIP") #ooo year really not fitting that well. We are still super overdispersed. 
qqrplot(mod.zip, main = "ZIP")
#This is not any better. Changes almost NOTHING



#numbers don't match quite right.....
cbind(nd, 
      Mean = predict(mod.zip, newdata=nd, type="response")
      #SE = predict(mod.zip, newdata=nd, type="response", se.fit=T)$se.fit
)
#THIS IS MUCH WORSE

#What about a zero inflated neg bin?

mod.zinb = zeroinfl(ImatureSperm ~ Time*Population|Time*Population, data = malesperm, dist = "negbin")
summary(mod.zinb) #theta is significant suggesting that the ZIP model doesn't fit because it's over dispersed. 
rootogram(mod.zinb) #ooo this looks very very nice. 
qqrplot(mod.zip2) #I'm not sure this applies anymore
#Disperson Parameter
mod.zinb$theta
zinb <- cbind(nd, 
      Mean = predict(mod.zinb, newdata=nd, type="response")
      #SE = predict(mod.zinb, newdata=nd, type="response", se.fit=T)$se.fit
)


dredge(mod.zinb)
AICc(mod.nb, mod.zinb) #I think this means that our zinb is the best model. 

#Comparing estimates to what's real, looks very much like the neg bin model is
#best, even though it's not perfect. That's funny because the rootogram is much
#better for the zero inflated neg bin. maybe I should plot both?

mam.zinb <- zeroinfl(ImatureSperm ~ Time*Population|Time*Population, data = malesperm, dist = "negbin")
summary(mam.zinb)
#CO has more immature sperm but sperm from all populations does decrease over time

zinbresult <- emmeans(mod.zinb, list(pairwise ~ Population+Time), adjust = "tukey")

zinbresultLabels <- as.data.frame(zinbresult$`emmeans of Population, Time`)
zinbresultLabels$Labels <- c("A", "B", "C", "A", "C", "AC")

ggplot()+
  geom_violin(data=malesperm, aes(x=Time, y=ImatureSperm, fill=Population), alpha=0.3,  position=position_dodge(1))+
  labs(x="Time (hrs)", y="Immature sperm")+
  theme_classic(base_family = "serif", base_size = 14)+
  geom_text(data=zinbresultLabels,  aes(x=Time, y=30, label = Labels, group=Population), position=position_dodge(1))+
  geom_point(data=zinbresultLabels, aes(x=Time, y=emmean, group=Population), position=position_dodge(1))+
  geom_linerange(data=zinbresultLabels, aes(x=Time, ymin=emmean+SE, ymax=emmean-SE, group=Population), position=position_dodge(1))+
  scale_fill_grey()
ggsave("~/Montogomerie Work/Drosophila Sperm/Plots/Immature Sperm Count_males only.jpeg", units="in", height=3, width=5, dev="jpeg")




