## ----setup, include=FALSE-------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----eval=TRUE, message=FALSE, results='hide',warning=FALSE---------------------------------------------
#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


packages <- c("unmarked", "reshape2", "dplyr", "ggplot2","AICcmodavg")

#run function to install packages
ipak(packages)


## -------------------------------------------------------------------------------------------------------
pcru <- csvToUMF(system.file("csv", "frog2001pcru.csv", package = "unmarked"),long = TRUE, type = "unmarkedFrameOccu")
head(pcru)
summary(pcru)


## -------------------------------------------------------------------------------------------------------
# Scale
obsCovs(pcru) <- scale(obsCovs(pcru))


## -------------------------------------------------------------------------------------------------------
# model 1: p(.)psi(.)
ex_mod_p1 <- occu(~ 1 ~ 1, data = pcru)
ex_mod_p1

# model 2: p(MinAfterSunset,Temperature)psi(.)
ex_mod_p2 <- occu(~ MinAfterSunset + Temperature ~ 1, data = pcru)
ex_mod_p2


## -------------------------------------------------------------------------------------------------------
backTransform(ex_mod_p1, "state")
backTransform(ex_mod_p1, "det")


backTransform(ex_mod_p2, "state")


## -------------------------------------------------------------------------------------------------------
backTransform(linearComb(ex_mod_p2, coefficients = c(0, 0, -1), type = "det"))


## -------------------------------------------------------------------------------------------------------
head(getData(ex_mod_p2))
newData <- data.frame(MinAfterSunset = 0, Temperature = -2:2)
head(newData)
predict(ex_mod_p2, type = "det", newdata = newData, appendData = TRUE)
predicted.data <- predict(ex_mod_p2, type = "det", newdata = newData, appendData = TRUE)


## -------------------------------------------------------------------------------------------------------
ggplot(predicted.data, aes(Temperature, Predicted)) + stat_smooth(method="glm", method.args = list(family="binomial"), level = 0.5)


## -------------------------------------------------------------------------------------------------------
fit.mlist <-fitList(fits=list('p.psi' = ex_mod_p1,'pMinAfterSunset_Temp_psi' = ex_mod_p2))


## -------------------------------------------------------------------------------------------------------
modSel(fit.mlist, nullmod = "p.psi")


## -------------------------------------------------------------------------------------------------------
coef(fit.mlist)
SE(fit.mlist)
frogPsi_m2 <-predict(fit.mlist, type="state") # occupancy
mean(frogPsi_m2$Predicted)


## -------------------------------------------------------------------------------------------------------
Naive = 96
Psi = 130*mean(frogPsi_m2$Predicted)
Psi - Naive


## -------------------------------------------------------------------------------------------------------
frogP_m2 <-predict(fit.mlist, type="det") # detection probability
hist(frogP_m2$Predicted)
summary(frogP_m2)


fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}


## ----warning=FALSE--------------------------------------------------------------------------------------
load("Data/elkforheb.RData") ##loads functions chisq and Nocc used in parboot statistics

pcru.pb <- parboot(ex_mod_p2, statistics = chisq, nsim = 50)
pcru.pb
plot(pcru.pb)


pcru.pb <- parboot(ex_mod_p2, Nocc, nsim = 50)
pcru.pb
plot(pcru.pb)
abline(v=96, col="red")



## -------------------------------------------------------------------------------------------------------
#load("Data/elkforheb.RData") 
ls()
str(ydata)
str(sp)
str(covar)

ggplot(covar, aes(easting, northing)) + geom_point()


## -------------------------------------------------------------------------------------------------------
new_covar <- covar %>%
  mutate(northness = cos(aspect90m)) %>%
  mutate(eastness = sin(aspect90m)) %>%
  rename(elev =dem90m) %>%
  mutate(elev2 = elev*elev) %>%
  mutate_if(is.integer, as.character) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_if(is.numeric, scale)

sel_cov <- new_covar %>%
  dplyr::select('northness', 'eastness', 'slope' = 'slope90m', 'elev',
                'elev2', 'd2road', 'tpi20', 'tpi100', 'tpi500', 'cc500', 'cc100', 'burns20','burns100', 'burns500', 'regen20', 'regen100', 'regen500', 'cuts20','cuts100', 'cuts500', 'NDVIAug500', 'NDVIAug20', 'NDVIJul500', 'NDVIJul100','NDVIJul20', 'NDVIAug100', 'dhicum20', 'dhicum100', 'dhicum500', 'dhimin20','dhimin100', 'dhimin500', 'dhiseas20', 'dhiseas100', 'dhiseas500', 'lure','protected2', 'protected3', 'trailtype',	'camera','ppltot', 'pplcat3a', 'pplcat3b', 'pplcat3c', 'pplcat5','motortot', 'motorcat3a', 'motorcat3b', 'motorcat3c', 'motorcat5','allhumantot', 'allhumancat3a', 'allhumancat3b', 'allhumancat3c','allhumancat5', 'pplcat2a', 'pplcat2b', 'pplcat2c', 'pplcat2d', 'motorcat2a','motorcat2b', 'motorcat2c', 'motorcat2d', 'allhumancat2a', 'allhumancat2b','allhumancat2c', 'allhumancat2d')

elk_cov <- sel_cov %>%
  dplyr::select('elev', 'slope', 'cuts100', 'pplcat3a',  'lure', 'camera') %>%
  mutate_if(is.numeric, scale)


## -------------------------------------------------------------------------------------------------------
#str(ydata) ## Note there are both NAs and 0s
ydata7 <- convert(ydata, 7)
ydata7 <- count2occu(ydata7) # makes sure occupancy data, not count

# Generate data in unmarked foramt
umf <- unmarkedFrameOccu(y = ydata7[,-1], siteCovs = sel_cov) 
summary(umf)



summary(umf)


## -------------------------------------------------------------------------------------------------------
fm_p1 <- occu(~1 ~1, data=umf)
#summary(fm_p1)
fm_p2 <- occu(~d2road ~1, data=umf)
fm_p3 <- occu(~tpi20 ~1, data=umf)
fm_p4 <- occu(~tpi100 ~1, data=umf)
fm_p5 <- occu(~tpi500 ~1, data=umf)
fm_p6 <- occu(~cc500 ~1, data=umf)
fm_p7 <- occu(~cc100 ~1, data=umf)
fm_p8 <- occu(~lure ~1, data=umf)
fm_p9 <- occu(~protected3 ~1, data=umf)
fm_p10 <- occu(~camera ~1, data=umf)
fm_p11 <- occu(~pplcat3a ~1, data=umf)
fm_p12 <- occu(~motorcat3a ~1, data=umf)
fm_p13 <- occu(~motorcat5 ~1, data=umf)
fm_p14 <- occu(~allhumancat5 ~1, data=umf)
fm_p15 <- occu(~pplcat2a ~1, data=umf)
fm_p16 <- occu(~motorcat2a ~1, data=umf)
fm_p17 <- occu(~allhumancat2a ~1, data=umf)


## -------------------------------------------------------------------------------------------------------
fms.all.p <- fitList(fits=list('psi(.)p(.)'=fm_p1, 'psi(.)p(d2road)'=fm_p2,'psi(.)p(tpi20)'=fm_p3,'psi(.)p(tpi100)'=fm_p4,'psi(.)p(tpi500)'=fm_p5,'psi(.)p(cc500.s)'=fm_p6,'psi(.)p(cc100)'=fm_p7,'psi(.)p(lure)'=fm_p8,'psi(.)p(protected3)'=fm_p9,'psi(.)p(camera)'=fm_p10,'psi(.)p(pplcat3a)'=fm_p11,'psi(.)p(motorcat3a)'=fm_p12,'psi(.)p(motorcat5)'=fm_p13,'psi(.)p(allhumancat5)'=fm_p14,'psi(.)p(pplcat2a)'=fm_p15,'psi(.)p(motorcat2a)'=fm_p16,'psi(.)p(allhumancat2a)'=fm_p17))

modSel(fms.all.p)


## -------------------------------------------------------------------------------------------------------
summary(fm_p9)


## -------------------------------------------------------------------------------------------------------
top_p=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~1, data=umf)
summary(top_p)


## -------------------------------------------------------------------------------------------------------
fm1=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~1, data=umf)
fm2=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~northness, data=umf)
fm3=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~eastness, data=umf)
fm4=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~slope, data=umf)
fm5=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~elev, data=umf)
fm6=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~elev +elev2, data=umf)
fm7=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~d2road, data=umf)
fm8=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~tpi20, data=umf)
fm9=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~tpi100, data=umf)
fm10=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~tpi500, data=umf)
fm12=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~cc500, data=umf)
fm13=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~cc100, data=umf)
fm14=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~burns20, data=umf)
fm15=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~burns100, data=umf)
fm16=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~burns500, data=umf)
fm17=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~regen20, data=umf)
fm18=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~regen100, data=umf)
fm19=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~regen500, data=umf)
fm20=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~cuts20, data=umf)
fm21=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~cuts100, data=umf)
fm22=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~cuts500, data=umf)
fm23=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~NDVIAug500, data=umf)
fm24=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~NDVIAug20, data=umf)
fm25=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~NDVIJul500, data=umf)
fm26=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~NDVIJul100, data=umf)
fm27=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~NDVIJul20, data=umf)
fm28=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~NDVIAug100, data=umf)
fm29=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhicum20, data=umf)
fm30=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhicum100, data=umf)
fm31=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhicum500, data=umf)
fm32=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhimin20, data=umf)
fm33=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhimin100, data=umf)
fm34=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhimin500, data=umf)
fm35=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhiseas20, data=umf)
fm36=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhiseas100, data=umf)
fm37=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~dhiseas500, data=umf)
fm38=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~protected2, data=umf)
fm39=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~protected3, data=umf)
fm40=occu(~protected3 + pplcat5 + tpi20 + camera +cc100~pplcat5, data=umf)


## -------------------------------------------------------------------------------------------------------
fms.all.psi <- fitList(fits=list('psi(1)p(...)'=fm1,'psi(northness)p(...)'=fm2,'psi(eastness)p(...)'=fm3,'psi(slope)p(...)'=fm4,'psi(elev)p(...)'=fm5,'psi(elev+elev2)p(...)'=fm6,'psi(d2road)p(...)'=fm7,'psi(tpi20)p(...)'=fm8,'psi(tpi100)p(...)'=fm9,'psi(tpi500)p(...)'=fm10, 'psi(cc500)p(...)'=fm12,'psi(cc100)p(...)'=fm13,'psi(burns20)p(...)'=fm14,'psi(burns100)p(...)'=fm15,'psi(burns500)p(...)'=fm16,'psi(regen20)p(...)'=fm17,'psi(regen100)p(...)'=fm18,'psi(regen500)p(...)'=fm19,'psi(cuts20)p(...)'=fm20,'psi(cuts100)p(...)'=fm21,'psi(cuts500)p(...)'=fm22,'psi(NDVIAug500)p(...)'=fm23,'psi(NDVIAug20)p(...)'=fm24,'psi(NDVIJul500)p(...)'=fm25,'psi(NDVIJul100)p(...)'=fm26,'psi(NDVIJul20)p(...)'=fm27,'psi(NDVIAug100)p(...)'=fm28,'psi(dhicum20)p(...)'=fm29,'psi(dhicum100)p(...)'=fm30,'psi(dhicum500)p(...)'=fm31,'psi(dhimin20)p(...)'=fm32,'psi(dhimin100)p(...)'=fm33,'psi(dhimin500)p(...)'=fm34,'psi(dhiseas20)p(...)'=fm35,'psi(dhiseas100)p(...)'=fm36,'psi(dhiseas500)p(...)'=fm37,'psi(protected2)p(...)'=fm38,'psi(protected3)p(...)'=fm39,'psi(pplcat5)p(...)'=fm40))
modSel(fms.all.psi)


## -------------------------------------------------------------------------------------------------------
# psi top 
occu(~protected3 + pplcat5 + tpi20 + camera +cc100 ~elev + elev2 + pplcat5 + dhiseas20 + NDVIJul20 + dhicum20 + cuts500 + burns500 + slope, data=umf) # neither burns, dhi's, nor elevation really significant

occu(~protected3 + pplcat5 + tpi20 + camera +cc100 ~elev + pplcat5 + NDVIJul20 + cuts500 + burns500 + slope, data=umf) # neither burns, dhi's, nor elevation really significant

# top model
top_fm=occu(~protected3 + pplcat5 + tpi20 + camera +cc100 ~elev + pplcat5 + NDVIJul20 + cuts500 + burns500 + slope, data=umf)
top_fm


## -------------------------------------------------------------------------------------------------------
## Detection Covariates
confint(top_fm, type="det", method = "profile") # to get CI around any beta coeff
## Psi
confint(top_fm, type="state", method = "profile") # to get CI around any beta coeff


## -------------------------------------------------------------------------------------------------------
(system.time(pb <- parboot(top_fm, statistic=chisq, nsim=100, report=10))) # if p>0.05, fail to reject null, and model fit is good
pb # p=
plot(pb)


## -------------------------------------------------------------------------------------------------------
?ranef

ELKranef = ranef(top_fm)
str(ELKranef)
hist(ELKranef@post)
head(ELKranef@post)

summary.y2=as.data.frame(t(apply(ydata7[,-1],1,na01.fnc)))
ELKoccu=cbind(data.frame(location=ydata7$location, easting=covar$easting, northing=covar$northing), summary.y2, data.frame(ELKpsi=ELKranef@post[,2,], ranefmode=bup(ELKranef , stat="mode"), ranefmean=bup(ELKranef , stat="mean")))
sum(ELKoccu$psinaive) # 250/698 = 
sum(ELKoccu$ranefmode) # 255/698 = 
sum(ELKoccu$ranefmean) # 276.6746/698 = 
Nocc(top_fm)

head(ELKoccu)


## -------------------------------------------------------------------------------------------------------
ggplot(ELKoccu, aes(easting, northing, colour =psinaive, size = psinaive)) + geom_point()


## -------------------------------------------------------------------------------------------------------
ggplot(ELKoccu, aes(easting, northing, colour =ELKpsi, size = ELKpsi)) + geom_point()

ggplot(ELKoccu, aes(easting, northing, colour =ELKpsi, size = as.factor(psinaive))) + geom_point()


## -------------------------------------------------------------------------------------------------------
estimate.of.Nocc <- Nocc(top_fm)
estimate.of.Nocc #same as sum(ranefmean)
system.time(pb.N <- parboot(top_fm, Nocc, nsim=100, report=10))   
# 100 Takes less  time (7min for 1000)
plot(pb.N)
abline(v=250, col="red")
summary(pb.N@t.star) #again same as sum(ranefmean)
quantile(pb.N@t.star, prob = c(0.025, 0.975)) #435.5439 494.5864


## -------------------------------------------------------------------------------------------------------
summary(top_fm)


## -------------------------------------------------------------------------------------------------------
fm_temp=occu(~lure + camera ~elev, data=umf)
newData=data.frame("elev"=0, "lure"=1, "camera"="flash")
predict(fm_temp, type = 'state', newdata = newData)


## -------------------------------------------------------------------------------------------------------
#rm(list=ls())
data<-read.csv("Data/blgr.csv")
head(data)
hist(data$Count1)
hist(data$Count2)
hist(data$Count3)


## -------------------------------------------------------------------------------------------------------
y<-data[,2:4]
n<-nrow(data)
#site level (individual) covariates
blgr.site<-data[,5:9]


## -------------------------------------------------------------------------------------------------------
time<-as.factor(rep(c(1,2,3),n))
blgr.obs<-data.frame(time)


## -------------------------------------------------------------------------------------------------------
blgr <- unmarkedFrameOccu(y = y, siteCovs = blgr.site,obsCovs=blgr.obs)
#summary of unmarked data frame
summary(blgr)
head(blgr)


## -------------------------------------------------------------------------------------------------------
?occuRN
rn1<-occuRN(~1 ~1,blgr, K=150)
summary(rn1)


## -------------------------------------------------------------------------------------------------------
backTransform(rn1,'det')
backTransform(rn1,"state")


## -------------------------------------------------------------------------------------------------------
re <- ranef(rn1)
plot(re)
ebup <- bup(re, stat="mean")
rn1EBUP <- sum(ebup)
rn1_CI <- confint(re,level=0.95)
rn1_CI


## -------------------------------------------------------------------------------------------------------
#time specific detection, constant occupancy 
rn2<-occuRN(~time ~1,blgr)
backTransform(rn2,"state")


## -------------------------------------------------------------------------------------------------------
rn3<-occuRN(~1 ~BQI.1.yes,blgr)
summary(rn3)
backTransform(linearComb(rn3, c(1, 0), type="state"))
backTransform(linearComb(rn3, c(0, 1), type="state"))
predict(rn3, type ="state")


## -------------------------------------------------------------------------------------------------------
rn4<-occuRN(~time ~BQI.1.yes,blgr)
summary(rn4)


## -------------------------------------------------------------------------------------------------------
rn5 <-occuRN(~time ~BQI.1.yes + Crop.history, blgr)
summary(rn5)

rn5AIC <-rn5@AIC
rn4AIC <- rn4@AIC
rn3AIC <- rn3@AIC
rn2AIC <- rn2@AIC
rn1AIC <- rn1@AIC


## -------------------------------------------------------------------------------------------------------
modelsAIC <- c(rn1AIC, rn2AIC, rn3AIC, rn4AIC, rn5AIC)
modelsAIC


## -------------------------------------------------------------------------------------------------------
predict(rn3, type ="state")

N_rn3 <- predict(rn3, type = "state")
str(N_rn3)
hist(N_rn3$Predicted)

re5 <- ranef(rn5)
plot(re5)


## -------------------------------------------------------------------------------------------------------
data2<-read.csv("Data/blgr.csv")
head(data2)


## -------------------------------------------------------------------------------------------------------
counts<-data2[,10:12]
hist(counts$Count1)
n<-nrow(data2)


## -------------------------------------------------------------------------------------------------------
blgr_count.site<-data2[,5:9]


## -------------------------------------------------------------------------------------------------------
time<-as.factor(rep(c(1,2,3),n))
blgr.obs<-data.frame(time)


## -------------------------------------------------------------------------------------------------------
blgr_count.obs<-data.frame(time)


## -------------------------------------------------------------------------------------------------------
blgr_count <- unmarkedFramePCount(y = counts, siteCovs = blgr_count.site,obsCovs=blgr_count.obs)
summary(blgr_count)


## -------------------------------------------------------------------------------------------------------
pc1<-pcount(~1 ~1,blgr_count, K = 500)
pc1

backTransform(pc1,"state")
backTransform(pc1,"det")

sum(counts$Count1)


## -------------------------------------------------------------------------------------------------------
Nmix.gof.test(pc1, nsim = 50, plot.hist = TRUE)


## -------------------------------------------------------------------------------------------------------
pc2<-pcount(~time ~1,blgr_count, K = 150)
backTransform(pc2,"state")
Nmix.gof.test(pc2, nsim = 50, plot.hist = TRUE)


## -------------------------------------------------------------------------------------------------------
pc3<-pcount(~time ~ Field.size,blgr_count, K = 150)
pc4<-pcount(~time ~BQI.1.yes, blgr_count, K = 150)
pc5<-pcount(~time ~Crop.history, blgr_count, K = 150)


## -------------------------------------------------------------------------------------------------------
pc5AIC <-pc5@AIC
pc4AIC <- pc4@AIC
pc3AIC <- pc3@AIC
pc2AIC <- pc2@AIC
pc1AIC <- pc1@AIC

modelsAIC <- c(pc1AIC, pc2AIC, pc3AIC, pc4AIC, pc5AIC)
pc5


## -------------------------------------------------------------------------------------------------------
backTransform(linearComb(pc5, coefficients = c(1, 0, 0, 0), type="state"))


## -------------------------------------------------------------------------------------------------------
backTransform(linearComb(pc5, coefficients = c(0, 0, 1, 0), type="state"))


## -------------------------------------------------------------------------------------------------------
NmixPred <- predict(pc5, type = "state")
hist(NmixPred$Predicted)


## -------------------------------------------------------------------------------------------------------
blgr_occm1 <- occu(~time ~Crop.history, data = blgr)
blgr_occm1
backTransform(linearComb(blgr_occm1, coefficients = c(1, 0, 0, 0),type="state"))
blgrPsi <- predict(blgr_occm1, type="state")

## Bind
blgr_Pred <- rbind(NmixPred$Predicted, N_rn3$Predicted)
plot(NmixPred$Predicted ~N_rn3$Predicted)


## ----eval=FALSE, include=FALSE--------------------------------------------------------------------------
## knitr::purl(input = "README.Rmd", output = "lab12.R", documentation = 1)

