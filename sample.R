##########
#
#   Analysis of summary trip data
#     March 26, 2016
#     Updated in Aug/Sept for Barcelona meeting 2017
#     Makes use of data prepped by other code that turned raw into trips
#   July 10, 2018 - this is what I used for Barcelona
#     Looks as if the significant effect was on HR, not Area, per se
#     Also looks as if there is too much Kurtosis - so can't really use it
#     If you get rid of Forage_Area outliers, you don't have the interaction effect
#
##########

#### remember to run the Prep Script to setup some other things you will need later

### first, setup all the libraries you need. Make sure you ahve the packages loaded
library(adehabitatLT)
library(adehabitatHR)
library(spatstat)

### other libraries needed
library(raster) 
library(rgdal)
library(gpclib) ### must be before maptools to give permission to use something in maptools
library(maptools)
library(ggplot2)
library(tkrplot)  ### having problems loading this

#### for distance, etc. from calibration code
library(geosphere)   ### this allows us to make SpatialPoints
library(circular)   ## for doing the circular stats

#### for basic statistical 
library(lme4)
library(psych)

### for plotting effects of a random effects model
library(effects)




### load data NOTE: 9/28/17 these are processed trips, not the raw files
All_Birds <- read.csv("~/Dropbox/Research/RStudio Files/Spatial_Geos/Data/TripFormat_Mar26.csv")

### updated with "U" for unknown sex birds
All_Birds <- read.csv("~/Dropbox/Research/RStudio Files/Spatial_Geos/Data/TripFormat_Mar28.csv")

#### make a density variable D = N/A - something wrong ith this - Sat night
All_Birds$FDensity = (All_Birds$Num_Points-2) / All_Birds$Forage_Area
All_Birds$Density = All_Birds$Num_Points / All_Birds$Forage_HR


#### prep the data - make sure they are factors, if not make them so
is.factor(All_Birds$Sex)

is.factor(All_Birds$Year)
All_Birds$Year = as.factor(All_Birds$Year)

is.factor(All_Birds$BBN)

is.factor(All_Birds$Burrow)
All_Birds$Burrow = as.factor(All_Birds$Burrow)

is.factor(All_Birds$Days_Post_Lay)
All_Birds$Days_Post_Lay = as.numeric(All_Birds$Days_Post_Lay)
is.numeric(All_Birds$Days_Post_Lay)
All_Birds$Speed_Return = as.numeric(All_Birds$Speed_Return)

### make sure the days post lay are correct
All_Birds$Days_Post_Lay = All_Birds$J_Trip_Start_Date - All_Birds$J_Lay_Date

### assign a short name to the dataset so we can leave the original unchanged
df <- All_Birds

#### Examine data
mNormCheck(df$Ptilo)
mNormCheck(df$J_Lay_Date)
mNormCheck(df$Max_Dist)
mNormCheck(df$Forage_HR)
mNormCheck(df$Speed_Forage)
mNormCheck(df$Speed_Leave)
mNormCheck(df$Foraging_Dist)
mNormCheck(df$Forage_Int) #### some have 0 because not enough points for the function
mNormCheck(df$fDensity)
##### simple tests



## get rid of outliers
df1 <- subset(All_Birds, Foraging_Dist < 3500)
df2 <- subset(df1, Sex != "U")
df2 <- subset(df, Sex != "U")
### get rid of the NA for Lay date
df3 <- subset(df2, BBN != "7901-11024") 
### call it df2 again for running models
df2 <- df3


#### mixed-effects models - use lme4 package to do this

#############
# Foraging Dist - summed all legs
#############
fdNULL <- lmer(Foraging_Dist ~ 1 + (1|Year) +  (1|BBN), data = df2)
fd01 <- lmer(Foraging_Dist ~ 1  + Sex + (1|Year) + (1|BBN), data = df2)
fd02 <- lmer(Foraging_Dist ~ 1  + Sex + J_Trip_Start_Date  + (1|Year) + (1|BBN), data = df2)   #### Best model here - both are highly significant
fd03 <- lmer(Foraging_Dist ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|Year) + (1|BBN), data = df2)
fd04 <- lmer(Foraging_Dist ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fd05 <- lmer(Foraging_Dist ~ 1  + Sex + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fd06 <- lmer(Foraging_Dist ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fd07 <- lmer(Foraging_Dist ~ 1  + J_Trip_Start_Date + (1|Year) + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(fdNULL, fd01, fd02, fd03, fd04, fd05, fd06, fd07)
summary(fd05) ## best one
mcheck(fd05)
plot(effect(term="Sex", mod=fd05, default.levels=2), multiline=TRUE,
            main = "", xlab = "", ylab = "Foraging Distance (km)"
            , cex = 2)  ## colors = c("blue", "red"))


plot(effect(term="Days_Post_Lay", mod=fd05, default.levels=2), multiline=TRUE,
     main = "", xlab = "Days Post Lay Date", ylab = "Foraging Distance (km)"
     , cex = 2)  ## colors = c("blue", "red"))

mcheck(fd06)
summary(fd06) ## best one almost
plot(effect(term="Sex:J_Trip_Start_Date", mod=fd03, default.levels=2), multiline=TRUE)
### plot them individually?
plot(allEffects(fd03), main = "", xlab = "Start Date (Day of Year)", ylab = "Foraging Distance (km)")


### just plot 
plot(effect(term="J_Trip_Start_Date", mod=fd02, default.levels=2), multiline=TRUE)
plot(effect(term="Days_Post_Lay", mod=fd06, default.levels=2), multiline=TRUE)

d#############
# Core Dist
#############

fcdNULL <- lmer(Forage_Core_Dist ~ 1 + (1|Year) + (1|BBN), data = df2)
fcd01 <- lmer(Forage_Core_Dist ~ 1  + Sex + (1|BBN), data = df2)
fcd02 <- lmer(Forage_Core_Dist ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fcd03 <- lmer(Forage_Core_Dist ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fcd04 <- lmer(Forage_Core_Dist ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fcd05 <- lmer(Forage_Core_Dist ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fcd06 <- lmer(Forage_Core_Dist ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(fcdNULL, fcd02, fcd03, fcd04, fcd05, fcd06)
summary(fcd05) ## best one - not much going on
summary(fcd02) ## second best

plot(effect(term="Days_Post_Lay", mod=fcd05, default.levels=2), multiline=TRUE,
     main = "", xlab = "Days Post Lay Date", ylab = "Mean Distance from KI (km)"
     , cex = 2)  ## colors = c("blue", "red"))

mcheck(fcd02)
summary(fcd03) ## best one almost
plot(effect(term="Sex:Days_Post_Lay", mod=fcd06, default.levels=2), multiline=TRUE)
plot(allEffects(fcd06), main = "", xlab = "Days Post Lay Date", ylab = "Core Distance from KI (km)")

#### do just core size

d#############
# Core Size
#############

fcdNULL <- lmer(Forage_Core ~ 1 + (1|Year) + (1|BBN), data = df2)
fcd01 <- lmer(Forage_Core ~ 1  + Sex + (1|BBN), data = df2)
fcd02 <- lmer(Forage_Core ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fcd03 <- lmer(Forage_Core ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fcd04 <- lmer(Forage_Core ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fcd05 <- lmer(Forage_Core ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fcd06 <- lmer(Forage_Core ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(fcdNULL, fcd02, fcd03, fcd04, fcd05, fcd06)
summary(fcd04) ## best one - not much going on
summary(fcd05) ## second best
plot(effect(term="Sex", mod=fcd02, default.levels=2), multiline=TRUE)
mcheck(fcd02)
summary(fcd03) ## best one almost

plot(allEffects(fcd03), main = "", xlab = "Start Date (Day of Year)", ylab = "Core Distance from KI (km)")

d#############
# Number of Pofnps on trip - lenght of trip
#############

fnpNULL <- lmer(Num_Points ~ 1 + (1|Year) + (1|BBN), data = df2)
fnp01 <- lmer(Num_Points ~ 1  + Sex + (1|BBN), data = df2)
fnp02 <- lmer(Num_Points ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fnp03 <- lmer(Num_Points ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fnp04 <- lmer(Num_Points ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fnp05 <- lmer(Num_Points ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fnp06 <- lmer(Num_Points ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(fnpNULL, fnp02, fnp03, fnp04, fnp05, fnp06)
summary(fnp05) ## best one - not much going on
summary(fnp04) ## second best
summary(fnp06) ## third and all 3 good
plot(effect(term="Sex", mod=fnp02, default.levels=2), multiline=TRUE)
mcheck(fnp02)
summary(fnp03) ## best one almost

plot(allEffects(fnp05), main = "", xlab = "Days Post Lay Date", ylab = "Time Away from KI (km)")


### where they go doesn't change between them, but how long they stay changes by sex with time

#############
# Foraging Speed
#############
fspNULL <- lmer(Speed_Forage ~ 1 + (1|Year) + (1|BBN), data = df2)
fsp01 <- lmer(Speed_Forage ~ 1  + Sex + (1|BBN), data = df2)
fsp02 <- lmer(Speed_Forage ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fsp03 <- lmer(Speed_Forage ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fsp04 <- lmer(Speed_Forage ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fsp05 <- lmer(Speed_Forage ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fsp06 <- lmer(Speed_Forage ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

anova(fspNULL, fsp02, fsp03, fsp04, fsp05, fsp06)
summary(fsp06) ## best one
mcheck(fsp06)  ### looks good
plot(effect(term="Sex:Days_Post_Lay", mod=fsp06, default.levels=2), multiline=TRUE)  ### interaction plot is amazing!
plot(allEffects(fsp06), main = "", xlab = "Days Post Lay Date", ylab = "Foraging Speed (km/hr)")

plot(allEffects(fsp03), main = "", xlab = "Start Date (Day of Year)", ylab = "Foraging Speed (km/hr)")

#############
# Home Range (95% mcp)
#############

fhrNULL <- lmer(Forage_HR ~ 1 + (1|Year) + (1|BBN), data = df2)
fhr01 <- lmer(Forage_HR ~ 1  + Sex + (1|BBN), data = df2)
fhr02 <- lmer(Forage_HR ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fhr03 <- lmer(Forage_HR ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fhr04 <- lmer(Forage_HR ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fhr05 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fhr06 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

anova(fhrNULL, fhr02, fhr03, fhr04, fhr05, fhr06)
summary(fhr06) ## best one
mcheck(fhr06)  ## bad residuals!
summary(fhr02) ## next best one
mcheck(fhr02)  ## bad residuals
plot(effect(term="Sex:Days_Post_Lay", mod=fhr06, default.levels=2), multiline=TRUE)  ### interaction plot - very cool result
plot(allEffects(fhr06), main = "", xlab = "Days Post Lay Date", ylab = "Home Range (km2)")


#############
# Forage HR
#############

faNULL <- lmer(Forage_HR ~ 1 + (1|Year) + (1|BBN), data = df2)
fa01 <- lmer(Forage_HR ~ 1  + Sex + (1|BBN), data = df2)
fa02 <- lmer(Forage_HR ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fa03 <- lmer(Forage_HR ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fa04 <- lmer(Forage_HR ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fa05 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fa06 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)
fa07 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN) + (1|Year), data = df2)


anova(faNULL, fa02, fa03, fa04, fa05, fa06)
summary(fa06) ## best one
mcheck(fa06)  ## bad residuals!
summary(fa02) ## next best one
summary(fa06) ## look at post lay date
mcheck(fa02)  ## bad residuals
plot(effect(term="Sex", mod=fa06, default.levels=2), multiline=TRUE, main = "", ylab = "Foraging Area (km2)")  ### interaction plot - very cool result
plot(allEffects(fa06), main = "", xlab = "Days Post Lay Date", ylab = "Foraging Area (km2)")



##
#############
# Forage Area - JULY 2018 - Noticed it is HOME RANGE, not AREA! duplicate in next section
#############

fareaNULL <- lmer(Forage_HR ~ 1 + (1|Year) + (1|BBN), data = df2)
farea01 <- lmer(Forage_HR ~ 1  + Sex + (1|BBN), data = df2)
farea02 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
farea03 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)
farea04 <- lmer(Forage_HR ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN) + (1|Year), data = df2)


anova(fareaNULL, farea02, farea03, farea04, farea01)
summary(farea03) ## best one
mcheck(farea03)  ## bad residuals!
summary(farea02) ## next best one
summary(farea06) ## look at post lay date
mcheck(farea02)  ## bad residuals
plot(effect(term="Sex", mod=farea03, default.levels=2), multiline=TRUE, main = "", ylab = "Foraging Area (km2)")  ### interaction plot - very cool result
plot(allEffects(farea06), main = "", xlab = "Days Post Lay Date", ylab = "Foraging Area (km2)")

### plot just the days effect
plot(effect(term="Days_Post_Lay", mod=farea03, default.levels=2), multiline=TRUE,
     main = "", xlab = "Days Post Lay Date", ylab = "Foraging Area (km2)"
     , cex = 2)  ## colors = c("blue", "red"))


##
#############
# Forage Area - JULY 2018 - redo above with Forage_Area, not HR
#             - added b/c revisiting the Kent Island only data
#############

#### 
# 2SD above mean = 136000, get rid of those for this analysis
##
df2a <- subset(df2, Forage_Area <= 136000)


FA18NULL <- lmer(Forage_Area ~ 1 + (1|Year) + (1|BBN), data = df2a)
FA1801 <- lmer(Forage_Area ~ 1  + Sex + (1|BBN), data = df2a)
FA1802 <- lmer(Forage_Area ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2a)
FA1803 <- lmer(Forage_Area ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2a)
FA1804 <- lmer(Forage_Area ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN) + (1|Year), data = df2a)


anova(FA18NULL, FA1802, FA1803, FA1804, FA1801)
summary(FA1803) ## best one
mcheck(FA1803)  ## bad residuals!
summary(FA1802) ## next best one
mcheck(FA1802)  ## bad residuals
plot(effect(term="Sex", mod=FA1803, default.levels=2), multiline=TRUE, main = "", ylab = "Foraging Area (km2)")  ### interaction plot - very cool result
plot(allEffects(FA1803), main = "", xlab = "Days Post Lay Date", ylab = "Foraging Area (km2)")

### plot just the days effect
plot(effect(term="Days_Post_Lay", mod=FA1803, default.levels=2), multiline=TRUE,
     main = "", xlab = "Days Post Lay Date", ylab = "Foraging Area (km2)"
     , cex = 2)  ## colors = c("blue", "red"))

########
#  END of July 28, revisit experience
#######################




#############
# Forage Intensity - Beautiful!
#############

### get rid of 0s
dfInt <- subset(df2, Forage_Int > 0)

firNULL <- lmer(Forage_Int ~ 1 + (1|Year) + (1|BBN), data = df2)
fir01 <- lmer(Forage_Int ~ 1  + Sex + (1|BBN), data = df2)
fir02 <- lmer(Forage_Int ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fir03 <- lmer(Forage_Int ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fir04 <- lmer(Forage_Int ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fir05 <- lmer(Forage_Int ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fir06 <- lmer(Forage_Int ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

anova(firNULL, fir02, fir03, fir04, fir05, fir06)
summary(fir06) ## best one
mcheck(fir06)  ## bad residuals!
summary(fir05) ## next best one
mcheck(fir02)  ## bad residuals
plot(effect(term="Sex:Days_Post_Lay", mod=fir06, default.levels=2), multiline=TRUE, ylab = "Foraging Intensity (km2)", main = "")

plot(allEffects(fir06), main = "", xlab = "Trip Start Date (Day of Year)", ylab = "Foraging Intensity (km2)")

#############
# Simplified Density (Area/N)
#############

### get rid of 0s

firNULL <- lmer(FDensity ~ 1 + (1|Year) + (1|BBN), data = df2)
fir01 <- lmer(FDensity ~ 1  + Sex + (1|BBN), data = df2)
fir02 <- lmer(FDensity ~ 1  + Sex + J_Trip_Start_Date  + (1|BBN), data = df2)   #### Best model here - both are highly significant
fir03 <- lmer(FDensity ~ 1  + Sex + J_Trip_Start_Date + J_Trip_Start_Date*Sex + (1|BBN), data = df2)
fir04 <- lmer(FDensity ~ 1  + Sex + J_Trip_Start_Date + Days_Post_Lay + (1|BBN), data = df2)
fir05 <- lmer(FDensity ~ 1  + Sex + Days_Post_Lay + (1|BBN), data = df2)
fir06 <- lmer(FDensity ~ 1  + Sex + Days_Post_Lay + Sex*Days_Post_Lay + (1|BBN), data = df2)

anova(firNULL, fir02, fir03, fir04, fir05, fir06)
summary(fir03) ## best one
mcheck(fir03)  ## bad residuals!
summary(fhr02) ## next best one
mcheck(fhr02)  ## bad residuals
plot(effect(term="Sex:J_Trip_Start_Date", mod=fir03, default.levels=2), multiline=TRUE)  ### interaction plot - very cool result
plot(allEffects(fir03), main = "", xlab = "J Date", ylab = "Forage Int")



#############
# Ptlo models
##############

#### don't work as is, because Ptilo measure is same for all trips for all birds
#### need to summarize data for each bird, then do ptilo - R do it?

### this gives mean values for columns 9-14
meanPtilo <- aggregate(df2[,2:24], list(df2$BBN), mean)

### however, I took those to Excel, filled in the non-mean values and go this
meanPtilo <- read.csv("~/Dropbox/Research/RStudio Files/Spatial_Geos/Data/ptilo01.csv")
meanPtilo$Year = as.factor(meanPtilo$Year)

ptNULL <- lmer(Ptilo ~ 1 + (1|Year), data = meanPtilo)
pt01 <- lmer(Ptilo ~ 1  + Sex + (1|Year), data = meanPtilo)
pt02 <- lmer(Ptilo ~ 1  + Sex + J_Lay_Date  + (1|Year), data = meanPtilo) 
pt03 <- lmer(Ptilo ~ 1  + Sex + J_Lay_Date + J_Lay_Date*Sex + (1|Year), data = meanPtilo)
pt04 <- lmer(Ptilo ~ 1  + J_Lay_Date  + (1|Year), data = meanPtilo) 

pt05 <- lmer(Ptilo ~ 1  + Foraging_Dist + (1|Year), data = meanPtilo)
pt06 <- lmer(Ptilo ~ 1  + Forage_Area + (1|Year), data = meanPtilo)
pt07 <- lmer(Ptilo ~ 1  + Forage_HR + (1|Year), data = meanPtilo)
pt08 <- lmer(Ptilo ~ 1  + Forage_Int + (1|Year), data = meanPtilo)
pt09 <- lmer(Ptilo ~ 1  + Sex + Foraging_Dist + (1|Year), data = meanPtilo)
pt10 <- lmer(Ptilo ~ 1  + Sex + Foraging_Dist + Sex*Foraging_Dist + (1|Year), data = meanPtilo)

anova(ptNULL, pt01, pt02, pt03, pt04, pt05, pt06, pt07, pt08, pt09, pt10)
summary(pt07) ## best one


plot(effect(term="Sex*Foraging_Dist", mod=pt10, default.levels=2), multiline=TRUE)



### simple correlations?



## however, speed and dist traveled are correlated, of course
plot(df2$Foraging_Dist ~ df2$Speed_Forage)

### jusut look at means
plot(df1$Foraging_Dist ~ df1$Sex )

### does leaving speed have anything to do with foraging speed? NO relationship
plot(df2$Speed_Forage ~ df2$Speed_Leave)

# do they leave fast, then slow down to forage? NULL is equal
t.test(df1$Speed_Leave, df1$Speed_Forage, paired = TRUE)   ### WAY NOT, they fly out fast, then slow down
t.test(df$Speed_Return, df$Speed_Forage, paired = TRUE)   ### WAY NOT, they fly out fast, then slow down
t.test(df$Speed_Return, df$Speed_Leave, paired = TRUE)   ### WAY NOT, they fly out fast, then slow down

### plot leaving speed by Sex
boxplot(Speed_Leave ~ Sex, data=df2, notch=TRUE, 
        col=(c("gold","darkgreen")),
        main="Embark Speed", xlab="X label")

boxplot(Foraging_Dist ~ Sex, data=df2, notch=TRUE, 
        col=(c("gold","darkgreen")),
        main="Foraging Dist", xlab="X label")

boxplot(Forage_HR ~ Sex, data=df2, notch=TRUE, 
        col=(c("gold","darkgreen")),
        main="Foraging Dist", xlab="X label")

boxplot(Forage_Area ~ Sex, data=df2, notch=TRUE, 
        col=(c("blue","red")),
        main="", xlab="", ylab = "Foraging Area (km2)")

###### Simple Tests
wilcox.test(df2$Speed_Forage ~ df2$Sex )
wilcox.test(df2$Forage_Area ~ df2$Sex )
wilcox.test(df2$Forage_Core ~ df2$Sex )
wilcox.test(df2$Forage_HR ~ df2$Sex )
wilcox.test(df2$Forage_Core_Dist ~ df2$Sex )
wilcox.test(df2$Foraging_Dist ~ df2$Sex )


###### models by individual sexes

#### just sub dff or dfm for df
### get only males
dfm <- subset(df2, Sex == "Male")

### get only females
dff <- subset(df2, Sex == "Female")

dfSex <- dfm
pTitle <- "Male"

#############
# Foraging Dist - summed all legs
#############
fdNULLSex <- lmer(Foraging_Dist ~ 1 + (1|Year) +  (1|BBN), data = dfSex)
fd01Sex <- lmer(Foraging_Dist ~ 1  + J_Trip_Start_Date  + (1|Year) + (1|BBN), data = dfSex)   #### Best model here - both are highly significant
fd02Sex <- lmer(Foraging_Dist ~ 1  + J_Trip_Start_Date + Days_Post_Lay + (1|Year) + (1|BBN), data = dfSex)
fd03Sex <- lmer(Foraging_Dist ~ 1  + Days_Post_Lay + (1|Year) + (1|BBN), data = dfSex)
fd04Sex <- lmer(Foraging_Dist ~ 1  + J_Trip_Start_Date + Days_Post_Lay + J_Trip_Start_Date*Days_Post_Lay + (1|Year) + (1|BBN), data = dfSex)

### compare them - can't use 01 /c not same sample size
anova(fdNULLSex, fd01Sex, fd02Sex, fd03Sex, fd04Sex)
summary(fd01Sex) ## best one
mcheck(fd01Sex)

### just plot 
plot(effect(term="J_Trip_Start_Date", mod=fd01Sex, default.levels=2), multiline=FALSE, main = pTitle)
plot(effect(term="Days_Post_Lay", mod=fd06, default.levels=2), multiline=FALSE)

#############
# Speed Leave - not correlated with Area
#############
fsplNULL <- lmer(Speed_Leave ~ 1 + (1|Year) +  (1|BBN), data = df2)
fspl01 <- lmer(Speed_Leave ~ 1 + Sex + (1|Year) +  (1|BBN), data = df2)
fspl02 <- lmer(Speed_Leave ~ 1 + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fspl03 <- lmer(Speed_Leave ~ 1 + Sex + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fspl04 <- lmer(Speed_Leave ~ 1 + Sex + + Days_Post_Lay + Sex*Days_Post_Lay + (1|Year) + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(fsplNULL, fspl01, fspl02, fspl03, fspl04)
summary(fspl01) ## best one
mcheck(fspl01)

### just plot 
pTitle = ""
yLabel = "Embark Speed (km/h)"
plot(effect(term="Sex", mod=fspl01, default.levels=2), multiline=FALSE, main = pTitle, ylab = yLabel)
plot(effect(term="Days_Post_Lay", mod=fspl01, default.levels=2), multiline=FALSE)

### show that they aren't correlated
summary(lm(FirstLeg ~ Forage_Area, data = df1))
spearmanTest(df2$FirstLeg, df2$Forage_Area) # p = 0.93

fsplCheck <- lmer(Speed_Leave ~ 1 + Forage_Area + (1|Year) +  (1|BBN), data = df2)
summary(fsplCheck)
plot(FirstLeg ~ Forage_Area, data = df2, xlab = "Foraging Area (km2)", ylab = "First Leg (km)")


#### distance first leg
df2$FirstLeg <- df2$Speed_Leave * 12
mNormCheck(df2$FirstLeg)


df2$Days_at_Sea <- df2$Num_Points - 2
### do analyses
#############
# Num_Points - not correlated with Area
#############
PtsNULL <- lmer(Days_at_Sea ~ 1 + (1|Year) +  (1|BBN), data = df2)
Pts01 <- lmer(Days_at_Sea ~ 1 + Sex + (1|Year) +  (1|BBN), data = df2)
Pts02 <- lmer(Days_at_Sea ~ 1 + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
Pts03 <- lmer(Days_at_Sea ~ 1 + Sex + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
Pts04 <- lmer(Days_at_Sea ~ 1 + Sex + + Days_Post_Lay + Sex*Days_Post_Lay + (1|Year) + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(PtsNULL, Pts01, Pts02, Pts03, Pts04)
summary(Pts02) ## best one
mcheck(Pts02)
### next best
summary(Pts03) ## best one
mcheck(Pts03)
### Plot it
yLabel = "Days Foraging"
plot(effect(term="Sex", mod=Pts03, default.levels=2), multiline=FALSE, main = pTitle, ylab = yLabel, xlab = "")

### not correlated wth area
plot(Days_at_Sea ~ Forage_Area, data = df2, xlab = "Foraging Area (km2)", ylab = "Days at Sea (x2)")


### do analyses
#############
# First Leg - not correlated with Area
#############
flegNULL <- lmer(FirstLeg ~ 1 + (1|Year) +  (1|BBN), data = df2)
fleg01 <- lmer(FirstLeg ~ 1 + Sex + (1|Year) +  (1|BBN), data = df2)
fleg02 <- lmer(FirstLeg ~ 1 + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fleg03 <- lmer(FirstLeg ~ 1 + Sex + Days_Post_Lay + (1|Year) + (1|BBN), data = df2)
fleg04 <- lmer(FirstLeg ~ 1 + Sex + + Days_Post_Lay + Sex*Days_Post_Lay + (1|Year) + (1|BBN), data = df2)

### compare them - can't use 01 /c not same sample size
anova(flegNULL, fleg01, fleg02, fleg03, fleg04)
summary(fleg01) ## best one
mcheck(fleg01)
### Plot it
yLabel = "First Leg (km)"
plot(effect(term="Sex", mod=fleg01, default.levels=2), multiline=FALSE, main = pTitle, ylab = yLabel, xlab = "")




### how are they all related
pairs(~ Forage_Area + Foraging_Dist +  Forage_Core_Dist + Max_Dist + Speed_Forage + Forage_Core,data=df2, 
      main="")

### how are they all related
pairs(~ Forage_Area + Foraging_Dist +  Forage_Core_Dist + Num_Points + Speed_Forage + Forage_Core,data=df2, 
      main="Simple Scatterplot Matrix")

### how are they all related
pairs(~ Forage_Area + Foraging_Dist +  Speed_Leave + Speed_Forage + Speed_Return, data=df2, 
      main="Simple Scatterplot Matrix")

##### Oct 13, 2016 in Moögingen at MPIO for talk

## need histogram of max distances to show effect of error
hist(All_Birds$Max_Dist, col = "grey", nclass = 20, xlab = "Distance from KI (km)", yla = "Foraging Trips", main = "")

### need correlations with main stats - get rid of outliers
plotData <- subset(All_Birds, Max_Dist < 1200)
plot(~ Forage_Area + Max_Dist + Dist_Daily + Forage_HR, data = plotData)

### now look at distribution of areas - maybe not use it
hist(All_Birds$Forage_Area, col = "grey", nclass = 20, xlab = "Forage Area (km2)", yla = "Foraging Trips", main = "")
hist(plotData$Forage_Area, col = "grey", nclass = 20, xlab = "Forage Area (km2)", yla = "Foraging Trips", main = "")


