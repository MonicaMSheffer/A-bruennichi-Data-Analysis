#this script should be run after the 00_climate_data_wrangling.R script, to load all of the data frames for final analysis for the manuscript. The script combines various primary data files, cleans up dates and names, and generally prepares them for formal analysis in later scripts

###load libraries
library(openxlsx)

###import data and check all column names, typos, filter out problem data, etc

##import data for analysis of adult phenotypic and life history traits across latitude
ddrad_females=read.csv("data/rad_phenotypes_completed_sortedLat.csv")
linda_data=read.xlsx("data/Linda_AllData_SpiderlingLegLength_MotherLegLength.xlsx")

#merge data into one data frame
index_mothers_ddrad=match(linda_data$Mother,ddrad_females$altID)
ddrad_phenotypesOffspring=cbind(linda_data,ddrad_females[index_mothers_ddrad,])

#make country names English
ddrad_phenotypesOffspring$Country[which(ddrad_phenotypesOffspring$Country=="Deutschland")]="Germany"
ddrad_phenotypesOffspring$Country[which(ddrad_phenotypesOffspring$Country=="Estland")]="Estonia"
ddrad_phenotypesOffspring$Country[which(ddrad_phenotypesOffspring$Country=="Frankreich")]="France"

#import processed climate data (see R script 00_climate_data_wrangling.R)
climate=read.csv("output/phenotype_environmentPCA_table.csv")
index_mother_climate=match(ddrad_phenotypesOffspring$ddradID,climate$ddradID)
ddrad_phenotypesOffspring=cbind(ddrad_phenotypesOffspring,climate[index_mother_climate,21:39])

###pigmentation data from Linda not ideal for analysis; better to do light/dark values in a binomial way
mothers_pigmentation=read.csv("data/SpiderButtResults.csv")
index_mothers_pigmentation=match(ddrad_phenotypesOffspring$Mother,mothers_pigmentation$MotherID)
ddrad_phenotypesOffspring=cbind(ddrad_phenotypesOffspring,mothers_pigmentation[index_mothers_pigmentation,3:13])

#now just get one value for mothers
ddrad_phenotypesMothers=ddrad_phenotypesOffspring[!duplicated(ddrad_phenotypesOffspring$Mother),]

##reciprocal transplant data import

#count data
count=read.csv("data/dat_post_count_wMotherLL.csv")
str(count) 
levels(as.factor(count$Opened.by)) #check that everyone wrote their name consistently - looks good

#body size data
size=read.csv("data/dat_post_size_wMotherLL.csv")
str(size)
levels(as.factor(size$Opened.by)) #looks good
levels(as.factor(size$Weighed.by)) #looks good
levels(as.factor(size$Photo.by)) #Julia Balk without space
size[size == "J.Balk"] <- "J. Balk"
levels(as.factor(size$Measured.by)) #Alex has without space
size[size == "A.Machnis"] <- "A. Machnis"

#supercooling data
scp=read.csv("data/dat_post_scp_wMotherLL.csv")
str(scp)
levels(as.factor(scp$Opened.by))
levels(as.factor(scp$Photo.by))
levels(as.factor(scp$Measured.by))
scp[scp == "L.Zander"] <- "L. Zander"
scp[scp == "A.Machnis"] <- "A. Machnis"
scp[scp == "J.Balk"] <- "J. Balk"
scp$Exclude[which(scp$Exclude!="yes")]="no" #some data points need to be excluded later based on weird peaks in the temperature profile or early differences in the protocol as the method was being developed

#chill coma data
ccr=read.csv("data/dat_post_ccr_wMotherLL.csv")
str(ccr)
levels(as.factor(ccr$Opened.by))
levels(as.factor(ccr$Photo.by))
levels(as.factor(ccr$Measured.by))
levels(as.factor(ccr$Video.watched.by))
ccr[ccr == "A.Machnis"] <- "A. Machnis"
ccr[ccr == "J.Balk"] <- "J. Balk"
ccr[ccr == "a.Machnis"] <- "A. Machnis"

#lower lethal temperature data
llt=read.csv("data/dat_post_llt_wMotherLL.csv")
str(llt)
levels(as.factor(llt$Opened.by)) #all look good

#Fix all to have dates as dates, not characters
count$Oviposition_date=as.Date(count$Oviposition_date,format="%d-%m-%y")
count$Planned.opening.Date=as.Date(count$Planned.opening.Date,format="%d-%m-%y")
count$Actual.opening.date=as.Date(count$Actual.opening.date,format="%d-%m-%y")
count$CollectionDate=as.Date(count$CollectionDate,format="%d-%m-%y")

size$Oviposition_date=as.Date(size$Oviposition_date,format="%d-%m-%y")
size$Planned.opening.Date=as.Date(size$Planned.opening.Date,format="%d-%m-%y")
size$Actual.opening.date=as.Date(size$Actual.opening.date,format="%d-%m-%y")

scp$Oviposition_date=as.Date(scp$Oviposition_date,format="%d-%m-%y")
scp$Planned.opening.Date=as.Date(scp$Planned.opening.Date,format="%d-%m-%y")
scp$Actual.opening.date=as.Date(scp$Actual.opening.date,format="%d-%m-%y")

ccr$Oviposition_date=as.Date(ccr$Oviposition_date,format="%d-%m-%y")
ccr$Planned.opening.Date=as.Date(ccr$Planned.opening.Date,format="%d-%m-%y")
ccr$Actual.opening.date=as.Date(ccr$Actual.opening.date,format="%d-%m-%y")

llt$Oviposition_date=as.Date(llt$Oviposition_date,format="%d-%m-%y")
llt$Planned.opening.Date=as.Date(llt$Planned.opening.Date,format="%d-%m-%y")
llt$Actual.opening.date=as.Date(llt$Actual.opening.date,format="%d-%m-%y")

#create variable for all sheets that is the number of days spent within winter treatment (started winter on October 1st 2018)

count$days_inTreatment=as.numeric(count$Actual.opening.date-as.Date("01/10/2018",format="%d-%m-%y"),units="days")
size$days_inTreatment=as.numeric(size$Actual.opening.date-as.Date("01/10/2018",format="%d-%m-%y"),units="days")
ccr$days_inTreatment=as.numeric(ccr$Actual.opening.date-as.Date("01/10/2018",format="%d-%m-%y"),units="days")
scp$days_inTreatment=as.numeric(scp$Actual.opening.date-as.Date("01/10/2018",format="%d-%m-%y"),units="days")

#calculate precise age upon opening for each egg sac
count$openingAge=as.numeric(count$Actual.opening.date-count$Oviposition_date,units="days")
mean(count$openingAge[which(count$round=="post")])
sd(count$openingAge[which(count$round=="post")])

#change winter treatment to "warm" (France) or "cold" (Estonia) for clarity
count$Winter.treatment[which(count$Winter.treatment=="France")]="warm"
count$Winter.treatment[which(count$Winter.treatment=="Estonia")]="cold"

scp$Winter.treatment[which(scp$Winter.treatment=="France")]="warm"
scp$Winter.treatment[which(scp$Winter.treatment=="Estonia")]="cold"

size$Winter.treatment[which(size$Winter.treatment=="France")]="warm"
size$Winter.treatment[which(size$Winter.treatment=="Estonia")]="cold"

ccr$Winter.treatment[which(ccr$Winter.treatment=="France")]="warm"
ccr$Winter.treatment[which(ccr$Winter.treatment=="Estonia")]="cold"

llt$Winter.treatment[which(llt$Winter.treatment=="France")]="warm"
llt$Winter.treatment[which(llt$Winter.treatment=="Estonia")]="cold"

#change origin to "core" (France) or "edge" (Latvia/Estonia) for clarity

count$Origin[which(count$Origin=="France")]="core"
count$Origin[which(count$Origin=="Estonia")]="edge"

scp$Origin[which(scp$Origin=="France")]="core"
scp$Origin[which(scp$Origin=="Estonia")]="edge"

size$Origin[which(size$Origin=="France")]="core"
size$Origin[which(size$Origin=="Estonia")]="edge"

ccr$Origin[which(ccr$Origin=="France")]="core"
ccr$Origin[which(ccr$Origin=="Estonia")]="edge"

llt$Origin[which(llt$Origin=="France")]="core"
llt$Origin[which(llt$Origin=="Estonia")]="edge"

