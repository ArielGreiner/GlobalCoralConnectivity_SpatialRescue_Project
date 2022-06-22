###code below written with help from Dr. Marco Andrello.
#setwd("~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco")

#read the shape file into R, append the tables to the shape file, remove useless columns
library(rgdal)
library(shapefiles)

#read the shape file into an object (which includes the dataframe that we want, and the polygons for plotting)
a <- readOGR(dsn = "~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/016polys_score_Mollweide_with50reefs.shp")
#access the dataframe
head(a@data)
#note: a@data$reefpct is the avg of the reef percentage scores of each of the 50 reef polygons, might be worth using that to get a more accurate
#calculation of %reef SA in Sally Wood's polygons...but there are issues with 50reefpolygons being in between multiple SallyWoodgridcells and we're 
#not sure if reefpct is actually reef percent area 
range(a@data$reefpct)

maindata <- a@data[,-c(1,3,4)]
head(maindata)

#add in %reef data
percentreef <- read.dbf(dbf.name = "~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/Reef_Cover.dbf")
head(percentreef$dbf) #dataframe, FID_ = TARGET_FID in maindata, missing the ones with no reef cover
percentreeftable <- percentreef$dbf
id <- match(x = maindata$TARGET_FID, table = percentreeftable$FID_) #assigns NAs to the values in TARGET_FID that don't show up in FID_
maindata2 <- cbind(maindata,percentreef$dbf[id,]) #attach the rows with percent cover data in the right spots, NAs elsewhere
#rows with NA: which(is.na(maindata2$FID_))
#testing: maindata2[307,]
names(maindata2)[15] <- "ReefArea" #rename AREA
names(maindata2)[16] <- "ReefPctArea" #rename PERCENTAGE

maindata2$FID_ <- NULL
names(maindata2)[13] <- "BCU_Name" #rename ReefName
names(maindata2)[12] <- "BCU_ID" #rename puid

#add in %50reef data (note: BCU refers to the 50 reefs)
percent50reef <- read.dbf(dbf.name = "~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/50Reef_Cover.dbf")
head(percent50reef$dbf) #dataframe, FID_ = TARGET_FID in maindata, missing the ones with no reef cover
percent50reeftable <- percent50reef$dbf
idfifty <- match(x = maindata$TARGET_FID, table = percent50reeftable$FID_) #assigns NAs to the values in TARGET_FID that don't show up in FID_
maindata3 <- cbind(maindata2,percent50reef$dbf[idfifty,])
head(maindata3)
maindata3$FID_ <- NULL
names(maindata3)[16] <- "BCUarea" #replace AREA
names(maindata3)[17] <- "BCUpctarea" #replace PERCENTAGE
#head(maindata3)

#add in %protected data
percentprotected <- read.dbf(dbf.name = "~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/MPA_Cover.dbf")
head(percentprotected$dbf) #dataframe, FID_ = TARGET_FID in maindata, missing the ones with no reef cover
idprt <- match(x = maindata$TARGET_FID, table = percentprotected$dbf$FID_) #assigns NAs to the values in TARGET_FID that don't show up in FID_
maindata4 <- cbind(maindata3,percentprotected$dbf[idprt,])
head(maindata4)
maindata4$FID_ <- NULL
names(maindata4)[18] <- "Protectedarea" #replace AREA
names(maindata4)[19] <- "pct_Protectedarea" #replace PERCENTAGE
head(maindata4)

#add in marine realms
marinerealm <- read.dbf(dbf.name = "~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/Marine Realms Costello/016polys_score_Realms_join_1to1_bis.dbf")
head(marinerealm$dbf) #dataframe, FID_ = TARGET_FID in maindata
#12397 rows <- so all of them were assigned realms bc join count = 1 for all also, etc
#which((maindata4$TARGET_FID - marinerealm$dbf$TARGET_FID) != 0) #integer(0) so same order 
maindata5 <- maindata4
maindata5$MarineRealm <- marinerealm$dbf$Realm
maindata5$MRealm_Name <- marinerealm$dbf$Name
head(maindata5)

#add in EEZ data
EEZ <- read.dbf(dbf.name = "~/Dropbox/University of Toronto/Research Related/Sally Wood Connectivity Matrix/Ariel_connectivity_SallyWood50reefEEZWDPA_fromMarco/EEZs/016polys_score_EEZ_join_1to1.dbf") #this one concatenated the names if there were multiple
head(EEZ$dbf) #dataframe, FID_ = TARGET_FID in maindata
#12397 rows <- so all of them were assigned EEZs, sometimes join count > 1
#which((maindata5$TARGET_FID - EEZ$dbf$TARGET_FID) != 0) #integer(0) so same order 
maindata6 <- maindata5
maindata6$NumEEZ <- EEZ$dbf$Join_Count
maindata6$EEZName <- EEZ$dbf$GeoName
maindata6$EEZ_Sovereign1 <- EEZ$dbf$Sovereign1
maindata6$EEZ_Sovereign2 <- EEZ$dbf$Sovereign2
maindata6$EEZ_Sovereign3 <- EEZ$dbf$Sovereign3
head(maindata6)

#remove rows with no reef cover?
id.noReef <- which(is.na(maindata6$ReefArea))
length(id.noReef)
# There are 105 cells without reefs.
# We have to remove these cells from the dataframe but also from the shapefile (if we want to append the dataframe to the shapefile)
# So we are going to first append the dataframe, then remove the empty cells

# Because I am experiencing issues in writing values of ReefArea, BCUarea and Protectedarea to the shapefile:
# I need to convert ReefArea, BCUarea and Protectedarea from m^2 to km^2
    # 1 - First check that areas are really expressed in square meters:
        (1/6*111) * (1/6)*111 # Area of a cell at the equator in square km
        # Ex cell 7401 is at the equator and has 19.12133% Reef Cover:
        maindata6[7401,]
        # It should be, in square km:
        (1/6*111) * (1/6)*111 * 0.1912133
        # That is: 
        60658966 * 1e-6
        # The difference is due to the position (not exactly at the equator), the length (not exactly 111 km). But what counts is the order of magnitude
        # It confirms that Areas are expressed in square meters can be converted to square km by multiplying them by 1e-6
    
    # 2 - Convert ReefArea, BCUarea and Protectedarea from m^2 to km^2
        maindata6$ReefArea <- maindata6$ReefArea * 1e-6
        maindata6$BCUarea <- maindata6$BCUarea * 1e-6
        maindata6$Protectedarea <- maindata6$Protectedarea * 1e-6
    
    
# Append complete dataframe to shape file
a@data <- maindata6

# Save and write the shapefile with all cells (including empty ones) if we ever need it...
save(a,file="cells12397.RData")
writeOGR(a,dsn=getwd(),layer="/cells12397",driver="ESRI Shapefile") #dsn is the folder it goes in, layer is the name of the shape file, will always use ESRI shapefile

# Remove empty cells
b <- a[-id.noReef,]
dim(b)

# Save and write the shapefile with reef cells only.
save(id.noReef,b,file="cells12292.RData")
writeOGR(b,dsn=getwd(),layer="/cells12292",driver="ESRI Shapefile")



