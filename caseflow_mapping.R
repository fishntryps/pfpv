
###########################################################
### this script maps malaria patient movement in Guyana ###
###########################################################

### movement is inferred via comparison of diagnosis localities/zones and patient responses about the location in which they stayed 2 weeks ago

### run time is seconds to minutes

### we start fresh by removing stored objects, load packages and set working/output directories in R 4.2.2

rm(list=ls())

library(tidyverse) 
library(RColorBrewer)    
library(patchwork)
library(ggrepel)
library(scatterpie)
library(MMWRweek)
library(zoo)
library(webshot2)
library(dplyr)
require(gridExtra)
library(grid)
library('ggnewscale')
library("FSA")

# Load dataframe of individual cases

main_path <-"C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana"
df_raw <- read.csv(paste(main_path,"data_clean.csv",sep="/"))

# here is the catalogue of geographic localities the ministry assigns case diagnosis and travel history responses to

geo_groups <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/COORDS_grouped_v4.tab", header=F)


### some tweaks are required such that locality names properly match the malaria case database

geo_groups$V1 <- gsub("_"," ", geo_groups$V1)
geo_groups$V2 <- gsub("_"," ", geo_groups$V2)
geo_groups$V1 <- gsub("\\.","", geo_groups$V1)
geo_groups$V1 <- tolower(geo_groups$V1)

### we will also be using centroids_via_COORDS_grouped_v4.txt created via...

# cent_epsg3786_lons <- c()
# cent_epsg3786_lats <- c()
# for (i in unique(geo_groups$V7)){
#   cent_epsg3786_lons <- c(cent_epsg3786_lons, mean(as.numeric(geo_groups$V2[geo_groups$V7==i])))
#   cent_epsg3786_lats <- c(cent_epsg3786_lats, mean(as.numeric(geo_groups$V3[geo_groups$V7==i])))
# }

# centroids <- as.data.frame(cbind(unique(geo_groups$V7),cent_epsg3786_lons, cent_epsg3786_lats))

# centroids[,2:3] convert to lon lat here https://mygeodata.cloud/cs2cs/

### more clean-up of geographic names

df_raw$Locality_infected <- gsub("\\.","", df_raw$Locality_infected)
df_raw$Locality_infected <- tolower(df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("b erbice","berbice", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("demarara$", "demerara", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("kamarangken settlemen.*", "kamarangkent settlement", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("vergenogen", "vergenoegen", df_raw$Locality_infected)

df_raw$Locality_sampled <- gsub("\\.","", df_raw$Locality_sampled)
df_raw$Locality_sampled <- tolower(df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("b erbice","berbice", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("demarara$", "demerara", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("kamarangken settlemen.*", "kamarangkent settlement", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("vergenogen", "vergenoegen", df_raw$Locality_sampled)

df_raw$geoinfect <- geo_groups$V7[match(df_raw$Locality_infected, geo_groups$V1)]
df_raw$geoinfect[df_raw$Locality_infected=="imported"] <- "imported"

df_raw$geosample <- geo_groups$V7[match(df_raw$Locality_sampled, geo_groups$V1)]
df_raw$geosample[is.na(df_raw$geosample)] <- ""
df_raw$geoinfect[is.na(df_raw$geoinfect)] <- ""

df_raw$District_infected <- df_raw$geoinfect
df_raw$District_sampled <- df_raw$geosample

# select only new cases, and passive detection

df <- df_raw %>% filter(Detection.method=="PASSIVE", Case == "NEW CASE")

# select only year of interest

df <- df %>% filter(Year_tested==2019)

# simplify ethnic groups

df <- df %>% mutate(EthnicGrouP = case_when(
  EthnicGrouP == "AFRO GUYANES" ~"Afro Guyanese",
  EthnicGrouP == "AMERINDIAN" ~ "Amerindian",
  EthnicGrouP == "EAST INDIAN" ~ "East Indian",
  EthnicGrouP == "MIXED" ~ "Mixed",
  TRUE ~ "Other" # This includes the few cases labelled as Chinese, European, and missing
))  %>% 
  # change order of Ethnic group manually 
  mutate(EthnicGrouP = factor(EthnicGrouP,levels=c("Other", 
                                                   "East Indian", "Afro Guyanese", "Amerindian", 
                                                   "Mixed")))
# Create a broader age group category
level.age <- c("<15", "15 - 55", "55+", "missing")

df <- df %>% 
  # add O for age less than one year in the AgeInYearMORE column
  mutate(AgeInYearMORE = case_when(
    !is.na(AgeInYearMORE) ~ AgeInYearMORE, # if older than 1 years old
    !is.na(LessThan1YEAR) ~ as.integer(0), # if younger than 1 year
    TRUE~ NA_integer_)) %>%  # if unknown
  # Create a new column with age group
  mutate(Age.group=case_when(
    AgeInYearMORE < 15 ~ level.age[1],
    AgeInYearMORE %in%c(15:54) ~ level.age[2],
    AgeInYearMORE >= 55 ~ level.age[3],
    TRUE~ "missing"
  ),
  Age.group = factor(Age.group, levels = level.age)
  )

### df is the basis for various epi exploration, e.g., proportion P. falciparum by age group

pvtemp <- df %>% filter(InfParas=="VIVAX") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group) %>% group_by(InfParas, Sex, Month_tested, Age.group) %>% summarise(n_cases=n())

pftemp <- df %>% filter(InfParas=="FALCIPARUM") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group) %>% group_by(InfParas, Sex, Month_tested, Age.group) %>% summarise(n_cases=n())

join <- pvtemp %>% left_join(y=pftemp, by=c("Month_tested", "Age.group", "Sex"))

join$prop_pf <- join$n_cases.y/(join$n_cases.x+join$n_cases.y)

ggplot(data=join[!is.na(join$prop_pf),] , aes(x=Age.group, y=prop_pf, fill=Sex)) + geom_hline(yintercept=0, lty="dashed", alpha=0) + geom_hline(yintercept=0.5, lty="dashed", alpha=1) +
  geom_boxplot(color="black", outlier.size = 1.5, size=0.7) + scale_fill_manual(values=c("indianred","navajowhite4")) + theme_classic() +
  theme(axis.title = element_blank()) + theme(axis.text = element_text(size=25, color="black")) + scale_y_continuous(breaks=c(0,0.25,0.5), labels=c("0%", "25%", "50%"))

wilcox.test(prop_pf ~ Age.group, data=join[join$Sex=="Female" & (join$Age.group=="<15" | join$Age.group=="15 - 55") ,])

kruskal.test(prop_pf ~ Age.group,
         data=join[join$Sex=="Male",])

dunnTest(prop_pf ~ Age.group,
         data=join[join$Sex=="Male",],
         method="bh")

dunnTest(prop_pf ~ Age.group,
         data=join[join$Sex=="Female",],
         method="bh")

### finer age group categorizatoin

level.age <- c("0-4", "5-9","10-14", "15-19","20-24", "25-29","30-34", "35-39","40-44", "45-49","50-54", "55-59","60+", "missing")

df <- df %>% 
  # add O for age less than one year in the AgeInYearMORE column
  mutate(AgeInYearMORE = case_when(
    !is.na(AgeInYearMORE) ~ AgeInYearMORE, # if older than 1 years old
    !is.na(LessThan1YEAR) ~ as.integer(0), # if younger than 1 year
    TRUE~ NA_integer_)) %>%  # if unknown
  # Create a new column with age group
  mutate(Age.group2=case_when(
    AgeInYearMORE %in% c(0:4) ~ level.age[1],
    AgeInYearMORE %in% c(5:9) ~ level.age[2],
    AgeInYearMORE %in% c(10:14) ~ level.age[3],
    AgeInYearMORE %in% c(15:19) ~ level.age[4],
    AgeInYearMORE %in% c(20:24) ~ level.age[5],
    AgeInYearMORE %in% c(25:29) ~ level.age[6],
    AgeInYearMORE %in% c(30:34) ~ level.age[7],
    AgeInYearMORE %in% c(35:39) ~ level.age[8],
    AgeInYearMORE %in% c(40:44) ~ level.age[9],
    AgeInYearMORE %in% c(45:49) ~ level.age[10],
    AgeInYearMORE %in% c(50:54) ~ level.age[11],
    AgeInYearMORE %in% c(55:59) ~ level.age[12],
    AgeInYearMORE >= 60 ~ level.age[13],
    TRUE~ "missing"
  ),
  Age.group = factor(Age.group, levels = level.age)
  )

pvtemp <- df %>% filter(InfParas=="VIVAX") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group2) %>% group_by(InfParas, Sex, Month_tested, Age.group2) %>% summarise(n_cases=n())

pftemp <- df %>% filter(InfParas=="FALCIPARUM") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group2) %>% group_by(InfParas, Sex, Month_tested, Age.group2) %>% summarise(n_cases=n())

join <- pvtemp %>% left_join(y=pftemp, by=c("Month_tested", "Age.group2", "Sex"))

join$prop_pf <- join$n_cases.y/(join$n_cases.x+join$n_cases.y)

ggplot(data=join[!is.na(join$prop_pf),] , aes(x=factor(Age.group2, levels=level.age), y=prop_pf, fill=Sex)) + geom_hline(yintercept=0, lty="dashed", alpha=0) + geom_hline(yintercept=1, lty="dashed", alpha=0) + geom_hline(yintercept=0.5, lty="dashed", alpha=1) +
  geom_boxplot(color="black", outlier.size = 1.5, size=0.7) + scale_fill_manual(values=c("indianred","navajowhite4")) + theme_classic() +
  theme(axis.title = element_blank()) + theme(axis.text = element_text(size=25, color="black")) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c("0%", "25%", "50%", "75%", "100%"))

#### back let's return to the case flow mapping objective...

# Flow data.frame per disctrict

df_district <- df %>% 
  # select year of interest

  filter(Year_tested == 2019) %>% 

  # select columns of interest

  dplyr::select("Region_sampled","Region_infected","District_sampled","District_infected","Date_tested","InfParas","Case","InfParas","Year_tested","vivax","falciparum","malariae","EpiWeek_tested") %>% 
  
  # remove the collected gametocyte only infections

  filter(InfParas!="GAMETOCYTES", Case=="NEW CASE") %>%
  
  # create a long dataframe for infections

  dplyr::select(-InfParas,-Case) %>% pivot_longer(cols=c("vivax","falciparum","malariae"), names_to = "Species") %>% 
  filter(!is.na(value)) %>% dplyr::select(-value) %>% 

  # specify if the case is local, foreign, unknown

  mutate(type = case_when((District_infected==District_sampled) & District_infected!="" ~ "Local",
                          (District_infected=="imported") ~ "imported",
                          (District_infected=="" | District_sampled=="") ~ "!Missing", 
                          (District_infected!=District_sampled & District_infected!="" & District_sampled!="" & District_infected!="imported") ~ "Outsourced"))

df_district$District_infected_gr <- gsub("$","_gr", df_district$District_infected)
df_district$District_sampled_gr <- gsub("$","_gr", df_district$District_sampled)

df_district <- df_district[df_district$type!="!Missing",]

### starting with P. falciparum, summarize outsourced case numbers based on the involved zones

targetspec="falciparum"

df_flow <- df_district %>%
  # sum yearly cases
  group_by(District_sampled,District_infected, Year_tested,Species,type, infectcol) %>% summarise(n_cases=n()) %>% 
  # keep only known local - imported, and positive flows 
  filter(n_cases > 0 & Species==targetspec) 

### here also some colors we will use for epidemiological zones

piesgroups <- c("Central_Coast","Chenapau","Cristinas_Border","East_of_GT","Greater_Annai","Greater_GT","Greater_Lethem","Head_Mazaruni","Head_Waini","Kaituma_and_Barima","Lower_Cuyuni","Lower_Essequibo","Lower_Mazaruni","Lower_Potaro","Mid_Berbice","Mid_Essequibo","Mid_Mazaruni_and_Issano_Rd","North_Delta","Pakaraima_South","Upper_Cuyuni","Upper_Demerara","Upper_Mazaruni","Upper_Rupununi","Waini","West_of_GT", "Head_Essequibo", "Mid_New", "Lower_New", "Greater_Apoteri")
piescols <- c("#ccd87c","#91e780","#df81c3","#4497d6","#c9779d","#e91046","#858fea","#b956cc","#d72124","#d57acf","#d1c025","#3c22de","#7f9edc","#9fdd32","#4dd795","#72f0d1","#71d577","#cd7c65","#d9bb78","#6212ec","#b565e9","#4dbddb","#0ec843","#df7c2b","#9fe86e", "#00ff00", "#00ffff", "#ff99ff", "#ffff66")
COLS <- as.data.frame(cbind(piesgroups,piescols))
df_district$infectcol <- COLS$piescols[match(df_district$District_infected, COLS$piesgroups)]
df_district$samplecol <- COLS$piescols[match(df_district$District_sampled, COLS$piesgroups)]

set.seed(123)  # set random generator state for the same output

geos <- read.table("centroids_via_COORDS_grouped_v4.txt", header=T)
geos$id <- rownames(geos)

dfOut <- df_flow %>% filter(type=="Outsourced")

q75Out <- quantile(dfOut$n_cases)[4]
sdOut <- sd(dfOut$n_cases)
dfOut$x <- geos$lon[match(dfOut$District_infected, geos$name)]
dfOut$y <- geos$lat[match(dfOut$District_infected, geos$name)]
dfOut$xend <- geos$lon[match(dfOut$District_sampled, geos$name)]
dfOut$yend <- geos$lat[match(dfOut$District_sampled, geos$name)]
dfOut$from<- geos$id[match(dfOut$District_infected, geos$name)]
dfOut$to<- geos$id[match(dfOut$District_sampled, geos$name)]
dfOut$weight <- dfOut$n_cases
dfOut$category <- dfOut$from
dfOut_nc_g50 <- dfOut %>% filter(n_cases >0 & n_cases > sdOut & n_cases <= 2*sdOut)  # keep only  positive flows
dfOut_nc_l50 <- dfOut %>% filter(n_cases >0 & n_cases <= sdOut)  # keep only  positive flows 
dfOut_nc_h50 <- dfOut %>% filter(n_cases >0 & n_cases > 2*sdOut)  # keep only  positive flows

dfLoc <- df_flow %>% filter(type=="Local")
dfLoc$x <- geos$lon[match(dfLoc$District_infected, geos$name)]
dfLoc$y <- geos$lat[match(dfLoc$District_infected, geos$name)]
dfLoc$xend <- geos$lon[match(dfLoc$District_sampled, geos$name)]
dfLoc$yend <- geos$lat[match(dfLoc$District_sampled, geos$name)]
dfLoc$from<- geos$id[match(dfLoc$District_infected, geos$name)]
dfLoc$to<- geos$id[match(dfLoc$District_sampled, geos$name)]
dfLoc$weight <- dfLoc$n_cases
dfLoc$category <- dfLoc$from

dfTot <- df_district  %>% filter(Species==targetspec) %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())
dfTotloc <- df_district  %>% filter(Species==targetspec & type=="Local") %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())
dfTotnotloc <- df_district  %>% filter(Species==targetspec & type=="Outsourced") %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())
dfTotimp <- df_district  %>% filter(Species==targetspec & type=="imported") %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())

dfTot$loc_cases <- dfTotloc$n_cases[match(dfTot$District_sampled, dfTotloc$District_sampled)]
dfTot$notloc_cases <- dfTotnotloc$n_cases[match(dfTot$District_sampled, dfTotnotloc$District_sampled)]
dfTot$imp_cases <- dfTotimp$n_cases[match(dfTot$District_sampled, dfTotimp$District_sampled)]
dfTot$loc_cases[is.na(dfTot$loc_cases)] <- 0
dfTot$notloc_cases[is.na(dfTot$notloc_cases)] <- 0
dfTot$imp_cases[is.na(dfTot$imp_cases)] <- 0
dfTot$x <- geos$lon[match(dfTot$District_sampled, geos$name)]
dfTot$y <- geos$lat[match(dfTot$District_sampled, geos$name)]

### categorize areas as endemic or non-endemic

endemics <- c("Central_Coast",
              "Chenapau",
              "Cristinas_Border",
              #"East_of_GT",
              "Greater_Annai",
              #"Greater_GT",
              "Greater_Lethem",
              "Head_Mazaruni",
              "Head_Waini",
              "Kaituma_and_Barima",
              "Lower_Cuyuni",
              "Lower_Essequibo",
              "Lower_Mazaruni",
              "Lower_Potaro",
              "Mid_Berbice",
              "Mid_Essequibo",
              "Mid_Mazaruni_and_Issano_Rd",
              "North_Delta",
              "Pakaraima_South",
              "Upper_Cuyuni",
              "Upper_Demerara",
              "Upper_Mazaruni",
              "Upper_Rupununi",
              "Waini",
              #"West_of_GT",
              "Head_Essequibo",
              "Mid_New",
              "Lower_New",
              "Greater_Apoteri")


# this function returns a data frame with interpolated points in between start and end point

interp_points_half <- function (data) {
  
  df <- data.frame(line_id=c(),x=c(),y=c())
  
  for (i in 1:nrow(data)) { 
    
    line <- data[i,]
    distAB <- ceiling(9*sqrt((line$xend - line$x)^2 + (line$yend - line$y)^2))
    # interpolate lats and longs in between the two
    longseq <- seq(
      as.numeric(line["x"]),
      as.numeric(line["xend"]),
      as.numeric((line["xend"] - line["x"])/distAB)
    )
    latseq <- seq(
      as.numeric(line["y"]),
      as.numeric(line["yend"]),
      as.numeric(line["yend"] - line["y"])/distAB
    )
    
    for (j in 1:distAB) {
      df <- rbind(df,data.frame(line_id=i,x=longseq[j],y=latseq[j],seg_num=j))
    }
  }
  
  df
}

# run data through function and correct direction of arrows

output <- interp_points_half(as.data.frame(dfOut)[,c(8:11)])
output$weight <- dfOut$weight[match(output$line_id, rownames(dfOut))]
output$District_sampled <- dfOut$District_sampled[match(output$line_id, rownames(dfOut))]
output$District_infected <- dfOut$District_infected[match(output$line_id, rownames(dfOut))]
output$angleInDegrees <- dfOut$angleInDegrees[match(output$line_id, rownames(dfOut))]
output_center <- output[output$seg_num==2,]

output$angleflip <- geos$lon[match(output$District_sampled, geos$name)]-geos$lon[match(output$District_infected, geos$name)] 
output$angleflip[output$angleflip > 0] <- 0
output$angleflip[output$angleflip < 0] <- 180 

output <- output %>% arrange(weight)

### colors to represent case flow numbers

output$segcolor <- "NA"
output$segcolor[output$weight>0 & output$weight<=2] <-  "lightyellow1"
output$segcolor[output$weight>2 & output$weight<=4] <-  brewer.pal(9, "YlOrRd") [2]
output$segcolor[output$weight>4 & output$weight<=6] <- brewer.pal(9, "YlOrRd") [3]
output$segcolor[output$weight>6 & output$weight<=8] <- brewer.pal(9, "YlOrRd") [4]
output$segcolor[output$weight>8 & output$weight<=12] <- brewer.pal(9, "YlOrRd") [5]
output$segcolor[output$weight>12 & output$weight<=24] <- brewer.pal(9, "YlOrRd") [6]
output$segcolor[output$weight>24 & output$weight<=48] <- brewer.pal(9, "YlOrRd") [7]
output$segcolor[output$weight>48 & output$weight<=72] <- brewer.pal(9, "YlOrRd") [8]
output$segcolor[output$weight>72 & output$weight<=100] <- brewer.pal(9, "YlOrRd") [9]
output$segcolor[output$weight>100] <-  "blue4"

plot(x=rep(1,10), y=1:10, col=c("lightyellow1", brewer.pal(9, "YlOrRd") [2],brewer.pal(9, "YlOrRd") [3],brewer.pal(9, "YlOrRd") [4], brewer.pal(9, "YlOrRd") [5],brewer.pal(9, "YlOrRd") [6],brewer.pal(9, "YlOrRd") [7],brewer.pal(9, "YlOrRd") [8],brewer.pal(9, "YlOrRd") [9], "blue4"), pch=15, cex=5)

### now we create the first map layer using GADM shape files for Guyana and bordering countries

source("C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools.R")
source("C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/fx_gadm_sp_loadCountries.R")

guy0.spldf <- gadm_sp_loadCountries("GUY",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
guy1.spldf <- gadm_sp_loadCountries("GUY",level = 1, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds1/")
guy2.spldf <- gadm_sp_loadCountries("GUY",level = 2, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds2/")
ven0.spldf <- gadm_sp_loadCountries("VEN",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
bra0.spldf <- gadm_sp_loadCountries("BRA",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
sur0.spldf <- gadm_sp_loadCountries("SUR",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")

maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))

p2 <- ggplot(geos) + ggtitle(paste("Plasmodium ", targetspec, " case flow", sep="")) + 
    
    geom_polygon(data=guy1.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    geom_polygon(data=bra0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    geom_polygon(data=ven0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    geom_polygon(data=sur0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    
    geom_path(aes(x=x,y=y,group=line_id), alpha=0.2, data=output[!output$District_sampled %in% endemics,], size=1, color=output$segcolor[!output$District_sampled %in% endemics]) +
    scale_color_gradient(low ="grey", high="black") + ggnewscale::new_scale_color() +
    geom_path(aes(x=x,y=y,group=line_id), alpha=0.2, data=output[output$District_sampled %in% endemics,], size=2.4, color=output$segcolor[output$District_sampled %in% endemics]) +
    #scale_size(range=c(0.1, 7.5)) +
     #scale_color_viridis_c("turbo") 
    #, size=3*sqrt(log10(output_center$weight[output_center$District_sampled %in% endemics]))
    
    geom_text(aes(x=x,y=y), label=">", alpha=1, data=output[output$District_sampled %in% endemics,], angle=output$angleInDegrees[output$District_sampled %in% endemics]+ output$angleflip[output$District_sampled %in% endemics],
              size=2.2, color=output$segcolor[output$District_sampled %in% endemics]) +
    geom_point(aes(x = lon, y = lat),           # draw nodes
              shape = 4, size = 1.6, data=geos[geos$name %in% unique(output$District_infected[output$weight > 0 & !output$District_infected %in% unique(dfTot$District_sampled)]),]) +
   
    coord_fixed(xlim = c(-61.5, -56.5), ylim = c(1.1, 8.35)) + maptheme 
  
for (i in dfTot$District_sampled[dfTot$District_sampled!=""]){
  p2 <- p2 + ggnewscale::new_scale_fill() +
    geom_scatterpie(aes(x=x, y=y, r=0.06*log10(n_cases)+0.02), cols=c("loc_cases", "notloc_cases", "imp_cases"),
                    data=dfTot[dfTot$District_sampled==i,], color=NA) +
    #scale_fill_manual(values=c("black",dfTot$infectcol[dfTot$District_infected==i])
    scale_fill_manual(values=c("black","grey", "cyan"))
}


Pf_p2 <- p2

### save underlying node and segment information for later correlation analysis to genetic data and to P. vivax

Pf_segs <- distinct(output[,5:7]) %>% arrange(weight)
colnames(Pf_segs) <- c("Pf_weight","Pf_District_sampled","Pf_District_infected")

Pf_locs <- dfTot[c("District_sampled", "n_cases","loc_cases")]
Pf_locs$weight <- Pf_locs$loc_cases/Pf_locs$n_cases


### now repeat mapping process for P. vivax

rm(list=setdiff(ls(), c("Pf_p2", "Pf_segs","Pf_locs")))

targetspec="vivax"

library(tidyverse) 
library(RColorBrewer)    
library(patchwork)
library(ggrepel)
library(scatterpie)
library(MMWRweek)
library(zoo)
library(webshot2)
library(dplyr)
require(gridExtra)
library(grid)
library('ggnewscale')

# Load dataframe of individual cases

main_path <-"C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana"
df_raw <- read.csv(paste(main_path,"data_clean.csv",sep="/"))

# here is the catalogue of geographic localities the ministry assigns case diagnosis and travel history responses to

geo_groups <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/COORDS_grouped_v4.tab", header=F)


### some tweaks are required such that locality names properly match the malaria case database

geo_groups$V1 <- gsub("_"," ", geo_groups$V1)
geo_groups$V2 <- gsub("_"," ", geo_groups$V2)
geo_groups$V1 <- gsub("\\.","", geo_groups$V1)
geo_groups$V1 <- tolower(geo_groups$V1)

### we will also be using centroids_via_COORDS_grouped_v4.txt created via...

# cent_epsg3786_lons <- c()
# cent_epsg3786_lats <- c()
# for (i in unique(geo_groups$V7)){
#   cent_epsg3786_lons <- c(cent_epsg3786_lons, mean(as.numeric(geo_groups$V2[geo_groups$V7==i])))
#   cent_epsg3786_lats <- c(cent_epsg3786_lats, mean(as.numeric(geo_groups$V3[geo_groups$V7==i])))
# }

# centroids <- as.data.frame(cbind(unique(geo_groups$V7),cent_epsg3786_lons, cent_epsg3786_lats))

# centroids[,2:3] convert to lon lat here https://mygeodata.cloud/cs2cs/

### more clean-up of geographic names

df_raw$Locality_infected <- gsub("\\.","", df_raw$Locality_infected)
df_raw$Locality_infected <- tolower(df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("b erbice","berbice", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("demarara$", "demerara", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("kamarangken settlemen.*", "kamarangkent settlement", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("vergenogen", "vergenoegen", df_raw$Locality_infected)

df_raw$Locality_sampled <- gsub("\\.","", df_raw$Locality_sampled)
df_raw$Locality_sampled <- tolower(df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("b erbice","berbice", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("demarara$", "demerara", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("kamarangken settlemen.*", "kamarangkent settlement", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("vergenogen", "vergenoegen", df_raw$Locality_sampled)

df_raw$geoinfect <- geo_groups$V7[match(df_raw$Locality_infected, geo_groups$V1)]
df_raw$geoinfect[df_raw$Locality_infected=="imported"] <- "imported"

df_raw$geosample <- geo_groups$V7[match(df_raw$Locality_sampled, geo_groups$V1)]
df_raw$geosample[is.na(df_raw$geosample)] <- ""
df_raw$geoinfect[is.na(df_raw$geoinfect)] <- ""

df_raw$District_infected <- df_raw$geoinfect
df_raw$District_sampled <- df_raw$geosample

# select only new cases, and passive detection

df <- df_raw %>% filter(Detection.method=="PASSIVE", Case == "NEW CASE")

# select only year of interest

df <- df %>% filter(Year_tested==2019)

# simplify ethnic groups

df <- df %>% mutate(EthnicGrouP = case_when(
  EthnicGrouP == "AFRO GUYANES" ~"Afro Guyanese",
  EthnicGrouP == "AMERINDIAN" ~ "Amerindian",
  EthnicGrouP == "EAST INDIAN" ~ "East Indian",
  EthnicGrouP == "MIXED" ~ "Mixed",
  TRUE ~ "Other" # This includes the few cases labelled as Chinese, European, and missing
))  %>% 
  # change order of Ethnic group manually 
  mutate(EthnicGrouP = factor(EthnicGrouP,levels=c("Other", 
                                                   "East Indian", "Afro Guyanese", "Amerindian", 
                                                   "Mixed")))
# Create a broader age group category
level.age <- c("<15", "15 - 55", "55+", "missing")

df <- df %>% 
  # add O for age less than one year in the AgeInYearMORE column
  mutate(AgeInYearMORE = case_when(
    !is.na(AgeInYearMORE) ~ AgeInYearMORE, # if older than 1 years old
    !is.na(LessThan1YEAR) ~ as.integer(0), # if younger than 1 year
    TRUE~ NA_integer_)) %>%  # if unknown
  # Create a new column with age group
  mutate(Age.group=case_when(
    AgeInYearMORE < 15 ~ level.age[1],
    AgeInYearMORE %in%c(15:54) ~ level.age[2],
    AgeInYearMORE >= 55 ~ level.age[3],
    TRUE~ "missing"
  ),
  Age.group = factor(Age.group, levels = level.age)
  )

### df is the basis for various epi exploration, e.g., proportion P. falciparum by age group

pvtemp <- df %>% filter(InfParas=="VIVAX") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group) %>% group_by(InfParas, Sex, Month_tested, Age.group) %>% summarise(n_cases=n())

pftemp <- df %>% filter(InfParas=="FALCIPARUM") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group) %>% group_by(InfParas, Sex, Month_tested, Age.group) %>% summarise(n_cases=n())

join <- pvtemp %>% left_join(y=pftemp, by=c("Month_tested", "Age.group", "Sex"))

join$prop_pf <- join$n_cases.y/(join$n_cases.x+join$n_cases.y)

ggplot(data=join[!is.na(join$prop_pf),] , aes(x=Age.group, y=prop_pf, fill=Sex)) + geom_hline(yintercept=0, lty="dashed", alpha=0) + geom_hline(yintercept=0.5, lty="dashed", alpha=1) +
  geom_boxplot(color="black", outlier.size = 1.5, size=0.7) + scale_fill_manual(values=c("indianred","navajowhite4")) + theme_classic() +
  theme(axis.title = element_blank()) + theme(axis.text = element_text(size=25, color="black")) + scale_y_continuous(breaks=c(0,0.25,0.5), labels=c("0%", "25%", "50%"))

wilcox.test(prop_pf ~ Age.group, data=join[join$Sex=="Female" & (join$Age.group=="<15" | join$Age.group=="15 - 55") ,])

kruskal.test(prop_pf ~ Age.group,
         data=join[join$Sex=="Male",])

dunnTest(prop_pf ~ Age.group,
         data=join[join$Sex=="Male",],
         method="bh")

dunnTest(prop_pf ~ Age.group,
         data=join[join$Sex=="Female",],
         method="bh")

### finer age group categorizatoin

level.age <- c("0-4", "5-9","10-14", "15-19","20-24", "25-29","30-34", "35-39","40-44", "45-49","50-54", "55-59","60+", "missing")

df <- df %>% 
  # add O for age less than one year in the AgeInYearMORE column
  mutate(AgeInYearMORE = case_when(
    !is.na(AgeInYearMORE) ~ AgeInYearMORE, # if older than 1 years old
    !is.na(LessThan1YEAR) ~ as.integer(0), # if younger than 1 year
    TRUE~ NA_integer_)) %>%  # if unknown
  # Create a new column with age group
  mutate(Age.group2=case_when(
    AgeInYearMORE %in% c(0:4) ~ level.age[1],
    AgeInYearMORE %in% c(5:9) ~ level.age[2],
    AgeInYearMORE %in% c(10:14) ~ level.age[3],
    AgeInYearMORE %in% c(15:19) ~ level.age[4],
    AgeInYearMORE %in% c(20:24) ~ level.age[5],
    AgeInYearMORE %in% c(25:29) ~ level.age[6],
    AgeInYearMORE %in% c(30:34) ~ level.age[7],
    AgeInYearMORE %in% c(35:39) ~ level.age[8],
    AgeInYearMORE %in% c(40:44) ~ level.age[9],
    AgeInYearMORE %in% c(45:49) ~ level.age[10],
    AgeInYearMORE %in% c(50:54) ~ level.age[11],
    AgeInYearMORE %in% c(55:59) ~ level.age[12],
    AgeInYearMORE >= 60 ~ level.age[13],
    TRUE~ "missing"
  ),
  Age.group = factor(Age.group, levels = level.age)
  )

pvtemp <- df %>% filter(InfParas=="VIVAX") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group2) %>% group_by(InfParas, Sex, Month_tested, Age.group2) %>% summarise(n_cases=n())

pftemp <- df %>% filter(InfParas=="FALCIPARUM") %>% dplyr::select(InfParas, Sex, Month_tested, Age.group2) %>% group_by(InfParas, Sex, Month_tested, Age.group2) %>% summarise(n_cases=n())

join <- pvtemp %>% left_join(y=pftemp, by=c("Month_tested", "Age.group2", "Sex"))

join$prop_pf <- join$n_cases.y/(join$n_cases.x+join$n_cases.y)

ggplot(data=join[!is.na(join$prop_pf),] , aes(x=factor(Age.group2, levels=level.age), y=prop_pf, fill=Sex)) + geom_hline(yintercept=0, lty="dashed", alpha=0) + geom_hline(yintercept=1, lty="dashed", alpha=0) + geom_hline(yintercept=0.5, lty="dashed", alpha=1) +
  geom_boxplot(color="black", outlier.size = 1.5, size=0.7) + scale_fill_manual(values=c("indianred","navajowhite4")) + theme_classic() +
  theme(axis.title = element_blank()) + theme(axis.text = element_text(size=25, color="black")) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c("0%", "25%", "50%", "75%", "100%"))

#### back let's return to the case flow mapping objective...

# Flow data.frame per disctrict

df_district <- df %>% 
  # select year of interest

  filter(Year_tested == 2019) %>% 

  # select columns of interest

  dplyr::select("Region_sampled","Region_infected","District_sampled","District_infected","Date_tested","InfParas","Case","InfParas","Year_tested","vivax","falciparum","malariae","EpiWeek_tested") %>% 
  
  # remove the collected gametocyte only infections

  filter(InfParas!="GAMETOCYTES", Case=="NEW CASE") %>%
  
  # create a long dataframe for infections

  dplyr::select(-InfParas,-Case) %>% pivot_longer(cols=c("vivax","falciparum","malariae"), names_to = "Species") %>% 
  filter(!is.na(value)) %>% dplyr::select(-value) %>% 

  # specify if the case is local, foreign, unknown

  mutate(type = case_when((District_infected==District_sampled) & District_infected!="" ~ "Local",
                          (District_infected=="imported") ~ "imported",
                          (District_infected=="" | District_sampled=="") ~ "!Missing", 
                          (District_infected!=District_sampled & District_infected!="" & District_sampled!="" & District_infected!="imported") ~ "Outsourced"))

df_district$District_infected_gr <- gsub("$","_gr", df_district$District_infected)
df_district$District_sampled_gr <- gsub("$","_gr", df_district$District_sampled)

df_district <- df_district[df_district$type!="!Missing",]

### summarize outsourced case numbers based on the involved zones

df_flow <- df_district %>%
  # sum yearly cases
  group_by(District_sampled,District_infected, Year_tested,Species,type, infectcol) %>% summarise(n_cases=n()) %>% 
  # keep only known local - imported, and positive flows 
  filter(n_cases > 0 & Species==targetspec) 

### here also some colors we will use for epidemiological zones

piesgroups <- c("Central_Coast","Chenapau","Cristinas_Border","East_of_GT","Greater_Annai","Greater_GT","Greater_Lethem","Head_Mazaruni","Head_Waini","Kaituma_and_Barima","Lower_Cuyuni","Lower_Essequibo","Lower_Mazaruni","Lower_Potaro","Mid_Berbice","Mid_Essequibo","Mid_Mazaruni_and_Issano_Rd","North_Delta","Pakaraima_South","Upper_Cuyuni","Upper_Demerara","Upper_Mazaruni","Upper_Rupununi","Waini","West_of_GT", "Head_Essequibo", "Mid_New", "Lower_New", "Greater_Apoteri")
piescols <- c("#ccd87c","#91e780","#df81c3","#4497d6","#c9779d","#e91046","#858fea","#b956cc","#d72124","#d57acf","#d1c025","#3c22de","#7f9edc","#9fdd32","#4dd795","#72f0d1","#71d577","#cd7c65","#d9bb78","#6212ec","#b565e9","#4dbddb","#0ec843","#df7c2b","#9fe86e", "#00ff00", "#00ffff", "#ff99ff", "#ffff66")
COLS <- as.data.frame(cbind(piesgroups,piescols))
df_district$infectcol <- COLS$piescols[match(df_district$District_infected, COLS$piesgroups)]
df_district$samplecol <- COLS$piescols[match(df_district$District_sampled, COLS$piesgroups)]

set.seed(123)  # set random generator state for the same output

geos <- read.table("centroids_via_COORDS_grouped_v4.txt", header=T)
geos$id <- rownames(geos)

dfOut <- df_flow %>% filter(type=="Outsourced")

q75Out <- quantile(dfOut$n_cases)[4]
sdOut <- sd(dfOut$n_cases)
dfOut$x <- geos$lon[match(dfOut$District_infected, geos$name)]
dfOut$y <- geos$lat[match(dfOut$District_infected, geos$name)]
dfOut$xend <- geos$lon[match(dfOut$District_sampled, geos$name)]
dfOut$yend <- geos$lat[match(dfOut$District_sampled, geos$name)]
dfOut$from<- geos$id[match(dfOut$District_infected, geos$name)]
dfOut$to<- geos$id[match(dfOut$District_sampled, geos$name)]
dfOut$weight <- dfOut$n_cases
dfOut$category <- dfOut$from
dfOut_nc_g50 <- dfOut %>% filter(n_cases >0 & n_cases > sdOut & n_cases <= 2*sdOut)  # keep only  positive flows
dfOut_nc_l50 <- dfOut %>% filter(n_cases >0 & n_cases <= sdOut)  # keep only  positive flows 
dfOut_nc_h50 <- dfOut %>% filter(n_cases >0 & n_cases > 2*sdOut)  # keep only  positive flows

dfLoc <- df_flow %>% filter(type=="Local")
dfLoc$x <- geos$lon[match(dfLoc$District_infected, geos$name)]
dfLoc$y <- geos$lat[match(dfLoc$District_infected, geos$name)]
dfLoc$xend <- geos$lon[match(dfLoc$District_sampled, geos$name)]
dfLoc$yend <- geos$lat[match(dfLoc$District_sampled, geos$name)]
dfLoc$from<- geos$id[match(dfLoc$District_infected, geos$name)]
dfLoc$to<- geos$id[match(dfLoc$District_sampled, geos$name)]
dfLoc$weight <- dfLoc$n_cases
dfLoc$category <- dfLoc$from

dfTot <- df_district  %>% filter(Species==targetspec) %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())
dfTotloc <- df_district  %>% filter(Species==targetspec & type=="Local") %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())
dfTotnotloc <- df_district  %>% filter(Species==targetspec & type=="Outsourced") %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())
dfTotimp <- df_district  %>% filter(Species==targetspec & type=="imported") %>% group_by(District_sampled, samplecol) %>% summarise(n_cases=n())

dfTot$loc_cases <- dfTotloc$n_cases[match(dfTot$District_sampled, dfTotloc$District_sampled)]
dfTot$notloc_cases <- dfTotnotloc$n_cases[match(dfTot$District_sampled, dfTotnotloc$District_sampled)]
dfTot$imp_cases <- dfTotimp$n_cases[match(dfTot$District_sampled, dfTotimp$District_sampled)]
dfTot$loc_cases[is.na(dfTot$loc_cases)] <- 0
dfTot$notloc_cases[is.na(dfTot$notloc_cases)] <- 0
dfTot$imp_cases[is.na(dfTot$imp_cases)] <- 0
dfTot$x <- geos$lon[match(dfTot$District_sampled, geos$name)]
dfTot$y <- geos$lat[match(dfTot$District_sampled, geos$name)]

### categorize areas as endemic or non-endemic

endemics <- c("Central_Coast",
              "Chenapau",
              "Cristinas_Border",
              #"East_of_GT",
              "Greater_Annai",
              #"Greater_GT",
              "Greater_Lethem",
              "Head_Mazaruni",
              "Head_Waini",
              "Kaituma_and_Barima",
              "Lower_Cuyuni",
              "Lower_Essequibo",
              "Lower_Mazaruni",
              "Lower_Potaro",
              "Mid_Berbice",
              "Mid_Essequibo",
              "Mid_Mazaruni_and_Issano_Rd",
              "North_Delta",
              "Pakaraima_South",
              "Upper_Cuyuni",
              "Upper_Demerara",
              "Upper_Mazaruni",
              "Upper_Rupununi",
              "Waini",
              #"West_of_GT",
              "Head_Essequibo",
              "Mid_New",
              "Lower_New",
              "Greater_Apoteri")


# this function returns a data frame with interpolated points in between start and end point

interp_points_half <- function (data) {
  
  df <- data.frame(line_id=c(),x=c(),y=c())
  
  for (i in 1:nrow(data)) { 
    
    line <- data[i,]
    distAB <- ceiling(9*sqrt((line$xend - line$x)^2 + (line$yend - line$y)^2))
    # interpolate lats and longs in between the two
    longseq <- seq(
      as.numeric(line["x"]),
      as.numeric(line["xend"]),
      as.numeric((line["xend"] - line["x"])/distAB)
    )
    latseq <- seq(
      as.numeric(line["y"]),
      as.numeric(line["yend"]),
      as.numeric(line["yend"] - line["y"])/distAB
    )
    
    for (j in 1:distAB) {
      df <- rbind(df,data.frame(line_id=i,x=longseq[j],y=latseq[j],seg_num=j))
    }
  }
  
  df
}

# run data through function and correct direction of arrows

output <- interp_points_half(as.data.frame(dfOut)[,c(8:11)])
output$weight <- dfOut$weight[match(output$line_id, rownames(dfOut))]
output$District_sampled <- dfOut$District_sampled[match(output$line_id, rownames(dfOut))]
output$District_infected <- dfOut$District_infected[match(output$line_id, rownames(dfOut))]
output$angleInDegrees <- dfOut$angleInDegrees[match(output$line_id, rownames(dfOut))]
output_center <- output[output$seg_num==2,]

output$angleflip <- geos$lon[match(output$District_sampled, geos$name)]-geos$lon[match(output$District_infected, geos$name)] 
output$angleflip[output$angleflip > 0] <- 0
output$angleflip[output$angleflip < 0] <- 180 

output <- output %>% arrange(weight)

### colors to represent case flow numbers

output$segcolor <- "NA"
output$segcolor[output$weight>0 & output$weight<=2] <-  "lightyellow1"
output$segcolor[output$weight>2 & output$weight<=4] <-  brewer.pal(9, "YlOrRd") [2]
output$segcolor[output$weight>4 & output$weight<=6] <- brewer.pal(9, "YlOrRd") [3]
output$segcolor[output$weight>6 & output$weight<=8] <- brewer.pal(9, "YlOrRd") [4]
output$segcolor[output$weight>8 & output$weight<=12] <- brewer.pal(9, "YlOrRd") [5]
output$segcolor[output$weight>12 & output$weight<=24] <- brewer.pal(9, "YlOrRd") [6]
output$segcolor[output$weight>24 & output$weight<=48] <- brewer.pal(9, "YlOrRd") [7]
output$segcolor[output$weight>48 & output$weight<=72] <- brewer.pal(9, "YlOrRd") [8]
output$segcolor[output$weight>72 & output$weight<=100] <- brewer.pal(9, "YlOrRd") [9]
output$segcolor[output$weight>100] <-  "blue4"

plot(x=rep(1,10), y=1:10, col=c("lightyellow1", brewer.pal(9, "YlOrRd") [2],brewer.pal(9, "YlOrRd") [3],brewer.pal(9, "YlOrRd") [4], brewer.pal(9, "YlOrRd") [5],brewer.pal(9, "YlOrRd") [6],brewer.pal(9, "YlOrRd") [7],brewer.pal(9, "YlOrRd") [8],brewer.pal(9, "YlOrRd") [9], "blue4"), pch=15, cex=5)

### now we create the first map layer using GADM shape files for Guyana and bordering countries

source("C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools.R")
source("C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/fx_gadm_sp_loadCountries.R")

guy0.spldf <- gadm_sp_loadCountries("GUY",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
guy1.spldf <- gadm_sp_loadCountries("GUY",level = 1, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds1/")
guy2.spldf <- gadm_sp_loadCountries("GUY",level = 2, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds2/")
ven0.spldf <- gadm_sp_loadCountries("VEN",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
bra0.spldf <- gadm_sp_loadCountries("BRA",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
sur0.spldf <- gadm_sp_loadCountries("SUR",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")

maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))

p2 <- ggplot(geos) + ggtitle(paste("Plasmodium ", targetspec, " case flow", sep="")) + 
    
    geom_polygon(data=guy1.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    geom_polygon(data=bra0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    geom_polygon(data=ven0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    geom_polygon(data=sur0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
    
    geom_path(aes(x=x,y=y,group=line_id), alpha=0.2, data=output[!output$District_sampled %in% endemics,], size=1, color=output$segcolor[!output$District_sampled %in% endemics]) +
    scale_color_gradient(low ="grey", high="black") + ggnewscale::new_scale_color() +
    geom_path(aes(x=x,y=y,group=line_id), alpha=0.2, data=output[output$District_sampled %in% endemics,], size=2.4, color=output$segcolor[output$District_sampled %in% endemics]) +
    #scale_size(range=c(0.1, 7.5)) +
     #scale_color_viridis_c("turbo") 
    #, size=3*sqrt(log10(output_center$weight[output_center$District_sampled %in% endemics]))
    
    geom_text(aes(x=x,y=y), label=">", alpha=1, data=output[output$District_sampled %in% endemics,], angle=output$angleInDegrees[output$District_sampled %in% endemics]+ output$angleflip[output$District_sampled %in% endemics],
              size=2.2, color=output$segcolor[output$District_sampled %in% endemics]) +
    geom_point(aes(x = lon, y = lat),           # draw nodes
              shape = 4, size = 1.6, data=geos[geos$name %in% unique(output$District_infected[output$weight > 0 & !output$District_infected %in% unique(dfTot$District_sampled)]),]) +
   
    coord_fixed(xlim = c(-61.5, -56.5), ylim = c(1.1, 8.35)) + maptheme 
  
for (i in dfTot$District_sampled[dfTot$District_sampled!=""]){
  p2 <- p2 + ggnewscale::new_scale_fill() +
    geom_scatterpie(aes(x=x, y=y, r=0.06*log10(n_cases)+0.02), cols=c("loc_cases", "notloc_cases", "imp_cases"),
                    data=dfTot[dfTot$District_sampled==i,], color=NA) +
    #scale_fill_manual(values=c("black",dfTot$infectcol[dfTot$District_infected==i])
    scale_fill_manual(values=c("black","grey", "cyan"))
}


Pv_p2 <- p2

Pv_segs <- distinct(output[,5:7]) %>% arrange(District_infected)
colnames(Pv_segs) <- c("Pv_weight","Pv_District_sampled","Pv_District_infected")


### various correlation analysis can now be applied to the data underlying Pf_p2 and Pv_p2 maps

Pv_segs$Pv_link <- paste(Pv_segs$Pv_District_infected,Pv_segs$Pv_District_sampled, sep="\nto ")
Pf_segs$Pf_link <- paste(Pf_segs$Pf_District_infected,Pf_segs$Pf_District_sampled, sep="\nto ")

combo_links <- Pv_segs[,c(1,4)] %>% left_join(y = Pf_segs[,c(1,4)], by = c("Pv_link" = "Pf_link"))

Pv_locs <- dfTot[c("District_sampled", "n_cases","loc_cases")]
Pv_locs$weight <- Pv_locs$loc_cases/Pv_locs$n_cases
combo_locs <- Pv_locs[,c(1,4)] %>% left_join(y = Pf_locs[,c(1,4)], by = c("District_sampled" = "District_sampled"))
combo_locs <- combo_locs[combo_locs$District_sampled!="",]

par(mfrow=c(1,1))

plot(combo_locs$weight.y, combo_locs$weight.x, pch=19, col="blue", xlab="\nP. falciparum: % sampled cases representing local infection", ylab="P. vivax: % sampled cases representing local infection\n", xlim=c(0,1), ylim=c(0,1), cex=3, las=1)

par(new=T)

ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(data=combo_locs[combo_locs$weight.y<combo_locs$weight.x,], aes(x=100*weight.y, y=100*weight.x), color="green", size=5) +
  geom_point(data=combo_locs[combo_locs$weight.y>combo_locs$weight.x,],aes(x=100*weight.y, y=100*weight.x), color="blue", size=5) +
  labs(x= "\nP. falciparum: % sampled cases representing local infection", y="P. vivax: % sampled cases representing local infection\n") +
  theme_classic() + theme(axis.title = element_text(size=15)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(data=combo_locs, aes(x=100*weight.y, y=100*weight.x, label=District_sampled), size=3, segment.color = NA)

par(mfrow=c(1,1))

plot(sqrt(combo_links$Pf_weight), sqrt(combo_links$Pv_weight), xaxt='n', yaxt='n', xlab="\nP. falciparum case flow", ylab="P. vivax case flow\n", pch=19, col="blue", cex=3)
axis(1, at=sqrt(c(10, 50, 100, 200)), labels=(c(10, 50, 100, 200)))
axis(2, at=sqrt(c(10, 50, 100, 200, 500)), labels=(c(10, 50, 100, 200, 500)), las=1)
cor.test(combo_links$Pf_weight,combo_links$Pv_weight)
text(sqrt(combo_links$Pf_weight), sqrt(combo_links$Pv_weight), labels=combo_links$Pv_link, cex=0.7, xaxt='n', yaxt='n', xlab="P. falciparum case flow", ylab="P. vivax case flow")

combo_links$Pf_weight2 <- combo_links$Pf_weight
combo_links$Pf_weight2[is.na(combo_links$Pf_weight==T)] <- 0
combo_links$Pv_weight2 <- combo_links$Pv_weight
combo_links$Pv_weight2[is.na(combo_links$Pv_weight2)==T] <- 0

without_flow_to_cap <- ggplot(data=combo_links[!grepl("Greater_GT|West_of_GT|Lower_Essequibo", combo_links$Pv_link),], aes(label=Pv_link, x=Pf_weight2/sum(combo_links$Pf_weight2, na.rm=T), y=Pv_weight2/sum(combo_links$Pv_weight2, na.rm=T))) + geom_smooth(method="lm") + geom_point(color="black", size=5, alpha=0.2) +
  labs(x= "\nP. falciparum: case flow as % of total case flow", y="P. vivax: case flow as % of total case flow\n") +
  theme_classic() + theme(axis.title = element_text(size=15)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(size=2, segment.color = NA, max.overlaps=15) 

with_flow_to_cap <- ggplot(data=combo_links, aes(label=Pv_link, x=Pf_weight2/sum(combo_links$Pf_weight2, na.rm=T), y=Pv_weight2/sum(combo_links$Pv_weight2, na.rm=T))) + geom_smooth(method="lm") + geom_point(color="black", size=5, alpha=0.2) +
  labs(x= "\nP. falciparum: case flow as % of total case flow", y="P. vivax: case flow as % of total case flow\n") +
  theme_classic() + theme(axis.title = element_text(size=15)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(size=2, segment.color = NA, max.overlaps=15) 

grid.arrange(without_flow_to_cap, with_flow_to_cap, ncol=1)

y=combo_links$Pv_weight2[!grepl("Greater_GT|West_of_GT|Lower_Essequibo", combo_links$Pv_link)]/sum(combo_links$Pv_weight2[!grepl("Greater_GT|West_of_GT|Lower_Essequibo", combo_links$Pv_link)], na.rm=T) 
x=combo_links$Pf_weight2[!grepl("Greater_GT|West_of_GT|Lower_Essequibo", combo_links$Pv_link)]/sum(combo_links$Pf_weight2[!grepl("Greater_GT|West_of_GT|Lower_Essequibo", combo_links$Pv_link)], na.rm=T) 
cor.test(y~x)

y=combo_links$Pv_weight2/sum(combo_links$Pv_weight2, na.rm=T) 
x=combo_links$Pf_weight2/sum(combo_links$Pf_weight2, na.rm=T) 

cor.test(y~x)

PfGTloc <- Pf_locs$loc_cases[Pf_locs$District_sampled=="Greater_GT"]
PfGTnonloc <- Pf_locs$n_cases[Pf_locs$District_sampled=="Greater_GT"]-Pf_locs$loc_cases[Pf_locs$District_sampled=="Greater_GT"]                  
                  
PvGTloc <- Pv_locs$loc_cases[Pv_locs$District_sampled=="Greater_GT"]
PvGTnonloc <- Pv_locs$n_cases[Pv_locs$District_sampled=="Greater_GT"]-Pv_locs$loc_cases[Pv_locs$District_sampled=="Greater_GT"]                  

GT_chi <- matrix(c(PfGTloc,
PfGTnonloc,
PvGTloc,
PvGTnonloc), nrow=2, byrow=T)

colnames(GT_chi) <- c("loc_n","nonloc_n")
rownames(GT_chi) <- c("falciparum","vivax")

chisq.test(GT_chi)
45/867
178/1611

mode <- c(1,2)[2]
thresh=c(">0.5","p99")[mode]
thresh_Pv=c("g50", "g11.5432")[mode]
thresh_Pf=c("g50", "g69.35646")[mode]
ncomp=50
version=c("v3","v4")[2]

### consider omitting these localities given less than 10 gen samples

# East_of_GT_falciparum 
# 1 
# Head_Mazaruni_vivax 
# 1 
# Pakaraima_South_vivax 
# 1 
# Mid_Berbice_vivax 
# 2 
# Central_Coast_falciparum 
# 3 
# East_of_GT_vivax 
# 3 
# Greater_GT_falciparum 
# 3 
# Greater_Lethem_vivax 
# 3 
# Head_Mazaruni_falciparum 
# 3 
# Upper_Demerara_vivax 
# 3 
# West_of_GT_falciparum 
# 3 
# Head_Waini_falciparum 
# 4 
# Head_Waini_vivax 
# 4 
# Upper_Rupununi_vivax 
# 4 
# Central_Coast_vivax 
# 5 
# North_Delta_falciparum 
# 5 
# North_Delta_vivax 
# 5 
# West_of_GT_vivax 
# 8

Pf_segs$Pf_epi_bi <- paste(pmin(Pf_segs$Pf_District_infected,Pf_segs$Pf_District_sampled), pmax(Pf_segs$Pf_District_infected,Pf_segs$Pf_District_sampled), sep="_to_")
Pf_gen <- read.table(paste("Pf_perc_", thresh_Pf, "_amoung_groups_",version,"_with_distances_and_indcount_and_cccount.txt", sep=""))
colnames(Pf_gen) <- gsub("^","Pf_", colnames(Pf_gen))
Pf_gen$Pf_gen_link <- paste(pmin(Pf_gen$Pf_V1, Pf_gen$Pf_V2), pmax(Pf_gen$Pf_V1, Pf_gen$Pf_V2), sep="_to_")
Pf_gen <- Pf_gen %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo|East_of_GT", Pf_gen_link))
Pf_bi <- Pf_segs[,c("Pf_epi_bi","Pf_weight")] %>% group_by(Pf_epi_bi) %>% summarise(Pf_sum=sum(Pf_weight))
Pf_bi <- Pf_bi %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo|East_of_GT", Pf_epi_bi))

Pv_segs$Pv_epi_bi <- paste(pmin(Pv_segs$Pv_District_infected,Pv_segs$Pv_District_sampled), pmax(Pv_segs$Pv_District_infected,Pv_segs$Pv_District_sampled), sep="_to_")
Pv_gen <- read.table(paste("Pv_perc_", thresh_Pv, "_amoung_groups_",version,"_with_distances_and_indcount_and_cccount.txt", sep=""))
colnames(Pv_gen) <- gsub("^","Pv_", colnames(Pv_gen))
Pv_gen$Pv_gen_link <- paste(pmin(Pv_gen$Pv_V1, Pv_gen$Pv_V2), pmax(Pv_gen$Pv_V1, Pv_gen$Pv_V2), sep="_to_")
Pv_gen <- Pv_gen %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo|Head_Mazaruni|Pakaraima_South", Pv_gen_link))
Pv_bi <- Pv_segs[,c("Pv_epi_bi","Pv_weight")] %>% group_by(Pv_epi_bi) %>% summarise(Pv_sum=sum(Pv_weight))
Pv_bi <- Pv_bi %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo|Head_Mazaruni|Pakaraima_South", Pv_epi_bi))

epi_vs_gen <- Pv_bi %>% left_join(y = Pf_bi, by = c("Pv_epi_bi" = "Pf_epi_bi")) %>% left_join(y = Pv_gen, by = c("Pv_epi_bi" = "Pv_gen_link")) %>% left_join(y = Pf_gen, by = c("Pv_epi_bi" = "Pf_gen_link")) %>% 
  filter(is.na(Pv_V4)==F & is.na(Pf_V4)==F & is.na(Pv_sum)==F & is.na(Pf_sum)==F) %>% filter(Pv_V4>ncomp & Pf_V4>ncomp) #####%>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo", Pv_epi_bi))

epi_vs_gen_Pv <- Pv_bi %>% left_join(y = Pv_gen, by = c("Pv_epi_bi" = "Pv_gen_link")) %>% filter(Pv_V4>ncomp) #####%>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo", Pv_epi_bi))
epi_vs_gen_Pf <- Pf_bi %>% left_join(y = Pf_gen, by = c("Pf_epi_bi" = "Pf_gen_link")) %>%  filter(Pf_V4>ncomp) #######%>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo", Pf_epi_bi))
colnames(epi_vs_gen_Pv) <- gsub("^Pv_", "", colnames(epi_vs_gen_Pv))
colnames(epi_vs_gen_Pf) <- gsub("^Pf_", "", colnames(epi_vs_gen_Pv))
epi_vs_gen_Pf$spec <- "falciparum"
epi_vs_gen_Pv$spec <- "vivax"
combo_indep <- rbind(epi_vs_gen_Pv,epi_vs_gen_Pf)
combo_indep$gen_as_frac_max <- "NA"
combo_indep$gen_as_frac_max[combo_indep$spec=="vivax"] <- combo_indep$V5[combo_indep$spec=="vivax"]/max(combo_indep$V5[combo_indep$spec=="vivax"])
combo_indep$gen_as_frac_max[combo_indep$spec=="falciparum"] <- combo_indep$V5[combo_indep$spec=="falciparum"]/max(combo_indep$V5[combo_indep$spec=="falciparum"])
combo_indep$totspecsum <- "NA"
combo_indep$totspecsum[combo_indep$spec=="falciparum"] <- sum(combo_indep$sum[combo_indep$spec=="falciparum"])
combo_indep$totspecsum[combo_indep$spec=="vivax"] <- sum(combo_indep$sum[combo_indep$spec=="vivax"])
combo_indep$totunspecsum <- sum(combo_indep$sum)
combo_indep$gen_as_frac_max <- as.numeric(combo_indep$gen_as_frac_max)
combo_indep$totspecsum <- as.numeric(combo_indep$totspecsum)
combo_indep$V4 <- as.numeric(combo_indep$V4)

all_seg_names <- table(c(unique(combo_indep$epi_bi[combo_indep$spec=="vivax"]), unique(combo_indep$epi_bi[combo_indep$spec=="falciparum"])))
combo_indep_shared_cf_detection <- combo_indep[combo_indep$epi_bi %in% names(all_seg_names[all_seg_names==2]),]

Pf_segs[grepl("Rupununi", Pf_segs$Pf_epi_bi),]
     
######## pooling case flow between species i.e. totunspecsum

### using separate epi-positive site pairs for each species

ggplot(aes(x=100*sum/totunspecsum, y=100*gen_as_frac_max, color=spec, size=V4), data=combo_indep) + geom_smooth(method="lm") + geom_point(alpha=0.7) +
    labs(x="\nCase flow (as % of total)", y=paste("Frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), size="# pairwise comparisons (gen)", color="Species") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + theme(axis.title = element_text(size=20)) +
  coord_cartesian(ylim = c(0,100)) + scale_y_continuous(trans='sqrt') + scale_size_continuous(range=c(3,12)) + facet_grid(.~spec) + theme_minimal() + theme(axis.text = element_text(size=12)) + theme(axis.title = element_text(size=20)) +
  geom_text_repel(data=combo_indep, aes(x=100*sum/totunspecsum, y=100*gen_as_frac_max, label=epi_bi), size=3, segment.color = NA)

y=combo_indep$sum[combo_indep$spec=="falciparum"]/combo_indep$totunspecsum[combo_indep$spec=="falciparum"]
x=combo_indep$gen_as_frac_max[combo_indep$spec=="falciparum"]
summary(lm(y~x))
y=combo_indep$sum[combo_indep$spec=="vivax"]/combo_indep$totunspecsum[combo_indep$spec=="vivax"]
x=combo_indep$gen_as_frac_max[combo_indep$spec=="vivax"]
summary(lm(y~x))

### NOT pooling case flow between species i.e. totspecsum

### using separate epi-positive site pairs for each species
ggplot(aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, color=spec, size=V4), data=combo_indep) + geom_smooth(method="lm") + geom_point(alpha=0.7) +
  labs(x="\nCase flow (as % of total)", y=paste("Frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), size="# pairwise comparisons (gen)", color="Species") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + theme(axis.title = element_text(size=20)) +
  coord_cartesian(ylim = c(0,100)) + scale_y_continuous(trans='sqrt') + scale_size_continuous(range=c(3,12)) + facet_grid(.~spec) + theme_minimal() + theme(axis.text = element_text(size=12)) + theme(axis.title = element_text(size=20)) +
  geom_text_repel(data=combo_indep, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, label=epi_bi), size=3, segment.color = NA)

combo_indep_fal <- combo_indep[combo_indep$spec=="falciparum",]
combo_indep_viv <- combo_indep[combo_indep$spec=="vivax",] 

cbind(seg_indices,1:length(seg_indices))
seg_indices <- unique(c(combo_indep_fal$epi_bi,combo_indep_viv$epi_bi))
seg_indices <- as.data.frame(cbind(seg_indices,1:length(seg_indices)))
colnames(seg_indices) <- c("seg_name","seg_index")

combo_indep_fal$num <- seg_indices$seg_index[match(combo_indep_fal$epi_bi, seg_indices$seg_name)]
combo_indep_viv$num <- seg_indices$seg_index[match(combo_indep_viv$epi_bi, seg_indices$seg_name)]

combo_indep_fal <- combo_indep_fal %>% mutate(geoBins = cut(V6/1000, breaks = c(0,50,100,150,500), labels = c("<50 km", "50-100 km", "100-150 km", ">150 km")))
combo_indep_viv <- combo_indep_viv %>% mutate(geoBins = cut(V6/1000, breaks = c(0,50,100,150,500), labels = c("<50 km", "50-100 km", "100-150 km", ">150 km")))

combo_indep_fal_PLOT <-   ggplot(data=combo_indep_fal, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max)) + geom_smooth(method="lm", color="black", lty="dashed",fullrange=T) + 
    geom_point(aes(fill=geoBins), shape=21, size=9) +
    #geom_point(shape=21, size=9, col="black") +
  labs(x="\nCase flow (as % of total)", y=paste("Frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), size="# pairwise comparisons (gen)", color="Species") +
  theme(axis.title = element_text(size=20)) +
  scale_x_continuous(breaks=c(0,5,10,15), labels=c("0%","5%","10%","15%")) + scale_fill_manual(values=rev(c("#fff9fb","#f79ebf","#D41159","#790a33"))) +
  scale_y_continuous(trans='sqrt', breaks=c(0,25,50,75,100), labels=c("0%","25%","50%","75%","100%")) + scale_size_continuous(range=c(3,15)) + theme_minimal() + theme(axis.text = element_text(size=22)) + theme(axis.title = element_text(size=20)) +
  geom_text_repel(data=combo_indep_fal, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, label = num), alpha=0, segment.color = NA)+ coord_cartesian(ylim = c(-0.8,100), xlim=c(-0.2,18))

combo_indep_viv_PLOT <- ggplot(data=combo_indep_viv, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max)) + geom_smooth(method="lm", color="black", lty="dashed",fullrange=T) +
  geom_point(aes(fill=geoBins), shape=21, size=9) +
  #geom_point(shape=21, size=9, col="black") +
  labs(x="\nCase flow (as % of total)", y=paste("Frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), size="# pairwise comparisons (gen)", color="Species") +
  theme(axis.title = element_text(size=20)) +
  scale_x_continuous(breaks=c(0,5,10,15), labels=c("0%","5%","10%","15%")) + scale_fill_manual(values=rev(c("#deeeff","#7cb9ff","#1A85FF","#003168"))) +
  scale_y_continuous(trans='sqrt', breaks=c(0,25,50,75,100), labels=c("0%","25%","50%","75%","100%")) + scale_size_continuous(range=c(3,15)) + theme_minimal() + theme(axis.text = element_text(size=22)) + theme(axis.title = element_text(size=20)) +
  geom_text_repel(data=combo_indep_viv, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, label = num), alpha=0, segment.color = NA)+ coord_cartesian(ylim = c(-0.8,100), xlim=c(-0.2,18))

temp_gf <- rbind(combo_indep_fal,combo_indep_viv)
temp_gf$geospec <- paste(temp_gf$geoBins, temp_gf$spec)

gen_vs_case_flow_plot <-   temp_gf %>%
  
  ggplot(., aes(x=100*sum/totspecsum, y=100*gen_as_frac_max)) + geom_smooth(method="lm", color="black", lty="dashed",fullrange=T) +
  geom_point(aes(fill=geospec), shape=21, size=8) +
  labs(x="\nCase flow (as % of total)", y=paste("Frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), size="# pairwise comparisons (gen)", color="Species") +
  theme(axis.title = element_text(size=20)) + 
  scale_fill_manual(breaks=c(">150 km falciparum","100-150 km falciparum", "50-100 km falciparum","<50 km falciparum",
                             ">150 km vivax","100-150 km vivax", "50-100 km vivax","<50 km vivax"), 
                              values=c("#fff9fb","#f79ebf","#D41159","#790a33","#deeeff","#7cb9ff","#1A85FF","#003168")) +
  scale_x_continuous(breaks=c(0,5,10,15), labels=c("0%","5%","10%","15%")) +
  scale_y_continuous(trans='sqrt', breaks=c(0,25,50,75,100), labels=c("0%","25%","50%","75%","100%")) + scale_size_continuous(range=c(3,15)) + theme_minimal() + theme(axis.text = element_text(size=22)) + theme(axis.title = element_text(size=20)) +
  geom_text_repel(max.overlaps=20, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, label = num), 
                  alpha=0, 
                  size=5, segment.color = NA) + coord_cartesian(ylim = c(-0.8,100), xlim=c(-0.2,18)) + 
  facet_grid(.~spec) + theme(panel.spacing = unit(1, "cm", data = NULL)) + theme(legend.position = "none")

cor.test(temp_gf$gen_as_frac_max[temp_gf$spec=="falciparum"], temp_gf$sum[temp_gf$spec=="falciparum"]/temp_gf$totspecsum[temp_gf$spec=="falciparum"])
cor.test(temp_gf$gen_as_frac_max[temp_gf$spec=="vivax"], temp_gf$sum[temp_gf$spec=="vivax"]/temp_gf$totspecsum[temp_gf$spec=="vivax"])

ggsave(paste("gen_vs_case_flow_", c("abs50","p99")[mode], ".pdf", sep=""),
       gen_vs_case_flow_plot,
       device = "pdf",
       units = "in",
       width=16,
       height=8,
       dpi = 600)

y=combo_indep$sum[combo_indep$spec=="falciparum"]/combo_indep$totspecsum[combo_indep$spec=="falciparum"]
x=combo_indep$gen_as_frac_max[combo_indep$spec=="falciparum"]
summary(lm(y~x))
y=combo_indep$sum[combo_indep$spec=="vivax"]/combo_indep$totspecsum[combo_indep$spec=="vivax"]
x=combo_indep$gen_as_frac_max[combo_indep$spec=="vivax"]
summary(lm(y~x))

### using only shared epi-positive site pairs for each species

ggplot(aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, color=spec, size=V4), data=combo_indep_shared_cf_detection) + geom_smooth(method="lm") + geom_point(alpha=0.7) +
  labs(x="\nCase flow (as % of total)", y=paste("Frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), size="# pairwise comparisons (gen)", color="Species") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + theme(axis.title = element_text(size=20)) +
  coord_cartesian(ylim = c(0,100)) + scale_y_continuous(trans='sqrt') + scale_size_continuous(range=c(3,12)) + facet_grid(.~spec) + theme_minimal() + theme(axis.text = element_text(size=12)) + theme(axis.title = element_text(size=20)) +
  geom_text_repel(data=combo_indep_shared_cf_detection, aes(x=100*sum/totspecsum, y=100*gen_as_frac_max, label=epi_bi), size=3, segment.color = NA)

y=combo_indep_shared_cf_detection$sum[combo_indep_shared_cf_detection$spec=="falciparum"]/combo_indep_shared_cf_detection$totspecsum[combo_indep_shared_cf_detection$spec=="falciparum"]
x=combo_indep_shared_cf_detection$gen_as_frac_max[combo_indep_shared_cf_detection$spec=="falciparum"]
summary(lm(y~x))

y=combo_indep_shared_cf_detection$sum[combo_indep_shared_cf_detection$spec=="vivax"]/combo_indep_shared_cf_detection$totspecsum[combo_indep_shared_cf_detection$spec=="vivax"]
x=combo_indep_shared_cf_detection$gen_as_frac_max[combo_indep_shared_cf_detection$spec=="vivax"]
summary(lm(y~x))

par(mar=c(8,8,0,0))

epi_vs_gen_relaxed_n <- Pv_bi %>% left_join(y = Pf_bi, by = c("Pv_epi_bi" = "Pf_epi_bi")) %>% left_join(y = Pv_gen, by = c("Pv_epi_bi" = "Pv_gen_link")) %>% left_join(y = Pf_gen, by = c("Pv_epi_bi" = "Pf_gen_link")) %>% 
  filter(is.na(Pv_V4)==F & is.na(Pf_V4)==F & is.na(Pv_sum)==F & is.na(Pf_sum)==F) %>% filter(Pv_V4>50 & Pf_V4>50) %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo", Pv_epi_bi))
epi_vs_gen_relaxed_n$Pf_gen_as_frac <- epi_vs_gen_relaxed_n$Pf_V5/max(epi_vs_gen_relaxed_n$Pf_V5)
epi_vs_gen_relaxed_n$Pv_gen_as_frac <- epi_vs_gen_relaxed_n$Pv_V5/max(epi_vs_gen_relaxed_n$Pv_V5)
plot(epi_vs_gen_relaxed_n$Pv_gen_as_frac ~ epi_vs_gen_relaxed_n$Pf_gen_as_frac, pch=19, col=rgb(0,0,1,alpha=0.5), cex=4, xlab="", ylab="")
text(epi_vs_gen_relaxed_n$Pv_gen_as_frac ~ epi_vs_gen_relaxed_n$Pf_gen_as_frac, label=epi_vs_gen_relaxed_n$Pv_epi_bi, pch=19, col="black", cex=0.6)
mtext(side=1, paste("\n\nP. falciparum: frequency of ", thresh, " IBD\n(as % max. value observed)", sep=""), line=4.25)
mtext(side=2, paste("P. vivax: frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), line=3.5)

ggplot(data=epi_vs_gen_relaxed_n, aes(x=Pf_gen_as_frac, y=Pv_gen_as_frac)) + geom_smooth(method="lm", col="black", size=0.3, linetype="dashed") + geom_point(size=6, alpha=0.7) + 
  geom_text_repel(data=epi_vs_gen_relaxed_n, aes(x=Pf_gen_as_frac, y=Pv_gen_as_frac, label=Pv_epi_bi), size=3, segment.color = NA) + 
  labs(x=paste("\nP. falciparum: frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep=""), y=paste("P. vivax: frequency of ", thresh, " IBD\n(as % max. value observed)\n", sep="")) +
  theme_minimal() + theme(axis.title.x = element_text(size=20, color="#D41159")) +
  theme(axis.title.y = element_text(size=20, color="#1A85FF")) +
  theme(axis.text = element_text(size=20))
  
y=epi_vs_gen_relaxed_n$Pf_gen_as_frac
x=epi_vs_gen_relaxed_n$Pv_gen_as_frac
cor.test(y,x)

epi_vs_gen_relaxed_n %>% pivot_wider("Pv_epi_bi")

# epi_vs_gen_Pv_LOC <- Pv_locs %>% mutate(District_sampled2 = paste(District_sampled,District_sampled, sep="_to_")) %>% left_join(y = Pv_gen, by = c("District_sampled2" = "Pv_gen_link")) %>% filter(Pv_V4>100 & n_cases>5) %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo", District_sampled2))
# epi_vs_gen_Pf_LOC <- Pf_locs %>% mutate(District_sampled2 = paste(District_sampled,District_sampled, sep="_to_")) %>% left_join(y = Pf_gen, by = c("District_sampled2" = "Pf_gen_link")) %>% filter(Pf_V4>100 & n_cases>5) %>% filter(!grepl("Greater_GT|West_of_GT|Lower_Essequibo", District_sampled2))
# colnames(epi_vs_gen_Pv_LOC) <- gsub("^Pv_", "", colnames(epi_vs_gen_Pv_LOC))
# colnames(epi_vs_gen_Pf_LOC) <- gsub("^Pf_", "", colnames(epi_vs_gen_Pf_LOC))
# epi_vs_gen_Pf_LOC$spec <- "falciparum"
# epi_vs_gen_Pv_LOC$spec <- "vivax"
# combo_indep_LOC <- rbind(epi_vs_gen_Pv_LOC,epi_vs_gen_Pf_LOC)
# combo_indep_LOC$gen_as_frac_max <- "NA"
# combo_indep_LOC$gen_as_frac_max[combo_indep_LOC$spec=="vivax"] <- combo_indep_LOC$V5[combo_indep_LOC$spec=="vivax"]/max(combo_indep_LOC$V5[combo_indep_LOC$spec=="vivax"])
# combo_indep_LOC$gen_as_frac_max[combo_indep_LOC$spec=="falciparum"] <- combo_indep_LOC$V5[combo_indep_LOC$spec=="falciparum"]/max(combo_indep_LOC$V5[combo_indep_LOC$spec=="falciparum"])
# combo_indep_LOC$gen_as_frac_max <- as.numeric(combo_indep_LOC$gen_as_frac_max)

ggplot(aes(x=100*weight, y=100*gen_as_frac_max, color=spec, size=V4), data=combo_indep_LOC) + geom_point(alpha=0.7) + 
  labs(x="\nLocal case retention (%)", y="Frequency of >0.5 IBD\n(as % max. value observed)\n", size="# pairwise comparisons (gen)", color="Species") +
  #labs(x="\nCase flow (as % of total)", y="Frequency of p99 IBD\n(as % max. value observed)\n", size="# pairwise comparisons (gen)", color="Species") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + geom_smooth(method="lm") + theme(axis.title = element_text(size=20)) +
  coord_cartesian(ylim = c(0,100)) + scale_y_continuous(trans='sqrt') + scale_size_continuous(range=c(3,12)) + theme_minimal() + theme(axis.text = element_text(size=12)) + theme(axis.title = element_text(size=20)) 

y=combo_indep_LOC$weight[combo_indep_LOC$spec=="falciparum"]
x=combo_indep_LOC$gen_as_frac_max[combo_indep_LOC$spec=="falciparum"]
summary(lm(y~x))

y=combo_indep_LOC$weight[combo_indep_LOC$spec=="vivax"]
x=combo_indep_LOC$gen_as_frac_max[combo_indep_LOC$spec=="vivax"]
summary(lm(y~x))

plot(epi_vs_gen$Pf_sum/sum(epi_vs_gen$Pf_sum), epi_vs_gen$Pf_V5/max(epi_vs_gen$Pf_V5), pch=19, col="#D41159", cex=3, ylim=c(0,1), xlim=c(0,max(max(epi_vs_gen$Pv_sum)/sum(epi_vs_gen$Pv_sum), max(epi_vs_gen$Pf_sum)/sum(epi_vs_gen$Pf_sum))))
par(new=T)
plot(epi_vs_gen$Pv_sum/sum(epi_vs_gen$Pv_sum), epi_vs_gen$Pv_V5/max(epi_vs_gen$Pv_V5), pch=19, col="#1A85FF", cex=3, ylim=c(0,1), xlim=c(0,max(max(epi_vs_gen$Pv_sum)/sum(epi_vs_gen$Pv_sum), max(epi_vs_gen$Pf_sum)/sum(epi_vs_gen$Pf_sum))))
text(epi_vs_gen$sum, epi_vs_gen$V5, label=epi_vs_gen$Pf_epi_bi, pch=19, col="#D41159", cex=0.5, ylim=c(0,1), xlim=c(0,max(max(epi_vs_gen$Pv_sum)/sum(epi_vs_gen$Pv_sum), max(epi_vs_gen$Pf_sum)/sum(epi_vs_gen$Pf_sum))))

temppf <- epi_vs_gen[,c("Pv_epi_bi", "Pf_V5", "Pf_sum", "Pf_V4")]
colnames(temppf) <- c("link","gen","epi", "n_gen")
temppf$spec <- "falciparum"
temppv <- epi_vs_gen[,c("Pv_epi_bi", "Pv_V5", "Pv_sum", "Pv_V4")]
colnames(temppv) <- c("link","gen","epi", "n_gen")
temppv$spec <- "vivax"
combo_gen_epi <- rbind(temppf, temppv)

ggplot(aes(x=100*epi/sum(epi), y=100*gen/max(gen), color=spec), data=combo_gen_epi) + geom_point(size=5) + 
  labs(x="\n case flow\nas % of total case flow", y="fraction comparisons highly related\nas % of max value observed\n") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + geom_smooth(method="lm") + theme(axis.title = element_text(size=20)) 

combo_gen_epi_unspec <- combo_gen_epi %>% group_by(link) %>% summarise(unspec_epi=sum(epi))  %>% left_join(y=combo_gen_epi, by = c("link" = "link"))

ggplot(aes(x=100*unspec_epi/sum(unspec_epi), y=100*gen/max(gen), color=spec), data=combo_gen_epi_unspec) + geom_point(size=5) + 
  labs(x="\n case flow\nas % of total case flow", y="fraction comparisons highly related\nas % of max value observed\n") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + geom_smooth(method="lm") + theme(axis.title = element_text(size=20))

combo_gen_epi_unspec <- combo_gen_epi %>% group_by(link) %>% summarise(unspec_epi=sum(epi))  %>% left_join(y=combo_gen_epi, by = c("link" = "link"))

combo_gen_epi_unspec$gen_as_frac_max <- "NA"
combo_gen_epi_unspec$gen_as_frac_max[combo_gen_epi_unspec$spec=="vivax"] <- combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="vivax"]/max(combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="vivax"])
combo_gen_epi_unspec$gen_as_frac_max[combo_gen_epi_unspec$spec=="falciparum"] <- combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="falciparum"]/max(combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="falciparum"])

combo_gen_epi_unspec$gen_as_frac_max <- as.numeric(combo_gen_epi_unspec$gen_as_frac_max)

ggplot(aes(x=100*unspec_epi/sum(unspec_epi), y=100*gen_as_frac_max, color=spec, size=n_gen), data=combo_gen_epi_unspec) + geom_smooth(method="lm") + geom_point() + 
  labs(x="\n case flow as % of total case flow", y="fraction pairwise comparisons highly related\nas % of max. fraction pairwise relatedness observed\n", size="# pairwise comparisons (gen)", color="Species") +
  scale_color_manual(values=c("#D41159","#1A85FF")) + facet_grid(.~spec) + theme(axis.text = element_text(size=12)) + theme(axis.title = element_text(size=20)) + 
  coord_cartesian(ylim = c(0,100)) + scale_y_continuous(trans='sqrt')

y=combo_gen_epi_unspec$unspec_epi[combo_gen_epi_unspec$spec=="falciparum"]/sum(combo_gen_epi_unspec$unspec_epi[combo_gen_epi_unspec$spec=="falciparum"])
x=combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="falciparum"]/max(combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="falciparum"])
summary(lm(y~x))

y=combo_gen_epi_unspec$unspec_epi[combo_gen_epi_unspec$spec=="vivax"]/sum(combo_gen_epi_unspec$unspec_epi[combo_gen_epi_unspec$spec=="vivax"])
x=combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="vivax"]/max(combo_gen_epi_unspec$gen[combo_gen_epi_unspec$spec=="vivax"])
summary(lm(y~x))
