
#####################################################################################################
### this script provides various summary plots and statistics from the 2019 malaria case database ###
#####################################################################################################

### run time is seconds to minutes

### we start fresh by removing stored objects, load packages and set working/output directories in R 4.2.2

rm(list=ls())

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(ggtext)
library(reemsplots2)
library(maps)
library(ggplot2)
library(RColorBrewer)
library(rgdal)
library(ggpattern)
library(magick)

main_path <-"C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana"

setwd(main_path)

### some colors to use later on

col.species <- c("#D41159","#1A85FF")
col.sex <- c( "grey24", "#0B775E")

### load dataframe of individual cases

df_raw <- read.csv(paste(main_path,"data_clean.csv",sep="/"))

### clean-up

geo_groups <- read.csv("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/COORDS_grouped_v3.csv", header=F)
geo_groups$V1 <- gsub("_"," ", geo_groups$V1)
geo_groups$V2 <- gsub("_"," ", geo_groups$V2)
geo_groups$V1 <- gsub("\\.","", geo_groups$V1)
geo_groups$V1 <- tolower(geo_groups$V1)
df_raw$Locality_infected <- gsub("\\.","", df_raw$Locality_infected)
df_raw$Locality_infected <- tolower(df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("b erbice","berbice", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("karaudanaw.*","karaudanawa", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("kato$","kato amer vil", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("mocomoco$","moco moco village", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("kariako mission$","kariako creek", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("shea$","shea amer vil$", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("yupukari amer vil$","yupukari", df_raw$Locality_infected)
#df_raw$Locality_infected <- gsub("paramakatoi$", "paramakatoi amer vil", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("demarara$", "demerara", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("kamarangken settlemen.*", "kamarangkent settlement", df_raw$Locality_infected)
df_raw$Locality_infected <- gsub("vergenogen", "vergenoegen", df_raw$Locality_infected)

df_raw$Locality_sampled <- gsub("\\.","", df_raw$Locality_sampled)
df_raw$Locality_sampled <- tolower(df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("b erbice","berbice", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("karaudanaw.*","karaudanawa", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("kato$","kato amer vil", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("mocomoco$","moco moco village", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("kariako mission$","kariako creek", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("shea$","shea amer vil$", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("yupukari amer vil$","yupukari", df_raw$Locality_sampled)
#df_raw$Locality_sampled <- gsub("paramakatoi$", "paramakatoi amer vil", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("demarara$", "demerara", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("kamarangken settlemen.*", "kamarangkent settlement", df_raw$Locality_sampled)
df_raw$Locality_sampled <- gsub("vergenogen", "vergenoegen", df_raw$Locality_sampled)

df_raw$geo <- geo_groups$V7[match(df_raw$Locality_infected, geo_groups$V1)]
df_raw$geo[df_raw$Locality_infected=="imported"] <- "imported"

# select only new cases, and passive detection

df <- df_raw %>% filter(Detection.method=="PASSIVE", Case == "NEW CASE")

# select only year of interest

df <- df %>% filter(Year_tested==2019)

### how many cases per species?

temp <- as.data.frame(table(df$InfParas))
sum(temp[grepl("FAL", temp$Var1),]$Freq)
sum(temp[grepl("VIV", temp$Var1),]$Freq)
sum(temp[grepl("MAL", temp$Var1),]$Freq)
sum(temp[grepl("", temp$Var1),]$Freq)

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
# create a broader age group category

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

df_agg <-  df %>% pivot_longer(cols = c("falciparum","vivax"),names_to = "Species") %>% 
  # aggregate 
  filter(!is.na(value)) %>% 
  group_by(Species, Region_infected, Sex, EthnicGrouP,Age.group)%>% summarise(n_cases=n()) 

monthlies <- df %>% pivot_longer(cols = c("falciparum","vivax"),names_to = "Species") %>% 
  # aggregate 
  filter(!is.na(value)) %>% 
  group_by(Species, Month_tested)%>% summarise(n_cases=n()) 

### cases by ethnicity

df_agg %>% group_by(Species,EthnicGrouP) %>% summarise(n_cases=sum(n_cases)) %>% 
  
  ggplot(.,aes(x=EthnicGrouP,y=n_cases,fill=Species))+ 
  geom_bar(position="stack",stat="identity",alpha=.9)+
  scale_fill_manual(values=col.species)+ ylab("Case count") + xlab("Ethnic group") +
  theme_minimal() + theme()

### cases by sex and age groups

age_sex_plot <-  df_agg %>% group_by(Species,Sex, Age.group) %>% summarise(n_cases=sum(n_cases)) %>% 
  
  ggplot(.,aes(x=Sex,y=n_cases,fill=Species))+ 
  geom_bar(position="stack",stat="identity",alpha=.9)+
  #scale_fill_manual(values=col.species)+ ylab("Case count") + xlab("Sex") +
  facet_grid(~Age.group) +
  #theme_minimal() + theme()
  
  scale_fill_manual(values=col.species)+ ylab("Case count (2019, GMOH)\n") + xlab("\nSex") +
  theme_minimal() + theme() + theme(axis.text.x = element_text(size=35)) +
  theme(axis.title.x = element_text(size=35)) + theme(axis.title.y = element_text(size=35)) +
  theme(strip.text = element_text(size=35)) +
  theme(axis.text.y = element_text(size=32)) + scale_x_discrete(labels=function(x) gsub("_", " ", x, fixed=TRUE))

### cases by ethnicity, sex, and age groups

df_agg %>% group_by(Species,EthnicGrouP,Age.group, Sex) %>% summarise(n_cases=sum(n_cases)) %>%
  ggplot(.,aes(x=Sex,y=n_cases,fill=factor(EthnicGrouP,levels=c("Afro Guyanese", "Amerindian", "East Indian","Mixed","Other") ))) + 
  geom_bar(position="stack",stat="identity",alpha=.9)+
  scale_fill_manual(name = "Ethnic group",values=c("#F98400","#333333","#666666","#CCCCCC","#EEEEEE"),
  )+ ylab("Case counts") +
  facet_grid(Species~Age.group) +
  theme_minimal()

spec_region_plot <- df_agg %>% group_by(Species,EthnicGrouP, Region_infected) %>% summarise(n_cases=sum(n_cases)) %>%
  filter(Region_infected %in% c("1","7","8")) %>%
  ggplot(.,aes(x=factor(EthnicGrouP,levels=c("Other","Afro Guyanese", "East Indian", "Amerindian","Mixed")),y=n_cases,fill=Species )) + 
  geom_bar(position="stack",stat="identity",alpha=.9)+
  facet_grid(Region_infected~.) +
  scale_fill_manual(values=col.species)+ ylab("Case count (2019, GMOH)\n") + xlab("\nEthnicity") +
  theme_minimal() + theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18, vjust=0.9)) +
  theme(axis.title.x = element_text(size=20)) + theme(axis.title.y = element_text(size=20)) +
  theme(strip.text.y = element_text(size=20, angle=0)) +
  theme(axis.text.y = element_text(size=17)) +  
  theme(panel.spacing = unit(0.65, "cm", data = NULL))

### proportion of Pf versus Pv, per age group, sex, and ethnicity

df_agg %>%  group_by(Species,EthnicGrouP,Age.group, Sex) %>% summarise(n_cases=sum(n_cases)) %>%
  # add columns of proportion infection in ethnicity x
  pivot_wider(names_from = EthnicGrouP, values_from = n_cases) %>% mutate(all_ethno = sum(`East Indian`,`Afro Guyanese`, Amerindian, Mixed, Other, na.rm=T)) %>% 
  mutate(`afro guyanese` = `Afro Guyanese`/all_ethno, amerindian = Amerindian/all_ethno, mixed = Mixed/all_ethno) %>% 
  pivot_longer(cols = c(`afro guyanese`,amerindian,mixed),names_to = "Propotion.ethnicity") %>% 
  
  ggplot(.,aes(x=Age.group,y=value,fill=Sex)) + 
  geom_bar(stat="identity",alpha=.9, position = "dodge")+
  scale_fill_manual(values = col.sex)+  ylab("Proportion infections found in ethnicity x") +
  facet_grid(Species ~ Propotion.ethnicity) +
  theme_minimal()

df$geosample <- geo_groups$V7[match(df$Locality_sampled, geo_groups$V1)]
df$geoinfect <- geo_groups$V7[match(df$Locality_infected, geo_groups$V1)]
df$geosample[is.na(df$geosample)] <- "unknown"
df$geoinfect[is.na(df$geoinfect)] <- "unknown"
df$geoinfect_to_sample <- paste(df$geoinfect,df$geosample, sep="_to_")

### load polyclonality information for each species and combine

Pf <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/Pf_CLEAN_Venez_numdate_hetfrac_COIL_EXT.txt", header=T)
Pv <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/Pv_CLEAN_Venez_numdate_hetfrac_COIL_EXT.txt", header=T)

Pf$spec <- "Pf"
Pv$spec <- "Pv"

x <- rbind(Pf,Pv)

# set up cut-off values 

qbreaks <- c(0,91,182,274,366,366+90,366+181,366+273,366+365)

# specify interval/bin labels

qtags <- c("2020-Q1","2020-Q2","2020-Q3","2020-Q4","2021-Q1","2021-Q2","2021-Q3","2021-Q4")

# bucketing values into bins

x$quarter <- cut(x$DaysFromStart2020, 
                 breaks=qbreaks, 
                 include.lowest=TRUE, 
                 right=FALSE, 
                 labels=qtags)

x$quarter <- factor(x$quarter, levels=qtags)

mbreaks <- c(0,31,60,91,121,152,182,213,244,274,305,335,366,366+31,366+59,366+90,366+120,366+151,366+181,366+212,366+243,366+273,366+304,366+334,366+365)
mtags <- c("Jan20","Feb20","Mar20","Apr20","May20","Jun20","Jul20","Aug20","Sep20","Oct20","Nov20","Dec20","Jan21","Feb21","Mar21","Apr21","May21","Jun21","Jul21","Aug21","Sep21","Oct21","Nov21","Dec21")

x$month <- cut(x$DaysFromStart2020, 
               breaks=mbreaks, 
               include.lowest=TRUE, 
               right=FALSE, 
               labels=mtags)

x$month <- factor(x$month, levels=mtags)

x$age_group <- cut(x$AgeInYearMORE , 
                   breaks=c(0,15,55,120), 
                   include.lowest=TRUE, 
                   right=FALSE, 
                   labels=c("under15","15to55","over55"))

x$age_group <- factor(x$age_group, levels=c("under15","15to55","over55"))

x$coil_p025 <- as.factor(x$coil_p025)

df_agg1 <-  df %>% pivot_longer(cols = c("falciparum","vivax"),names_to = "Species") %>% 
  # aggregate 
  filter(!is.na(value)) %>% 
  group_by(Species, Region_infected, Sex, EthnicGrouP,Age.group, District_infected, geo) %>% summarise(n_cases=n()) 

xx <- df_agg1 %>% group_by(geo) %>% filter(geo!="") %>% summarise(n_cases=sum(n_cases))
xxlevels <- xx$geo[order(xx$n_cases)]
xxlevels <- gsub("imported", "Imported", xxlevels)

pp <- df_agg1 %>% group_by(Species, geo) %>% filter(geo!="") %>% summarise(n_cases=sum(n_cases))

pp$geo <- gsub("imported", "Imported", pp$geo)
neworder <- c("Lower_Potaro",
              "Lower_Mazaruni",
              "Lower_Cuyuni",
              "Kaituma_and_Barima",
              "Cristinas_Border",
              "Mid_Mazaruni_and_Issano_Rd",
              "Upper_Cuyuni",
              "Waini",
              "Lower_Essequibo",
              "Upper_Rupununi",
              "Central_Coast",
              "Greater_GT",
              "North_Delta",
              "Upper_Mazaruni",
              "Mid_Essequibo",
              "Greater_Lethem",
              "West_of_GT",
              "Head_Waini",
              "Head_Mazaruni",
              "Greater_Annai",
              "Chenapau",
              "Pakaraima_South",
              "Upper_Demerara",
              "Mid_Berbice",
              "East_of_GT",
              "Greater_Apoteri",
              "Head_Essequibo")

pp$geo <- factor(pp$geo, levels=c(neworder, "Imported"))
pp <- pp[order(pp$geo),]

temp1 <- pp %>% filter(Species=="falciparum") 
temp12 <- pp %>% filter(Species=="vivax") %>% left_join(y=temp1, by=c("geo"="geo"))
temp12$n_cases.y[is.na(temp12$n_cases.y)] <- 0
temp12$Species.y <- "falciparum"
temp12$tot_n <- temp12$n_cases.x + temp12$n_cases.y  


temp12$geo <- factor(temp12$geo, levels=c(neworder, "Imported"))
cumul_vivax <- temp12[order(temp12$geo), 1:3] %>% mutate(cv=cumsum(n_cases.x))
colnames(cumul_vivax) <- c( "Species","geo", "n_cases", "cv")
cumul_order <- cumul_vivax$geo
cumul_falciparum <- temp12[,c(4,2,5)] %>% mutate(cv=cumsum(n_cases.y))
colnames(cumul_falciparum) <- c( "Species","geo", "n_cases", "cv")

CUMUL <- rbind(cumul_vivax, cumul_falciparum)

label_colors <- gsub("^Head_Waini$|^North_Delta$|^Waini$|^Kaituma_and_Barima$|^Central_Coast$","#ffb9fd", cumul_order)
label_colors <- gsub("^Head_Mazaruni$|^Upper_Mazaruni$|^Mid_Mazaruni_and_Issano_Rd$|^Lower_Mazaruni$|^Lower_Essequibo$|^Lower_Cuyuni$|^Upper_Cuyuni$|^Cristinas_Border$","#54e587", label_colors)
label_colors <- gsub("^Lower_Potaro$|^Chenapau$|^Pakaraima_South$|^Mid_Essequibo$","#8cebff", label_colors)
label_colors <- gsub("Head_Essequibo$|^Greater_Apoteri$|^East_of_GT$|^Mid_Berbice$|^Upper_Demerara$|^Greater_Annai$|^West_of_GT$|^Greater_Lethem$|^Greater_GT$|^Upper_Rupununi$", "#FFFF99", label_colors)

CUMUL$labelcolor <- CUMUL$geo

 CUMUL$labelcolor <- gsub("^Head_Waini$|^North_Delta$|^Waini$|^Kaituma_and_Barima$|^Central_Coast$","#ffb9fd", CUMUL$labelcolor)
 CUMUL$labelcolor <- gsub("^Head_Mazaruni$|^Upper_Mazaruni$|^Mid_Mazaruni_and_Issano_Rd$|^Lower_Mazaruni$|^Lower_Essequibo$|^Lower_Cuyuni$|^Upper_Cuyuni$|^Cristinas_Border$","#54e587",  CUMUL$labelcolor)
 CUMUL$labelcolor <- gsub("^Lower_Potaro$|^Chenapau$|^Pakaraima_South$|^Mid_Essequibo$","#8cebff",  CUMUL$labelcolor)
 CUMUL$labelcolor <- gsub("Head_Essequibo$|^Greater_Apoteri$|^East_of_GT$|^Mid_Berbice$|^Upper_Demerara$|^Greater_Annai$|^West_of_GT$|^Greater_Lethem$|^Greater_GT$|^Upper_Rupununi$", "#FFFF99",  CUMUL$labelcolor)


 CUMUL$geo  <- gsub("_", " ", CUMUL$geo)
 cumul_order <- gsub("_", " ", cumul_order)
 cumul_order <- gsub("Mid Mazaruni and ", "Mid Maz. and ", cumul_order)
 CUMUL$geo <- gsub("Mid Mazaruni and ", "Mid Maz. and ", CUMUL$geo)
 cumul_order <- gsub(" and ", " + ", cumul_order)
 CUMUL$geo <- gsub(" and ", " + ", CUMUL$geo)


cumulative_case_plot <- ggplot(CUMUL[CUMUL$geo!="Imported",], 
                           aes(x = factor(geo, levels=cumul_order[cumul_order!="Imported"]), y = as.numeric(cv), group = factor(Species), color = factor(Species))) + 
  geom_line(size=1.4, alpha=1) +
  geom_point(size=9.3, alpha=1, color="white") + theme_minimal() +
  geom_point(size=9.3, alpha=0.6) + theme_minimal() +
  
  scale_color_manual(values=c("#D41159","#1A85FF")) +
  labs(y="\n\nCumulative cases\n", x="Epidemiological zone") +
  #coord_flip(ylim=c(-10,10000), clip = "off") +
  theme(axis.text.y = element_text(size=30, color="black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_text(size=22)) + 
  theme(axis.title.y = element_text(size=22)) + #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks=c(-15000, seq(0,10000,1000)), label=c("",seq(0,10000,1000)), expand = expansion(mult = c(0.4, 0.1))) +
  theme(plot.margin = unit(c(1,1,1,7), "lines")) + 
  geom_richtext(label.padding=grid::unit(c(0,5,0,1), "pt"), aes(y = as.numeric(cv)-600), size=10.5, angle=90, label.size=0.0001, alpha=c(rep(0, 27), rep(0.4,27)),
                label=c(rep("",27),CUMUL$geo[CUMUL$geo!="Imported"][1:27]), color=NA, label.color=NA, fill=CUMUL$labelcolor[CUMUL$labelcolor!="Imported"], hjust=1, vjust=0.5) +
  geom_richtext(aes(y = as.numeric(cv)-600), size=10.5, angle=90, label.size=0.0001,
                label=c(rep("",27),CUMUL$geo[CUMUL$geo!="Imported"][1:27]), color="black", label.color=NA, fill=NA, hjust=1, vjust=0.5, nudge_x = 0.08) +
  geom_richtext(aes(y = as.numeric(cv)-40), size=7.5, label=c(1:27,1:27), color="black", label.color=NA, fill=NA) + theme(panel.grid.major.x = element_blank()) + 
  theme(plot.margin = unit(c(0.18,0,0,0), "cm"))
 
### map of epidemiological zones used in cumulative_case_plot above

mysites <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/COORDS_grouped_v4.tab", header=F)
piesgroups <- c("Central_Coast","Chenapau","Cristinas_Border","East_of_GT","Greater_Annai","Greater_GT","Greater_Lethem","Head_Mazaruni","Head_Waini","Kaituma_and_Barima","Lower_Cuyuni","Lower_Essequibo","Lower_Mazaruni","Lower_Potaro","Mid_Berbice","Mid_Essequibo","Mid_Mazaruni_and_Issano_Rd","North_Delta","Pakaraima_South","Upper_Cuyuni","Upper_Demerara","Upper_Mazaruni","Upper_Rupununi","Waini","West_of_GT", "Mid_New","Lower_New","Greater_Apoteri", "Venezuela")
piescols <- c("#ccd87c","#91e780","#df81c3","#4497d6","#c9779d","#e91046","#858fea","#b956cc","#d72124","#d57acf","#d1c025","#3c22de","#7f9edc","#9fdd32","#4dd795","#72f0d1","#71d577","#cd7c65","#d9bb78","#6212ec","#b565e9","#4dbddb","#0ec843","#df7c2b","#9fe86e", "blue2", "blue3", "blue4", "yellow")
COLS <- as.data.frame(cbind(piesgroups,piescols))
mysites$geocol <- COLS$piescols[match(mysites$V7, COLS$piesgroups)]

source("C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools.R")
source("C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/fx_gadm_sp_loadCountries.R")

guy0.spldf <- gadm_sp_loadCountries("GUY",level = 0, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds0/")
guy1.spldf <- gadm_sp_loadCountries("GUY",level = 1, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds1/")
guy2.spldf <- gadm_sp_loadCountries("GUY",level = 2, basefile = "C:/Users/philipp/OneDrive - Harvard University/Documents/github/AMPLseq_Colombia2/GADMTools/world.rds.files/rds2/")

guy1.spldf[["spdf"]][["NAME_1"]]
guy2.spldf[["spdf"]][["NAME_2"]]

COLS <- as.data.frame(cbind(piesgroups,piescols))

myplot <- 
  ggplot() +  geom_polygon(data=guy0.spldf$spdf, aes(x=long, y=lat, group = group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=guy1.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) + theme_minimal() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) +
  
  geom_polygon(data=guy1.spldf[["spdf"]][guy1.spldf[["spdf"]][["NAME_1"]] == "Barima-Waini",], aes(x=long, y=lat, group=group), fill = "#ffb9fd", alpha=0.5, color = "black", size=.3) +
  
  geom_polygon(data=guy1.spldf[["spdf"]][guy1.spldf[["spdf"]][["NAME_1"]] == "Cuyuni-Mazaruni",], aes(x=long, y=lat, group=group), fill = "#54e587", alpha=0.5, color = "black", size=.3) +
  
  geom_polygon(data=guy1.spldf[["spdf"]][guy1.spldf[["spdf"]][["NAME_1"]] == "Potaro-Siparuni",], aes(x=long, y=lat, group=group), fill = "#8cebff", alpha=0.5, color = "black", size=.3) +
  coord_quickmap() 

sites <- levels(factor(mysites$V7))

sites <- cumul_order[cumul_order!="Imported"]
sites1 <- c("Head_Waini","North_Delta","Waini","Kaituma_and_Barima","Central_Coast")
sites2 <- c("Head_Mazaruni","Upper_Mazaruni","Mid_Mazaruni_and_Issano_Rd","Lower_Mazaruni","Lower_Essequibo","Lower_Cuyuni","Upper_Cuyuni","Cristinas_Border")
sites3 <- c("Lower_Potaro","Chenapau","Pakaraima_South","Mid_Essequibo")
sites4 <- c("Head_Essequibo","Greater_Apoteri","East_of_GT","Mid_Berbice","Upper_Demerara","Greater_Annai","West_of_GT","Greater_Lethem","Greater_GT","Upper_Rupununi")

mysites <- rbind(mysites,
                 c(NA,NA,NA,NA,NA,NA,"Head_Essequibo",-58.56778+0.1, 1.763056),
                 c(NA,NA,NA,NA,NA,NA,"Head_Essequibo",-58.56778-0.1, 1.763056),
                 c(NA,NA,NA,NA,NA,NA,"Head_Essequibo",-58.56778, 1.763056+0.1),
                 c(NA,NA,NA,NA,NA,NA,"Head_Essequibo",-58.56778, 1.763056-0.1))
mysites$V8 <- as.numeric(mysites$V8)
mysites$V9 <- as.numeric(mysites$V9)

for (i in sites){
  zzz <- mysites[mysites$V7==i,][,8:9]

  myplot <- myplot + 
    geom_polygon_pattern(data=zzz[chull(zzz$V8, zzz$V9),],
                         aes(x = V8, y = V9), 
                         linewidth = 0.1,
                         pattern= "circle",
                         fill            = NA, 
                         pattern_spacing = 0.005, 
                         pattern_density = 0.005, 
                         pattern_fill    = 'brown') +
    coord_equal() 
  
  rm(zzz)
}

myorder <- c("Lower_Potaro","Lower_Mazaruni","Lower_Cuyuni",
             "Kaituma_and_Barima","Cristinas_Border",
             "Mid_Mazaruni_and_Issano_Rd","Upper_Cuyuni",
             "Waini","Lower_Essequibo","Upper_Rupununi","Central_Coast",
             "Greater_GT","North_Delta","Upper_Mazaruni","Mid_Essequibo",
             "Greater_Lethem","West_of_GT","Head_Waini","Head_Mazaruni",
             "Greater_Annai","Chenapau","Pakaraima_South",
             "Upper_Demerara","Mid_Berbice","East_of_GT")
myorder <- rev(cumul_order[cumul_order!="Imported"])

myorder <- as.data.frame(myorder)
colnames(myorder) <- "name"
myorder$id2 <- 1:27
geos <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/Pf_centroids_degrees_v3.txt")
colnames(geos) <- c("name","lon", "lat")
geos$id <- rownames(geos)

geos2 <- geos %>% left_join(y = myorder, by = c("name" = "name"))

geos2 <- geos2 %>% filter(name %in% myorder$name)

geos2$labelcol <- geos2$name
geos2$labelcol <- gsub("^Head_Waini$|^North_Delta$|^Waini$|^Kaituma_and_Barima$|^Central_Coast$","#ffb9fd", geos2$labelcol)
geos2$labelcol <- gsub("^Head_Mazaruni$|^Upper_Mazaruni$|^Mid_Mazaruni_and_Issano_Rd$|^Lower_Mazaruni$|^Lower_Essequibo$|^Lower_Cuyuni$|^Upper_Cuyuni$|^Cristinas_Border$","#54e587", geos2$labelcol)
geos2$labelcol <- gsub("^Lower_Potaro$|^Chenapau$|^Pakaraima_South$|^Mid_Essequibo$","#8cebff", geos2$labelcol)
geos2$labelcol <- gsub("Head_Essequibo$|^Greater_Apoteri$|^East_of_GT$|^Mid_Berbice$|^Upper_Demerara$|^Greater_Annai$|^West_of_GT$|^Greater_Lethem$|^Greater_GT$|^Upper_Rupununi$", "grey", geos2$labelcol)

numlabel_zoneplot <- myplot + geom_text(data = geos2, aes(x =lon, y = lat, label=id2), size=6, color=geos2$labelcol)
numlabel_zoneplot <- myplot + geom_text(data = geos2, aes(x =lon, y = lat, label=id2), size=9, color="black")
