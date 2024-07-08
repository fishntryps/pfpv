
###################################################################################
### this script correlates polyclonality to case counts and estimated incidence ###
###################################################################################

### run time is seconds to minutes

### clear workspace and load libraries in R 4.2.2

rm(list=ls())

library(tidyverse) 
library(RColorBrewer)
library(patchwork)
library(swfscMisc)
library(raster)
library(dplyr)
library(geosphere)
library(raster)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel)

# load dataframe of individual cases and locality catalogue

main_path <-"C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana"
setwd(main_path)

df_raw <- read.csv(paste(main_path,"data_clean.csv",sep="/"))

geo_groups <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/COORDS_grouped_v4.tab", header=F)

### some clean-up required

geo_groups <- geo_groups[!geo_groups$V7 %in% "Venezuela",]
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

df$geosample <- geo_groups$V7[match(df$Locality_sampled, geo_groups$V1)]
df$geoinfect <- geo_groups$V7[match(df$Locality_infected, geo_groups$V1)]
df$geosample[is.na(df$geosample)] <- "unknown"
df$geoinfect[is.na(df$geoinfect)] <- "unknown"
df$geoinfect_to_sample <- paste(df$geoinfect,df$geosample, sep="_to_")


### here we find the projected population counts (LandScan) for zones within within Guyana

mysites <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/COORDS_grouped_v4.tab", header=F)

### plot LandScan raster

par(mfrow=c(1,1))

 r <- raster("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/landscan-global-2021.tif")
 e <- as(extent(-63.5, -55.5, 1.2, 9.2), 'SpatialPolygons')
 crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"
 r <- crop(r, e)

par(mar=c(0.2,0.2,0.2,0.2))
plot(r)
map <- rworldmap::getMap(resolution = "high")
map <- broom::tidy(map)
polygon(map[map$id=="Guyana",], lwd=0.2, lty="dashed")

piesgroups <- c("Central_Coast","Chenapau","Cristinas_Border","East_of_GT","Greater_Annai","Greater_GT","Greater_Lethem","Head_Mazaruni","Head_Waini","Kaituma_and_Barima","Lower_Cuyuni","Lower_Essequibo","Lower_Mazaruni","Lower_Potaro","Mid_Berbice","Mid_Essequibo","Mid_Mazaruni_and_Issano_Rd","North_Delta","Pakaraima_South","Upper_Cuyuni","Upper_Demerara","Upper_Mazaruni","Upper_Rupununi","Waini","West_of_GT", "Head_Essequibo", "Mid_New", "Lower_New", "Greater_Apoteri")
piescols <- c("#ccd87c","#91e780","#df81c3","#4497d6","#c9779d","#e91046","#858fea","#b956cc","#d72124","#d57acf","#d1c025","#3c22de","#7f9edc","#9fdd32","#4dd795","#72f0d1","#71d577","#cd7c65","#d9bb78","#6212ec","#b565e9","#4dbddb","#0ec843","#df7c2b","#9fe86e", "#00ff00", "#00ffff", "#ff99ff", "#ffff66")
COLS <- as.data.frame(cbind(piesgroups,piescols))

### this loop plots circles around zones of interest and computes population size within each circle

forpolygons <- geo_groups

 forpolygons_quadverts <- forpolygons[!forpolygons$V7 %in% c("Head_Essequibo","Lower_New","Mid_New"),] %>% filter(V7 %in% c("Kaituma_and_Barima","Lower_Cuyuni","Lower_Essequibo","Lower_Mazaruni","Lower_Potaro","Mid_Essequibo","Mid_Mazaruni_and_Issano_Rd","Upper_Cuyuni","Upper_Mazaruni","Greater_GT"))
 totpop <- list()
 totcells <- list()
 mykm2 <- c()
 sites <- levels(factor(forpolygons_quadverts$V7))
 for (i in sites){
   print(i)
   #zzz <- forpolygons_quadverts[forpolygons_quadverts$V7==i,][,8:9]
   # polygon(zzz$V8, zzz$V9)
   zzz <- as.data.frame(circle.polygon(mean(forpolygons_quadverts$V8[forpolygons_quadverts$V7==i]), mean(forpolygons_quadverts$V9[forpolygons_quadverts$V7==i]), 24.47)) # 24.47 arcmin = 45 km
   colnames(zzz) <- c("V8","V9")
   #polygon(circle.polygon(mean(forpolygons_quadverts$V8[forpolygons_quadverts$V7==i]), mean(forpolygons_quadverts$V9[forpolygons_quadverts$V7==i]), 15.33))
   poly1 <- sp::Polygon(makePoly(zzz[chull(zzz$V8, zzz$V9),]))
   firstPoly <- sp::Polygons(list(poly1), ID = "A")
   firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
   #r_poly <- raster::mask(r, rgeos::gBuffer(firstSpatialPoly, width=0.05))
   r_poly <- raster::mask(r, firstSpatialPoly)
   totpop <- append(totpop, sum(r_poly@data@values, na.rm=T))
   totcells <- append(totcells, length(r_poly@data@values[!is.na(r_poly@data@values)]))
   plot(rgeos::gBuffer(firstSpatialPoly, width=0.03),add=T, col=rgb(1,0,0, alpha=0.4), lty=0)
   #polygon(x=zzz[chull(zzz$V8, zzz$V9),1], y=zzz[chull(zzz$V8, zzz$V9),2],
     #      col=rgb(as.vector(col2rgb(COLS$piescols[COLS$piesgroups==i])/255)[1], as.vector(col2rgb(COLS$piescols[COLS$piesgroups==i])/255)[2], as.vector(col2rgb(COLS$piescols[COLS$piesgroups==i])/255)[3], alpha=0.6), lty=1)
   
   mykm2 <- c(mykm2, length(r_poly@data@values[!is.na(r_poly@data@values)]) * 0.86)
   rm(poly1)
   rm(firstPoly)
   rm(firstSpatialPoly)
   rm(r_poly)
   rm(zzz)
 }
 
totpop <- as.data.frame(as.matrix(totpop))
colnames(totpop) <- "totpop"
totpop$geo <- sites
totpop$totpop <- as.vector(unlist(totpop$totpop))

totcells <- as.data.frame(as.matrix(totcells))
colnames(totcells) <- "totcells"
totcells$geo <- sites
totcells$totcells <- as.vector(unlist(totcells$totcells))


### load in polyclonality information for each species and combine

Pf <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/Pf_CLEAN_Venez_numdate_hetfrac_COIL_EXT.txt", header=T)
Pv <- read.table("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana/Pv_CLEAN_Venez_numdate_hetfrac_COIL_EXT.txt", header=T)

Pf$spec <- "Pf"
Pv$spec <- "Pv"

x <- rbind(Pf,Pv)

x$age_group <- cut(x$AgeInYearMORE , 
                   breaks=c(0,15,55,120), 
                   include.lowest=TRUE, 
                   right=FALSE, 
                   labels=c("under15","15to55","over55"))

x$age_group <- factor(x$age_group, levels=c("under15","15to55","over55"))

x$coil_p025 <- as.factor(x$coil_p025)

x$WhereInfecTED <- gsub("_"," ", x$WhereInfecTED)
x$WhereInfecTED <- gsub("\\.","", x$WhereInfecTED)
x$WhereInfecTED <- tolower(x$WhereInfecTED)

x$Species <- gsub("Pf","falciparum",x$spec)
x$Species <- gsub("Pv","vivax",x$Species)
x$geo <- geo_groups$V7[match(x$WhereInfecTED, geo_groups$V1)]
x$geo[is.na(x$geo)] <- "unknown"

### we will only address zones for which genomic sample size is at least 20

x <- x[x$geospec %in%  names(table(x$geospec)[table(x$geospec)>=20]),]

x_agg1 <-  x %>% group_by(Species, coil_p025, geo) %>% filter(geo!="unknown") %>% summarise(n_cases=n()) 
epicases_geo <- df_agg1 %>% group_by(Species, geo) %>% filter(geo!="") %>% summarise(n_cases=sum(n_cases))

### here we compute polyclonality

xagg2 <- x_agg1 %>% pivot_wider(names_from = coil_p025, values_from = n_cases, values_fill = 0) 
#need <- c(as.character("4"),as.character("5"))
need <- c(as.character("3"), as.character("4"),as.character("5"))
xagg2[need[!(need %in% colnames(xagg2))]] = 0
xagg2$poly_any_frac <- 1-((xagg2$`1`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`))
xagg2$poly_2_frac <- (xagg2$`2`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)
xagg2$poly_3_frac <- (xagg2$`3`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)  
xagg2$poly_4_frac <- (xagg2$`4`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)  
xagg2$poly_5_frac <- (xagg2$`5`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)  
xagg2$gen_cases <- as.numeric(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)
gen_epi <- left_join(xagg2, epicases_geo, by=c('Species', 'geo'))

gen_epi$meancoi <- ((1 * gen_epi$`1`) + (2 * gen_epi$`2`) + (3 * gen_epi$`3`) + (4 * gen_epi$`4`) + (5 * gen_epi$`5`)) / (gen_epi$`1` + gen_epi$`2` + gen_epi$`3` + gen_epi$`4` + gen_epi$`5`)
sdlist <- list()
for (i in 1:length(gen_epi$meancoi)){
  sdlist <- append(sdlist, sd(c(rep(1, gen_epi$`1`[i]),rep(2, gen_epi$`2`[i]),rep(3, gen_epi$`3`[i]),rep(4, gen_epi$`4`[i]),rep(5, gen_epi$`5`[i]))))
}

p025list <- list()
for (i in 1:length(gen_epi$meancoi)){
  p025list <- append(p025list, quantile(c(rep(1, gen_epi$`1`[i]),rep(2, gen_epi$`2`[i]),rep(3, gen_epi$`3`[i]),rep(4, gen_epi$`4`[i]),rep(5, gen_epi$`5`[i])), probs = 0.025))
}


p975list <- list()
for (i in 1:length(gen_epi$meancoi)){
  p975list <- append(p975list, quantile(c(rep(1, gen_epi$`1`[i]),rep(2, gen_epi$`2`[i]),rep(3, gen_epi$`3`[i]),rep(4, gen_epi$`4`[i]),rep(5, gen_epi$`5`[i])), probs = 0.975))
}

pminlist <- list()
for (i in 1:length(gen_epi$meancoi)){
  pminlist <- append(pminlist, min(c(rep(1, gen_epi$`1`[i]),rep(2, gen_epi$`2`[i]),rep(3, gen_epi$`3`[i]),rep(4, gen_epi$`4`[i]),rep(5, gen_epi$`5`[i]))))
}

pmaxlist <- list()
for (i in 1:length(gen_epi$meancoi)){
  pmaxlist <- append(pmaxlist, max(c(rep(1, gen_epi$`1`[i]),rep(2, gen_epi$`2`[i]),rep(3, gen_epi$`3`[i]),rep(4, gen_epi$`4`[i]),rep(5, gen_epi$`5`[i]))))
}


gen_epi$sdcoi <- as.vector(unlist(sdlist))
gen_epi$p025 <- as.vector(unlist(p025list))
gen_epi$p975 <- as.vector(unlist(p975list))
gen_epi$min <- as.vector(unlist(pminlist))
gen_epi$max <- as.vector(unlist(pmaxlist))

### here we bootstrap 100x

new.matrix.i <- vector("list", 100)
for (i in 1:100){
new.matrix.i[[i]] <- x %>% group_by(geospec) %>% slice_sample(n=100, replace=T)}

next.matrix.i <- vector("list", 100)
for (i in 1:100){
s <- new.matrix.i[[i]]

x_agg1 <-  s %>% group_by(Species, coil_p025, geo) %>% filter(geo!="unknown") %>% summarise(n_cases=n()) 
epicases_geo <- df_agg1 %>% group_by(Species, geo) %>% filter(geo!="") %>% summarise(n_cases=sum(n_cases))


xagg2 <- x_agg1 %>% pivot_wider(names_from = coil_p025, values_from = n_cases, values_fill = 0) 
#need <- c(as.character("4"),as.character("5"))
need <- c(as.character("3"), as.character("4"),as.character("5"))
xagg2[need[!(need %in% colnames(xagg2))]] = 0
xagg2$poly_any_frac <- 1-((xagg2$`1`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`))
xagg2$poly_2_frac <- (xagg2$`2`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)
xagg2$poly_3_frac <- (xagg2$`3`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)  
xagg2$poly_4_frac <- (xagg2$`4`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)  
xagg2$poly_5_frac <- (xagg2$`5`)/(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)  
xagg2$gen_cases <- as.numeric(xagg2$`1`+ xagg2$`2`+ xagg2$`3` + xagg2$`4` + xagg2$`5`)
subsample_gen_epi <- left_join(xagg2, epicases_geo, by=c('Species', 'geo'))

next.matrix.i[[i]] <- subsample_gen_epi}

subsamples <- as.data.frame(next.matrix.i[1])[,c(1,2,8,14)]
for (i in 2:100){
subsamples <- cbind(subsamples, as.data.frame(next.matrix.i[[i]])[8])}
subsamples$mean_poly_any_frac <- apply(subsamples[, c(3,5:103)],1,mean)
subsamples$sd <- apply(subsamples[, c(3,5:103)],1,sd)
subsamples$lower <- apply(subsamples[, c(3,5:103)],1, function(x) quantile(x, probs=0.025))
subsamples$upper <- apply(subsamples[, c(3,5:103)],1, function(x) quantile(x, probs=0.975))

subsamples_with_meancoi <- as.data.frame(next.matrix.i[1])[,c(1,2,3,4,5,6,7)]
subsamples_with_meancoi$meancoi <- ((1 * subsamples_with_meancoi$X1) + (2 * subsamples_with_meancoi$X2) + (3 * subsamples_with_meancoi$X3) + (4 * subsamples_with_meancoi$X4) + (5 * subsamples_with_meancoi$X5)) / (subsamples_with_meancoi$X1 + subsamples_with_meancoi$X2 + subsamples_with_meancoi$X3 + subsamples_with_meancoi$X4 + subsamples_with_meancoi$X5)
subsamples_with_meancoi <- subsamples_with_meancoi[,c(1,2,8)]
for (i in 2:100){
  temp <- as.data.frame(next.matrix.i[i])[,c(1,2,3,4,5,6,7)]
  temp$meancoi <- ((1 * temp$X1) + (2 * temp$X2) + (3 * temp$X3) + (4 * temp$X4) + (5 * temp$X5)) / (temp$X1 + temp$X2 + temp$X3 + temp$X4 + temp$X5)
  subsamples_with_meancoi <- cbind(subsamples_with_meancoi, temp[,8])}
subsamples_with_meancoi$mean_meancoi <- apply(subsamples_with_meancoi[, c(3,102)],1,mean)
subsamples_with_meancoi$sd_meancoi <- apply(subsamples_with_meancoi[, c(3:102)],1,sd)
subsamples_with_meancoi$lower_meancoi <- apply(subsamples_with_meancoi[, c(3:102)],1, function(x) quantile(x, probs=0.025))
subsamples_with_meancoi$upper_meancoi <- apply(subsamples_with_meancoi[, c(3:102)],1, function(x) quantile(x, probs=0.975))

subsamples <- cbind(subsamples, subsamples_with_meancoi[,c(103:106)])

piesgroups <- c("Central_Coast","Chenapau","Cristinas_Border","East_of_GT","Greater_Annai","Greater_GT","Greater_Lethem","Head_Mazaruni","Head_Waini","Kaituma_and_Barima","Lower_Cuyuni","Lower_Essequibo","Lower_Mazaruni","Lower_Potaro","Mid_Berbice","Mid_Essequibo","Mid_Mazaruni_and_Issano_Rd","North_Delta","Pakaraima_South","Upper_Cuyuni","Upper_Demerara","Upper_Mazaruni","Upper_Rupununi","Waini","West_of_GT")
piescols <- c("#ccd87c","#91e780","#df81c3","#4497d6","#c9779d","#e91046","#858fea","#b956cc","#d72124","#d57acf","#d1c025","#3c22de","#7f9edc","#9fdd32","#4dd795","#72f0d1","#71d577","#cd7c65","#d9bb78","#6212ec","#b565e9","#4dbddb","#0ec843","#df7c2b","#9fe86e")
COLS <- as.data.frame(cbind(piesgroups,piescols))
subsamples$geocol <- COLS$piescols[match(subsamples$geo, COLS$piesgroups)]


subsamplestemp_pf <- subsamples[,c(1,2,4,104:112)] %>% filter(Species=="falciparum")
subsamplestemp_pv <- subsamples[,c(1,2,4,104:112)] %>% filter(Species=="vivax")

subsamplestemp_pf$init_poly_any_frac <- gen_epi$poly_any_frac[gen_epi$Species=="falciparum"][match(subsamplestemp_pf$geo, gen_epi$geo[gen_epi$Species=="falciparum"])]
subsamplestemp_pf$init_meancoi <- gen_epi$meancoi[gen_epi$Species=="falciparum"][match(subsamplestemp_pf$geo, gen_epi$geo[gen_epi$Species=="falciparum"])]
subsamplestemp_pf$gen_cases <- gen_epi$gen_cases[gen_epi$Species=="falciparum"][match(subsamplestemp_pf$geo, gen_epi$geo[gen_epi$Species=="falciparum"])]



subsamplestemp_pv$init_poly_any_frac <- gen_epi$poly_any_frac[gen_epi$Species=="vivax"][match(subsamplestemp_pv$geo, gen_epi$geo[gen_epi$Species=="vivax"])]
subsamplestemp_pv$init_meancoi <- gen_epi$meancoi[gen_epi$Species=="vivax"][match(subsamplestemp_pv$geo, gen_epi$geo[gen_epi$Species=="vivax"])]
subsamplestemp_pv$gen_cases <- gen_epi$gen_cases[gen_epi$Species=="vivax"][match(subsamplestemp_pv$geo, gen_epi$geo[gen_epi$Species=="vivax"])]


myvec <- vector()
for (i in levels(factor((subsamplestemp_pv$geo)))){
  rgb <- subsamplestemp_pv$geocol[subsamplestemp_pv$geo==i]
  myvec <- c(myvec, rgb)
}
subsamplestemp_pv$geocol <- as.factor(subsamplestemp_pv$geocol)
subsamplestemp_pv$n_cases_div_totpop <- subsamplestemp_pv$n_cases/totpop$totpop[match(subsamplestemp_pv$geo, totpop$geo)]
subsamplestemp_pv$n_cases_div_totpop_and_cells <- subsamplestemp_pv$n_cases/(totpop$totpop[match(subsamplestemp_pv$geo, totpop$geo)] * totcells$totcells[match(subsamplestemp_pv$geo, totcells$geo)])


### now we can plot polyclonality against estimated incidence (dividing by 12 months)




ggplot(both,aes(x=(1/12)*n_cases_div_totpop, y=init_poly_any_frac, shape=Species)) + 
  #geom_text_repel(aes(label=geo), adj=0, cex=3) +
  geom_smooth(method="lm", lwd=0.5, lty="dashed", col=c(rep(col.species[1],80), rep(col.species[2],80)), se=F) +
  geom_errorbar(aes(ymin=init_poly_any_frac - sd, ymax=init_poly_any_frac + sd, color=Species), width=(1/12)*max(subsamplestemp_pv$n_cases_div_totpop)/80) + ylim(0,0.42) + 
  xlim(0,max((1/12)*subsamplestemp_pv$n_cases_div_totpop)) +
  scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10), labels=c(breaks = c("0.00", "0.02", "0.04", "0.06", "0.08", "0.10")), limits=c(0,0.1)) +
  #scale_x_continuous(breaks = c(0.00000, 0.00005, 0.00010, 0.00015, 0.00020), labels=c(breaks = c("0.00000", "0.00005", "0.00010", "0.00015", "0.00020")), limits=c(0,0.000203)) +
  #scale_x_continuous(breaks = pretty(subsamplestemp_pv$n_cases_div_totpop, n = 6)) +
  #scale_x_continuous(breaks = pretty(subsamplestemp_pv$n_cases_div_totpop, n = 5)) +
  scale_y_continuous(breaks = c(0.0, 0.10, 0.20, 0.30, 0.40), labels=c("0.00", "0.10", "0.20", "0.30", "0.40"), limits=c(0,0.42)) +
  #xlim(0, 0.025) +
  geom_point(aes(size=gen_cases, color=Species)) + scale_size(range=c(2,14)) +
  labs(color="Transmission zone", size="Genetic sample size") +
  ylab("Poly-clonality rate\n") + xlab("\n# cases divided by LandScan population projection (2019)") + scale_color_manual(values=myvec) + #xlim(-0.0001,0.008) +
  discrete_scale('shape', 'shape_d', function (n) c(16,16)[seq_len(n)]) + scale_color_manual(values=c("#D41159","#1A85FF")) +
  #geom_point(aes(size=gen_cases_div_totpop, x=n_cases_div_totpop, y=init_poly_any_frac), color="black", shape=10, data=both[both$Species=="vivax",]) + scale_size(range=c(2,13)) +
  #geom_point(aes(size=gen_cases_div_totpop, x=n_cases_div_totpop, y=init_poly_any_frac), color="black", shape=1, data=both[both$Species=="falciparum",]) + scale_size(range=c(2,13)) +
  theme(axis.title.x = element_text(size=15)) + 
  theme(axis.title.y = element_text(size=15)) +
  theme(legend.title=element_text(size=12, face="bold")) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.y = element_text(size=14, color="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1,size=14, color="black")) + 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Kaituma_and_Barima"],
                   y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Kaituma_and_Barima"], label= "1", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Lower_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Cuyuni"], label= "2", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Lower_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Essequibo"], label= "3", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Lower_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Mazaruni"], label= "4", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Lower_Potaro"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Potaro"], label= "5", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Mid_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Mid_Essequibo"], label= "6", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Mid_Mazaruni_and_Issano_Rd"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Mid_Mazaruni_and_Issano_Rd"], label= "7", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Upper_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Upper_Cuyuni"], label= "8", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="falciparum" & both$geo=="Upper_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Upper_Mazaruni"], label= "9", size=3.5, col="black") +

annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Kaituma_and_Barima"],
         y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Kaituma_and_Barima"], label= "1", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Lower_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Cuyuni"], label= "2", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Lower_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Essequibo"], label= "3", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Lower_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Mazaruni"], label= "4", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Lower_Potaro"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Potaro"], label= "5", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Mid_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Mid_Essequibo"], label= "6", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Mid_Mazaruni_and_Issano_Rd"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Mid_Mazaruni_and_Issano_Rd"], label= "7", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Upper_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Upper_Cuyuni"], label= "8", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Upper_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Upper_Mazaruni"], label= "9", size=3.5, col="black") +
  annotate("text", x=(1/12)*both$n_cases_div_totpop[both$Species=="vivax" & both$geo=="Greater_GT"],
         y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Greater_GT"], label= "10", size=3.5, col="black")

### correlations are not significant

cor.test(both$n_cases_div_totpop[both$Species %in% "falciparum"], both$init_poly_any_frac[both$Species %in% "falciparum"])
cor.test(both$n_cases_div_totpop[both$Species %in% "vivax"], both$init_poly_any_frac[both$Species %in% "vivax"])

ggplot(both,aes(x=n_cases, y=init_poly_any_frac, shape=Species)) + 
  #geom_text_repel(aes(label=geo), adj=0, cex=3) +
  geom_smooth(method="lm", lwd=0.5, lty="dashed", col=c(rep(col.species[1],80), rep(col.species[2],80)), se=F) +
  geom_errorbar(aes(ymin=init_poly_any_frac - sd, ymax=init_poly_any_frac + sd, color=Species), width=max(subsamplestemp_pv$n_cases)/80) + ylim(0,0.42) + 
  xlim(0,max(subsamplestemp_pv$n_cases)) +
  #scale_x_continuous(breaks = pretty(subsamplestemp_pv$n_cases, n = 5), limits=c(0,0.0143)) +
  scale_x_continuous(breaks = pretty(subsamplestemp_pv$n_cases, n = 5), labels=c("0.00", "0.00", "0.00", "0.00", "0.00", "0.00"), limits=c(0,1.05*max(subsamplestemp_pv$n_cases))) +
  scale_y_continuous(breaks = c(0.0, 0.10, 0.20, 0.30, 0.40), labels=c("0.00", "0.10", "0.20", "0.30", "0.40"), limits=c(0,0.42)) +
  #xlim(0, 0.025) +
  geom_point(aes(size=gen_cases, color=Species)) + scale_size(range=c(2,14)) +
  labs(color="Transmission zone", size="Genetic sample size") +
  ylab("Poly-clonality rate\n") + xlab("\n# cases (2019)") + scale_color_manual(values=myvec) + #xlim(-0.0001,0.008) +
  discrete_scale('shape', 'shape_d', function (n) c(16,16)[seq_len(n)]) + scale_color_manual(values=c("#D41159","#1A85FF")) +
  #geom_point(aes(size=gen_cases, x=n_cases, y=init_poly_any_frac), color="black", shape=10, data=both[both$Species=="vivax",]) + scale_size(range=c(2,13)) +
  #geom_point(aes(size=gen_cases, x=n_cases, y=init_poly_any_frac), color="black", shape=1, data=both[both$Species=="falciparum",]) + scale_size(range=c(2,13)) +
  theme(axis.title.x = element_text(size=15)) + 
  theme(axis.title.y = element_text(size=15)) +
  theme(legend.title=element_text(size=12, face="bold")) +
  theme(legend.text=element_text(size=10)) +
  theme(axis.text.y = element_text(size=14, color="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1,size=14, color="black")) + 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Kaituma_and_Barima"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Kaituma_and_Barima"], label= "1", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Lower_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Cuyuni"], label= "2", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Lower_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Essequibo"], label= "3", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Lower_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Mazaruni"], label= "4", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Lower_Potaro"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Lower_Potaro"], label= "5", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Mid_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Mid_Essequibo"], label= "6", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Mid_Mazaruni_and_Issano_Rd"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Mid_Mazaruni_and_Issano_Rd"], label= "7", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Upper_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Upper_Cuyuni"], label= "8", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="falciparum" & both$geo=="Upper_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="falciparum" & both$geo=="Upper_Mazaruni"], label= "9", size=3.5, col="black") +
    
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Kaituma_and_Barima"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Kaituma_and_Barima"], label= "1", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Lower_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Cuyuni"], label= "2", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Lower_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Essequibo"], label= "3", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Lower_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Mazaruni"], label= "4", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Lower_Potaro"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Lower_Potaro"], label= "5", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Mid_Essequibo"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Mid_Essequibo"], label= "6", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Mid_Mazaruni_and_Issano_Rd"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Mid_Mazaruni_and_Issano_Rd"], label= "7", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Upper_Cuyuni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Upper_Cuyuni"], label= "8", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Upper_Mazaruni"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Upper_Mazaruni"], label= "9", size=3.5, col="black") +
  annotate("text", x=both$n_cases[both$Species=="vivax" & both$geo=="Greater_GT"],
           y=both$init_poly_any_frac[both$Species=="vivax" & both$geo=="Greater_GT"], label= "10", size=3.5, col="black")