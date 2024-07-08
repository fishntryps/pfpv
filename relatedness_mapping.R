
############################################################################################################
### this script maps frequency of >0.50 pairwise relatedness (IBD) among epidemiological zones of Guyana ###
############################################################################################################

### run time is seconds to minutes

### we start fresh by removing stored objects, load packages and set working directory in R 4.2.2

rm(list=ls())

library(assertthat)
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)
library(ggmap)
require(gridExtra)
library(grid)
library(scatterpie)
library(tidyverse) 
library(RColorBrewer)    
library(patchwork)
library(ggrepel)
library(MMWRweek)
library(zoo)
library(circlize)
library(webshot2)
library(dplyr)

setwd("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana")

#### we will address P. falciparum cases first, loading file containing frequency of >0.50 pairwise relatedness within and between epidemiological zones (computed via script: Pf_g50.sh)

gen <- read.table("Pf_perc_g50_amoung_groups_v4_with_distances_and_indcount_and_cccount.txt", header=F)

### analysis will exclude Venezuelan infections

gen <- gen %>% filter(V1!="Venezuela" & V2!="Venezuela")
gen$geosample <- gen$V1
gen$geoinfect <- gen$V2
gen$geosample[is.na(gen$geosample)] <- ""
gen$geoinfect[is.na(gen$geoinfect)] <- ""

df <- gen
df$District_infected <- pmin(df$geoinfect,df$geosample)
df$District_sampled <- pmax(df$geoinfect,df$geosample)

### classify comparisons as within or between-zone (carrying over 'local' and 'outsourced' terminology from case flow mapping)

df_district <- df %>% mutate(type = case_when(District_infected==District_sampled ~ "Local",
                                              (District_infected=="" | District_sampled=="") ~ "!Missing", 
                                              (District_infected!=District_sampled & District_infected!="" & District_sampled!="") ~ "Outsourced"))

df_district$District_infected_gr <- gsub("$","_gr", df_district$District_infected)
df_district$District_sampled_gr <- gsub("$","_gr", df_district$District_sampled)

df_district <- df_district[df_district$type!="!Missing",]

### we will be plotting percent pairs exhibiting >0.50 IBD and use total number of comparisons to size plotted elements

df_district$perc_clo <- df_district$V5*100
df_district$num_comps <- df_district$V4
df_district$unclonal_comps <- df_district$V4-df_district$V3

set.seed(123)  # set random generator state for the same output

### this contains centroid coordinates for each of the epidemiological zones

geos <- read.table("centroids_via_COORDS_grouped_v4.txt", header=T)
geos$id <- rownames(geos)

### we make a separate dataframe for 'outsourced' comparisons, also adding coordinates corresponding to each zone

dfOut <- df_district %>% filter(type=="Outsourced")
q75Out <- quantile(dfOut$perc_clo)[4]
sdOut <- sd(dfOut$perc_clo)
dfOut$x <- geos$lon[match(dfOut$District_infected, geos$name)]
dfOut$y <- geos$lat[match(dfOut$District_infected, geos$name)]
dfOut$xend <- geos$lon[match(dfOut$District_sampled, geos$name)]
dfOut$yend <- geos$lat[match(dfOut$District_sampled, geos$name)]
dfOut$from<- geos$id[match(dfOut$District_infected, geos$name)]
dfOut$to<- geos$id[match(dfOut$District_sampled, geos$name)]
dfOut$weight <- dfOut$perc_clo
dfOut$category <- dfOut$from
dfOut$maxcc <- pmax(dfOut$V9,dfOut$V10)

locsizes <- read.table("Pf_sampsizes_grouped_v4_myIBD_with_locallyclonalnums.txt", header=T)
dfOut$locsize1 <- locsizes$numsamps[match(dfOut$District_infected, locsizes$loc)]
dfOut$locsize2 <- locsizes$numsamps[match(dfOut$District_sampled, locsizes$loc)]
dfOut$minlocsize <- pmin(dfOut$locsize1,dfOut$locsize2)

### we assign colors based on the frequency of >0.50 IBD 

dfOut$segcolor <- "NA"
#dfOut$segcolor[dfOut$weight>0 & dfOut$weight<=0.5] <-  "ivory"
#dfOut$segcolor[dfOut$weight>0.5 & dfOut$weight<=1] <-  "seashell"
dfOut$segcolor[dfOut$weight>0 & dfOut$weight<=0.5] <-  "grey94"
dfOut$segcolor[dfOut$weight>0.5 & dfOut$weight<=1] <-  "grey85"
dfOut$segcolor[dfOut$weight>1 & dfOut$weight<=1.5] <- brewer.pal(9, "YlOrRd") [1]
dfOut$segcolor[dfOut$weight>1.5 & dfOut$weight<=2] <- brewer.pal(9, "YlOrRd") [2]
dfOut$segcolor[dfOut$weight>2 & dfOut$weight<=2.5] <- brewer.pal(9, "YlOrRd") [3]
dfOut$segcolor[dfOut$weight>2.5 & dfOut$weight<=3] <- brewer.pal(9, "YlOrRd") [4]
dfOut$segcolor[dfOut$weight>3 & dfOut$weight<=3.5] <- brewer.pal(9, "YlOrRd") [5]
dfOut$segcolor[dfOut$weight>3.5 & dfOut$weight<=4] <- brewer.pal(9, "YlOrRd") [6]
dfOut$segcolor[dfOut$weight>4 & dfOut$weight<=4.5] <- brewer.pal(9, "YlOrRd") [7]
dfOut$segcolor[dfOut$weight>4.5 & dfOut$weight<=5] <-  brewer.pal(9, "YlOrRd") [8]
dfOut$segcolor[dfOut$weight>5 & dfOut$weight<=6.5] <-  brewer.pal(9, "YlOrRd") [9]
dfOut$segcolor[dfOut$weight>6.5 & dfOut$weight<=8] <-  "darkviolet"
dfOut$segcolor[dfOut$weight>8 & dfOut$weight<=9.5] <-  "darkorchid4"
dfOut$segcolor[dfOut$weight>9.5 & dfOut$weight<=11] <-  "blue1"
dfOut$segcolor[dfOut$weight>11] <-  "midnightblue"

### we will only plot zone comparisons containing >50 observations

dfOut_g0_50comps <- dfOut %>% filter(num_comps > 50 & perc_clo > 0 & dfOut$minlocsize >0)  # keep only  positive flows 

### here are some other filter that we could optionally use

dfOut_nc_g50 <- dfOut %>% filter(num_comps > 50 & perc_clo > 0 & dfOut$minlocsize >1)  # keep only  positive flows 
dfOut_nc_g50_le <- dfOut %>% filter(num_comps > 50 & perc_clo > 0 & dfOut$minlocsize >1)  # keep only  positive flows 
dfOut_nc_g50_numcc_g1 <- dfOut_nc_g50[dfOut_nc_g50$V9>1 | dfOut_nc_g50$V10>1,]

### we make a separate dataframe for 'local' comparisons, also adding coordinates corresponding to each zone

dfLoc <- df_district %>% filter(type=="Local")
q75Loc <- quantile(dfLoc$perc_clo)[4]
dfLoc$x <- geos$lon[match(dfLoc$District_infected, geos$name)]
dfLoc$y <- geos$lat[match(dfLoc$District_infected, geos$name)]
dfLoc$xend <- geos$lon[match(dfLoc$District_sampled, geos$name)]
dfLoc$yend <- geos$lat[match(dfLoc$District_sampled, geos$name)]
dfLoc$from<- geos$id[match(dfLoc$District_infected, geos$name)]
dfLoc$to<- geos$id[match(dfLoc$District_sampled, geos$name)]
dfLoc$weight <- dfLoc$perc_clo
dfLoc$category <- dfLoc$from
dfLoc$maxcc <- pmax(dfLoc$V9,dfLoc$V10)

# cat groups_v2.txt | while read i;
# do awk '{print $1"\t"$12"xxx"$2"\t"$14}' myIBD_grouped_v2.txt | perl -pe 's/xxx/\n/g' | sort | uniq |
#   egrep -c "[[:space:]]$i$" | awk '{print "'$i'\t"$0}'; done > Pf_sampsizes_grouped_v2_myIBD.txt
# cut -f1 Pf_sampsizes_grouped_v2_myIBD.txt | while read i;
# do awk '$3 > 0.6935646 && $12=="'$i'" && $14=="'$i'" {print $1"xxx"$2}' myIBD_grouped_v2.txt |
#   perl -pe 's/xxx/\n/g' | sort | uniq | wc -l; done | paste Pf_sampsizes_grouped_v2_myIBD.txt - |
#   sed '1iloc\tnumsamps\tnumlocallyclonal' > Pf_sampsizes_grouped_v2_myIBD_with_locallyclonalnums.txt

dfLoc$locsize <- locsizes$numsamps[match(dfLoc$District_sampled, locsizes$loc)]
dfLoc$numclonal <- locsizes$numlocallyclonal[match(dfLoc$District_sampled, locsizes$loc)]
dfLoc$numunclonal <- dfLoc$locsize - dfLoc$numclonal
extrasites1 <- dfOut[dfOut$District_infected %in% locsizes$loc[locsizes$numsamps==1],]
extrasites2 <- dfOut[dfOut$District_sampled %in% locsizes$loc[locsizes$numsamps==1],]
extrasites2 <- extrasites2[!extrasites2$District_sampled %in% extrasites1$District_infected,]
dfLoc$perclotimeslocsize <- dfLoc$locsize * dfLoc$V5
dfLoc$perunclotimeslocsize <- dfLoc$locsize - dfLoc$perclotimeslocsize
dfLoc_nc_g50 <- dfLoc %>% filter(num_comps > 50)  # keep only  positive flows 

dfLoc$fillcolor <- "white"
dfLoc$fillcolor[dfLoc$weight>0 & dfLoc$weight<=0.5] <-  "grey94"
dfLoc$fillcolor[dfLoc$weight>0.5 & dfLoc$weight<=1] <-  "grey85"
dfLoc$fillcolor[dfLoc$weight>1 & dfLoc$weight<=1.5] <- brewer.pal(9, "YlOrRd") [1]
dfLoc$fillcolor[dfLoc$weight>1.5 & dfLoc$weight<=2] <- brewer.pal(9, "YlOrRd") [2]
dfLoc$fillcolor[dfLoc$weight>2 & dfLoc$weight<=2.5] <- brewer.pal(9, "YlOrRd") [3]
dfLoc$fillcolor[dfLoc$weight>2.5 & dfLoc$weight<=3] <- brewer.pal(9, "YlOrRd") [4]
dfLoc$fillcolor[dfLoc$weight>3 & dfLoc$weight<=3.5] <- brewer.pal(9, "YlOrRd") [5]
dfLoc$fillcolor[dfLoc$weight>3.5 & dfLoc$weight<=4] <- brewer.pal(9, "YlOrRd") [6]
dfLoc$fillcolor[dfLoc$weight>4 & dfLoc$weight<=4.5] <- brewer.pal(9, "YlOrRd") [7]
dfLoc$fillcolor[dfLoc$weight>4.5 & dfLoc$weight<=5] <-  brewer.pal(9, "YlOrRd") [8]
dfLoc$fillcolor[dfLoc$weight>5 & dfLoc$weight<=6.5] <-  brewer.pal(9, "YlOrRd") [9]
dfLoc$fillcolor[dfLoc$weight>6.5 & dfLoc$weight<=8] <-  "darkviolet"
dfLoc$fillcolor[dfLoc$weight>8 & dfLoc$weight<=9.5] <-  "darkorchid4"
dfLoc$fillcolor[dfLoc$weight>9.5 & dfLoc$weight<=11] <-  "blue1"
dfLoc$fillcolor[dfLoc$weight>11] <-  "midnightblue"

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

genmap3 <- ggplot(geos) + ggtitle("clonal connections") + 
  
  geom_polygon(data=guy1.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=bra0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=ven0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=sur0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  
    geom_curve(aes(x = x, y = y, xend = xend, yend = yend), size=4, color=dfOut_g0_50comps$segcolor,
             data = dfOut_g0_50comps, curvature = 0, alpha = 0.8) + scale_size(range=c(0.05,7)) +
  #geom_curve(aes(x = x, y = y, xend = xend, yend = yend), size=1, color = "black",
  #           data = dfOut_nc_g50_numcc_g1, curvature = 0, alpha = 1) +
  geom_point(aes(x = x, y = y),           # draw nodes
             shape = 23, size = 1.4, fill="black", data=extrasites1) +
  geom_point(aes(x = xend, y = yend),           # draw nodes
             shape = 23, size = 1.4, fill="black", data=extrasites2) +
  #geom_point(aes(x = x, y = y),           # draw nodes
  #shape = 21, fill=dfLoc$infectcol, stroke=0.5*sqrt(dfLoc$perc_clo), size = 8*(log10(dfLoc$locsize)), data=dfLoc) +
  coord_fixed(xlim = c(-61.5, -56.5), ylim = c(1.1, 8.35)) + maptheme 

genmap3 <- genmap3 + geom_point(aes(x = x, y = y, size=locsize), shape=21, fill=dfLoc$fillcolor, data=dfLoc) + scale_size(range=c(3,9))

genmap3

Pf_genmap3 <- genmap3 

### we can now create the same map for P. vivax, clearing workspace except for Pf_genmap3

rm(list=ls()[ls()!="Pf_genmap3"])

library(assertthat)
library(dplyr)
library(purrr)
library(igraph)
library(ggplot2)
library(ggmap)
require(gridExtra)
library(grid)
library(scatterpie)
library(tidyverse) 
library(RColorBrewer)    
library(patchwork)
library(ggrepel)
library(MMWRweek)
library(zoo)
library(circlize)
library(webshot2)
library(dplyr)

setwd("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana")

#### load file containing frequency of >0.50 pairwise relatedness within and between epidemiological zones (computed via script: Pv_g50.sh)

gen <- read.table("Pv_perc_g50_amoung_groups_v4_with_distances_and_indcount_and_cccount.txt", header=F)

### analysis will exclude Venezuelan infections

gen <- gen %>% filter(V1!="Venezuela" & V2!="Venezuela")
gen$geosample <- gen$V1
gen$geoinfect <- gen$V2
gen$geosample[is.na(gen$geosample)] <- ""
gen$geoinfect[is.na(gen$geoinfect)] <- ""

df <- gen
df$District_infected <- pmin(df$geoinfect,df$geosample)
df$District_sampled <- pmax(df$geoinfect,df$geosample)

### classify comparisons as within or between-zone (carrying over 'local' and 'outsourced' terminology from case flow mapping)

df_district <- df %>% mutate(type = case_when(District_infected==District_sampled ~ "Local",
                                              (District_infected=="" | District_sampled=="") ~ "!Missing", 
                                              (District_infected!=District_sampled & District_infected!="" & District_sampled!="") ~ "Outsourced"))

df_district$District_infected_gr <- gsub("$","_gr", df_district$District_infected)
df_district$District_sampled_gr <- gsub("$","_gr", df_district$District_sampled)

df_district <- df_district[df_district$type!="!Missing",]

### we will be plotting percent pairs exhibiting >0.50 IBD and use total number of comparisons to size plotted elements

df_district$perc_clo <- df_district$V5*100
df_district$num_comps <- df_district$V4
df_district$unclonal_comps <- df_district$V4-df_district$V3

set.seed(123)  # set random generator state for the same output

### this contains centroid coordinates for each of the epidemiological zones

geos <- read.table("centroids_via_COORDS_grouped_v4.txt", header=T)
geos$id <- rownames(geos)

### we make a separate dataframe for 'outsourced' comparisons, also adding coordinates corresponding to each zone

dfOut <- df_district %>% filter(type=="Outsourced")
q75Out <- quantile(dfOut$perc_clo)[4]
sdOut <- sd(dfOut$perc_clo)
dfOut$x <- geos$lon[match(dfOut$District_infected, geos$name)]
dfOut$y <- geos$lat[match(dfOut$District_infected, geos$name)]
dfOut$xend <- geos$lon[match(dfOut$District_sampled, geos$name)]
dfOut$yend <- geos$lat[match(dfOut$District_sampled, geos$name)]
dfOut$from<- geos$id[match(dfOut$District_infected, geos$name)]
dfOut$to<- geos$id[match(dfOut$District_sampled, geos$name)]
dfOut$weight <- dfOut$perc_clo
dfOut$category <- dfOut$from
dfOut$maxcc <- pmax(dfOut$V9,dfOut$V10)

locsizes <- read.table("Pv_sampsizes_grouped_v4_myIBD_with_locallyclonalnums.txt", header=T)
dfOut$locsize1 <- locsizes$numsamps[match(dfOut$District_infected, locsizes$loc)]
dfOut$locsize2 <- locsizes$numsamps[match(dfOut$District_sampled, locsizes$loc)]
dfOut$minlocsize <- pmin(dfOut$locsize1,dfOut$locsize2)

### we assign colors based on the frequency of >0.50 IBD 

dfOut$segcolor <- "NA"
#dfOut$segcolor[dfOut$weight>0 & dfOut$weight<=0.5] <-  "ivory"
#dfOut$segcolor[dfOut$weight>0.5 & dfOut$weight<=1] <-  "seashell"
dfOut$segcolor[dfOut$weight>0 & dfOut$weight<=0.5] <-  "grey94"
dfOut$segcolor[dfOut$weight>0.5 & dfOut$weight<=1] <-  "grey85"
dfOut$segcolor[dfOut$weight>1 & dfOut$weight<=1.5] <- brewer.pal(9, "YlOrRd") [1]
dfOut$segcolor[dfOut$weight>1.5 & dfOut$weight<=2] <- brewer.pal(9, "YlOrRd") [2]
dfOut$segcolor[dfOut$weight>2 & dfOut$weight<=2.5] <- brewer.pal(9, "YlOrRd") [3]
dfOut$segcolor[dfOut$weight>2.5 & dfOut$weight<=3] <- brewer.pal(9, "YlOrRd") [4]
dfOut$segcolor[dfOut$weight>3 & dfOut$weight<=3.5] <- brewer.pal(9, "YlOrRd") [5]
dfOut$segcolor[dfOut$weight>3.5 & dfOut$weight<=4] <- brewer.pal(9, "YlOrRd") [6]
dfOut$segcolor[dfOut$weight>4 & dfOut$weight<=4.5] <- brewer.pal(9, "YlOrRd") [7]
dfOut$segcolor[dfOut$weight>4.5 & dfOut$weight<=5] <-  brewer.pal(9, "YlOrRd") [8]
dfOut$segcolor[dfOut$weight>5 & dfOut$weight<=6.5] <-  brewer.pal(9, "YlOrRd") [9]
dfOut$segcolor[dfOut$weight>6.5 & dfOut$weight<=8] <-  "darkviolet"
dfOut$segcolor[dfOut$weight>8 & dfOut$weight<=9.5] <-  "darkorchid4"
dfOut$segcolor[dfOut$weight>9.5 & dfOut$weight<=11] <-  "blue1"
dfOut$segcolor[dfOut$weight>11] <-  "midnightblue"

### we will only plot zone comparisons containing >50 observations

dfOut_g0_50comps <- dfOut %>% filter(num_comps > 50 & perc_clo > 0 & dfOut$minlocsize >0)  # keep only  positive flows 

### here are some other filter that we could optionally use

dfOut_nc_g50 <- dfOut %>% filter(num_comps > 50 & perc_clo > 0 & dfOut$minlocsize >1)  # keep only  positive flows 
dfOut_nc_g50_le <- dfOut %>% filter(num_comps > 50 & perc_clo > 0 & dfOut$minlocsize >1)  # keep only  positive flows 
dfOut_nc_g50_numcc_g1 <- dfOut_nc_g50[dfOut_nc_g50$V9>1 | dfOut_nc_g50$V10>1,]

### we make a separate dataframe for 'local' comparisons, also adding coordinates corresponding to each zone

dfLoc <- df_district %>% filter(type=="Local")
q75Loc <- quantile(dfLoc$perc_clo)[4]
dfLoc$x <- geos$lon[match(dfLoc$District_infected, geos$name)]
dfLoc$y <- geos$lat[match(dfLoc$District_infected, geos$name)]
dfLoc$xend <- geos$lon[match(dfLoc$District_sampled, geos$name)]
dfLoc$yend <- geos$lat[match(dfLoc$District_sampled, geos$name)]
dfLoc$from<- geos$id[match(dfLoc$District_infected, geos$name)]
dfLoc$to<- geos$id[match(dfLoc$District_sampled, geos$name)]
dfLoc$weight <- dfLoc$perc_clo
dfLoc$category <- dfLoc$from
dfLoc$maxcc <- pmax(dfLoc$V9,dfLoc$V10)

# cat groups_v2.txt | while read i;
# do awk '{print $1"\t"$12"xxx"$2"\t"$14}' myIBD_grouped_v2.txt | perl -pe 's/xxx/\n/g' | sort | uniq |
#   egrep -c "[[:space:]]$i$" | awk '{print "'$i'\t"$0}'; done > Pv_sampsizes_grouped_v2_myIBD.txt
# cut -f1 Pv_sampsizes_grouped_v2_myIBD.txt | while read i;
# do awk '$3 > 0.6935646 && $12=="'$i'" && $14=="'$i'" {print $1"xxx"$2}' myIBD_grouped_v2.txt |
#   perl -pe 's/xxx/\n/g' | sort | uniq | wc -l; done | paste Pv_sampsizes_grouped_v2_myIBD.txt - |
#   sed '1iloc\tnumsamps\tnumlocallyclonal' > Pv_sampsizes_grouped_v2_myIBD_with_locallyclonalnums.txt

dfLoc$locsize <- locsizes$numsamps[match(dfLoc$District_sampled, locsizes$loc)]
dfLoc$numclonal <- locsizes$numlocallyclonal[match(dfLoc$District_sampled, locsizes$loc)]
dfLoc$numunclonal <- dfLoc$locsize - dfLoc$numclonal
extrasites1 <- dfOut[dfOut$District_infected %in% locsizes$loc[locsizes$numsamps==1],]
extrasites2 <- dfOut[dfOut$District_sampled %in% locsizes$loc[locsizes$numsamps==1],]
extrasites2 <- extrasites2[!extrasites2$District_sampled %in% extrasites1$District_infected,]
dfLoc$perclotimeslocsize <- dfLoc$locsize * dfLoc$V5
dfLoc$perunclotimeslocsize <- dfLoc$locsize - dfLoc$perclotimeslocsize
dfLoc_nc_g50 <- dfLoc %>% filter(num_comps > 50)  # keep only  positive flows 

dfLoc$fillcolor <- "white"
dfLoc$fillcolor[dfLoc$weight>0 & dfLoc$weight<=0.5] <-  "grey94"
dfLoc$fillcolor[dfLoc$weight>0.5 & dfLoc$weight<=1] <-  "grey85"
dfLoc$fillcolor[dfLoc$weight>1 & dfLoc$weight<=1.5] <- brewer.pal(9, "YlOrRd") [1]
dfLoc$fillcolor[dfLoc$weight>1.5 & dfLoc$weight<=2] <- brewer.pal(9, "YlOrRd") [2]
dfLoc$fillcolor[dfLoc$weight>2 & dfLoc$weight<=2.5] <- brewer.pal(9, "YlOrRd") [3]
dfLoc$fillcolor[dfLoc$weight>2.5 & dfLoc$weight<=3] <- brewer.pal(9, "YlOrRd") [4]
dfLoc$fillcolor[dfLoc$weight>3 & dfLoc$weight<=3.5] <- brewer.pal(9, "YlOrRd") [5]
dfLoc$fillcolor[dfLoc$weight>3.5 & dfLoc$weight<=4] <- brewer.pal(9, "YlOrRd") [6]
dfLoc$fillcolor[dfLoc$weight>4 & dfLoc$weight<=4.5] <- brewer.pal(9, "YlOrRd") [7]
dfLoc$fillcolor[dfLoc$weight>4.5 & dfLoc$weight<=5] <-  brewer.pal(9, "YlOrRd") [8]
dfLoc$fillcolor[dfLoc$weight>5 & dfLoc$weight<=6.5] <-  brewer.pal(9, "YlOrRd") [9]
dfLoc$fillcolor[dfLoc$weight>6.5 & dfLoc$weight<=8] <-  "darkviolet"
dfLoc$fillcolor[dfLoc$weight>8 & dfLoc$weight<=9.5] <-  "darkorchid4"
dfLoc$fillcolor[dfLoc$weight>9.5 & dfLoc$weight<=11] <-  "blue1"
dfLoc$fillcolor[dfLoc$weight>11] <-  "midnightblue"

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

genmap3 <- ggplot(geos) + ggtitle("clonal connections") + 
  
  geom_polygon(data=guy1.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=bra0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=ven0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  geom_polygon(data=sur0.spldf[["spdf"]], aes(x=long, y=lat, group=group), fill = "white", color = "black", size=.3) +
  
    geom_curve(aes(x = x, y = y, xend = xend, yend = yend), size=4, color=dfOut_g0_50comps$segcolor,
             data = dfOut_g0_50comps, curvature = 0, alpha = 0.8) + scale_size(range=c(0.05,7)) +
  #geom_curve(aes(x = x, y = y, xend = xend, yend = yend), size=1, color = "black",
  #           data = dfOut_nc_g50_numcc_g1, curvature = 0, alpha = 1) +
  geom_point(aes(x = x, y = y),           # draw nodes
             shape = 23, size = 1.4, fill="black", data=extrasites1) +
  geom_point(aes(x = xend, y = yend),           # draw nodes
             shape = 23, size = 1.4, fill="black", data=extrasites2) +
  #geom_point(aes(x = x, y = y),           # draw nodes
  #shape = 21, fill=dfLoc$infectcol, stroke=0.5*sqrt(dfLoc$perc_clo), size = 8*(log10(dfLoc$locsize)), data=dfLoc) +
  coord_fixed(xlim = c(-61.5, -56.5), ylim = c(1.1, 8.35)) + maptheme 

genmap3 <- genmap3 + geom_point(aes(x = x, y = y, size=locsize), shape=21, fill=dfLoc$fillcolor, data=dfLoc) + scale_size(range=c(3,9))

genmap3

Pv_genmap3 <- genmap3
