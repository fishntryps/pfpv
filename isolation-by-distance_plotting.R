
##################################################################################################################################
### this script creates isolation-by-distance plots for P. falciparum and P. vivax using various measures of relatedness (IBD) ###
##################################################################################################################################

### run time is seconds to minutes

### load packages and set working directory in R 4.2.2

library(mosaic)
setwd("C:/Users/philipp/OneDrive - Harvard University/Documents/Harvard/Pf_Guyana")

### set up 3 panel plot 

par(mfrow=c(1,3))
par(mar=c(5,5,2,2))

### load P.falciparum file containing fraction IBD, spatial distances, and Julian-type sampling dates for each sample pair
### either use all pairs or exclude pairs involving superclones

pf <- read.table("Pf_test_for_R_only_p025_mono_no_ge10clones_noninit.txt", header=F)

# pf <- read.table("Pf_test_for_R_only_p025_mono_noninit.txt", header=F)

### optionally exclude infections from Venezuela

# pf <- pf[!grepl("^CEM|Venez|^SPT|^PW",pf$V1) & !grepl("^CEM|Venez|^SPT|^PW",pf$V2),]

### we only need the columns representing fraction IBD, spatial distance in meters, and sampling dates 

pf <- pf[,c(3,10,11,12)]

### spatial distance in km

pf$V10 <- pf$V10/1000

### temporal distance in days (not used in this script but kept in case we wanted to assess IBD relationships to temporal distance instead of spatial distance)

pf$daydist <- sqrt((pf$V11 - pf$V12)^2)

colnames(pf) <- c("ibd","geodist","d1","d2","daydist")

### create a function to calculate fraction sammple pairs with IBD above 0.5 within sliding window of spatial distance

starts <- seq(0, max(pf$geodist), step)
Fun <- function(st, width, pf, fun) {fun(subset(pf, geodist >= st & geodist <= st + width - 1)$ibd)}
perc_clon <- function(x) length(x[x>0.5])/length(x)
width <- 31
step <- 10
pf_plot <- data.frame(starts, ends = starts + width - 1, perc_clon = sapply(starts, Fun, width, pf, perc_clon), N = sapply(starts, Fun, width, pf, length))

### this part just creates an empty plot that we will use as a background for the upcoming plot loop

plot(pf_plot$starts,pf_plot$perc_clon, cex=0, xlim=c(0,300), ylim=c(0,0.20), las=1, ylab="fraction comparisons", cex.lab=0.9, xlab="spatial distance (km)")

### this instructs next plot to use the same panel as before

par(new=T)

### do the same except looping across different IBD thresholds and plotting all in same left panel plot

k=0
pflist <- c()
  
for (i in c(0.375, 0.4, 0.45, 0.5, 0.6, 0.7,0.9)){
  
  k=k+1
  starts <- seq(1, max(pf$geodist), step)
  perc_clon <- function(x) length(x[x>=i])/length(x)
  pf_plot <- data.frame(starts, ends = starts + width - 1, perc_clon = sapply(starts, Fun, width, pf, perc_clon), N = sapply(starts, Fun, width, pf, length))
  
  ### exclude large distance classes where sample sizes become unreliable

  pf_plot <- pf_plot[pf_plot$starts<=291,]
  
  pf_plot$q05 <- NA
  pf_plot$q95 <- NA

### this part additionally generates confidence intervals by bootstrapping

  for (row in 1:dim(pf_plot)[1]){
    sub <- pf_plot[row,]
    myprob <- sub$perc_clon
    myN <- sub$N
    boot_dsn <- do(1000) * rflip(myN, prob=myprob)
    pf_plot[row,"q05"] <- quantile(boot_dsn$prop, probs=0.05, na.rm=T)
    pf_plot[row,"q95"] <- quantile(boot_dsn$prop, probs=0.95, na.rm=T)
  }

  colors <- c(
  "#fdd7dc",
  "#f8aeb9",
  "#ef8598",
  "#dc3c68",
  "#b80d4c",
  "#830634",
  "#3a0113")
  
  for (row in 1:dim(pf_plot)[1]){
    sub <- pf_plot[row,]
    points(x=sub$starts, y=sub$q05, pch="-",col=colors[k])
    points(x=sub$starts, y=sub$q95, pch="-",col=colors[k])
    segments(sub$starts, sub$q05, x1 = sub$starts, y1 = sub$q95, col=colors[k])}
  par(new=T)
  plot(pf_plot$starts,pf_plot$perc_clon, type='l', col=colors[k], cex=1.2, xlim=c(0,300), ylim=c(0,0.20), las=1, ylab="", cex.lab=1.5, xlab="", xaxt='n', yaxt='n')
  par(new=T)
  plot(pf_plot$starts,pf_plot$perc_clon, pch=21, bg=colors[k], lwd=0.02, cex=1.2, xlim=c(0,300), ylim=c(0,0.20), las=1, ylab="", cex.lab=1.5, xlab="", xaxt='n', yaxt='n')
  par(new=T)
  pf_lm <- summary(lm(pf_plot$perc_clon[pf_plot$starts<75] ~ pf_plot$starts[pf_plot$starts<75]))
  pflist <- c(pflist, pf_lm$coefficients[2])
}

par(new=F)

### now for the next panel (middle) repeat all of the above for P. vivax, optionally excluding infections from Venezuela or infections representing relapses, etc.

pv <- read.table("Pv_test_for_R_only_p025_mono_noninit.txt", header=F)

pv <- pv[!grepl("^CEM|Venez|^SPT|^PW",pv$V1) & !grepl("^CEM|Venez|^SPT|^PW",pv$V2),]

# relapses <- c("G4G1013",
#              "G4G1019",
#              "G4G1040",
#              "G4G1426",
#              "G4G691",
#              "G4G697",
#              "G4G454",
#              "G7B582",
#              "G4G424",
#              "G7B200",
#              "G7B180")

pv <- pv %>% filter(!V1 %in% relapses & !V2 %in% relapses) 

pv <- pv[,c(3,10,11,12)]
pv$V10 <- pv$V10/1000
pv$daydist <- sqrt((pv$V11 - pv$V12)^2)
colnames(pv) <- c("ibd","geodist","d1","d2","daydist")
pv <- pv[!is.na(pv$geodist),]
Fun <- function(st, width, pv, fun) {fun(subset(pv, geodist >= st & geodist <= st + width - 1)$ibd)}
width <- 31
step <- 10

pv_plot <- data.frame(starts, ends = starts + width - 1, perc_clon = sapply(starts, Fun, width, pv, perc_clon), N = sapply(starts, Fun, width, pv, length))

plot(pv_plot$starts,pv_plot$perc_clon, cex=0, xlim=c(0,300), ylim=c(0,0.004), las=1, ylab="fraction comparisons\n", cex.lab=0.9, xlab="spatial distance (km)")

par(new=T)

k=0
pvlist <- c()

for (i in c(0.375, 0.4, 0.45, 0.5, 0.6, 0.7,0.9)){
  
  k=k+1
  starts <- seq(1, max(pv$geodist), step)
  perc_clon <- function(x) length(x[x>=i])/length(x)
  pv_plot <- data.frame(starts, ends = starts + width - 1, perc_clon = sapply(starts, Fun, width, pv, perc_clon), N = sapply(starts, Fun, width, pv, length))
  
  pv_plot <- pv_plot[pv_plot$starts<=291,]
  
  pv_plot$q05 <- NA
  pv_plot$q95 <- NA
  for (row in 1:dim(pv_plot)[1]){
    sub <- pv_plot[row,]
    myprob <- sub$perc_clon
    myN <- sub$N
    boot_dsn <- do(1000) * rflip(myN, prob=myprob)
    pv_plot[row,"q05"] <- quantile(boot_dsn$prop, probs=0.05, na.rm=T)
    pv_plot[row,"q95"] <- quantile(boot_dsn$prop, probs=0.95, na.rm=T)
  }

  colors <- c(
    "#d5e8ff",
    "#97c5ff",
    "#6cadff",
    "#3c93ff",
    "#1573de",
    "#0c509e",
    "#022148")
  
  for (row in 1:dim(pv_plot)[1]){
    sub <- pv_plot[row,]
    points(x=sub$starts, y=sub$q05, pch="-",col=colors[k])
    points(x=sub$starts, y=sub$q95, pch="-",col=colors[k])
    segments(sub$starts, sub$q05, x1 = sub$starts, y1 = sub$q95, col=colors[k])}
  par(new=T)
  plot(pv_plot$starts,pv_plot$perc_clon, type='l', col=colors[k], cex=1.2, xlim=c(0,300), ylim=c(0,0.004), las=1, ylab="", cex.lab=1.5, xlab="", xaxt='n', yaxt='n')
  par(new=T)
  plot(pv_plot$starts,pv_plot$perc_clon, pch=21, bg=colors[k], , lwd=0.02, cex=1.2, xlim=c(0,300), ylim=c(0,0.004), las=1, ylab="", cex.lab=1.5, xlab="", xaxt='n', yaxt='n')
  par(new=T)
  pv_lm <- summary(lm(pv_plot$perc_clon[pv_plot$starts<75] ~ pv_plot$starts[pv_plot$starts<75]))
  pvlist <- c(pvlist, pv_lm$coefficients[2])
}

### optionally make a color key

# rect(xleft=215, xright=225, ytop=0.0036, ybottom=0.0035, col="#d5e8ff")
# rect(xleft=215, xright=225, ytop=0.0035, ybottom=0.0034, col="#97c5ff")
# rect(xleft=215, xright=225, ytop=0.0034, ybottom=0.0033, col="#6cadff")
# rect(xleft=215, xright=225, ytop=0.0033, ybottom=0.0032, col="#3c93ff")
# rect(xleft=215, xright=225, ytop=0.0032, ybottom=0.0031, col="#1573de")
# rect(xleft=215, xright=225, ytop=0.0031, ybottom=0.0030, col="#0c509e")
# rect(xleft=215, xright=225, ytop=0.0029, ybottom=0.0028, col="#022148")
# 
# rect(xleft=200, xright=210, ytop=0.0036, ybottom=0.0035, col="#fdd7dc")
# rect(xleft=200, xright=210, ytop=0.0035, ybottom=0.0034, col="#f8aeb9")
# rect(xleft=200, xright=210, ytop=0.0034, ybottom=0.0033, col="#ef8598")
# rect(xleft=200, xright=210, ytop=0.0033, ybottom=0.0032, col="#dc3c68")
# rect(xleft=200, xright=210, ytop=0.0032, ybottom=0.0031, col="#b80d4c")
# rect(xleft=200, xright=210, ytop=0.0031, ybottom=0.0030, col="#830634")
# rect(xleft=200, xright=210, ytop=0.0029, ybottom=0.0028, col="#3a0113")


### repeat the above process but calculating mean absolute IBD instead of fraction sammple pairs above IBD thresholds

width= 31
exclrelapse <- "Yes"
exclPflarge <- "No"
step <- 10

pf <- read.table("Pf_test_for_R_only_p025_mono_no_ge10clones_noninit.txt", header=F)

# pf <- read.table("Pf_test_for_R_only_p025_mono_noninit.txt", header=F)

### optionally exclude infections from Venezuela

# pf <- pf[!grepl("^CEM|Venez|^SPT|^PW",pf$V1) & !grepl("^CEM|Venez|^SPT|^PW",pf$V2),]

pf <- pf[,c(3,10,11,12)]
pf$V10 <- pf$V10/1000
pf <- pf[pf$V10<1000,]
pf$daydist <- sqrt((pf$V11 - pf$V12)^2)
colnames(pf) <- c("ibd","geodist","d1","d2","daydist")
pf <- pf[!is.na(pf$geodist),]

Fun <- function(st, width, pf, fun) {fun(subset(pf, geodist >= st & geodist <= st + width - 1)$ibd)}

starts <- seq(0, max(pf$geodist), step)

MEAN_ibd <- function(x) mean(x)

pf_plot <- data.frame(starts, ends = starts + width - 1, MEAN_ibd = sapply(starts, Fun, width, pf, MEAN_ibd), N = sapply(starts, Fun, width, pf, length))

pf_plot

pf_plot$q05 <- NA
pf_plot$q95 <- NA
for (row in 1:dim(pf_plot)[1]){
  sub <- pf_plot[row,]
  myprob <- sub$MEAN_ibd
  myN <- sub$N
  boot_dsn <- do(1000) * rflip(myN, prob=myprob)
  pf_plot[row,"q05"] <- quantile(boot_dsn$prop, probs=0.05, na.rm=T)
  pf_plot[row,"q95"] <- quantile(boot_dsn$prop, probs=0.95, na.rm=T)
}

pv <- read.table("Pv_test_for_R_only_p025_mono_noninit.txt", header=F)

pv <- pv[!grepl("^CEM|Venez|^SPT|^PW",pv$V1) & !grepl("^CEM|Venez|^SPT|^PW",pv$V2),]

# relapses <- c("G4G1013",
#              "G4G1019",
#              "G4G1040",
#              "G4G1426",
#              "G4G691",
#              "G4G697",
#              "G4G454",
#              "G7B582",
#              "G4G424",
#              "G7B200",
#              "G7B180")


# pv <- pv %>% filter(!V1 %in% relapses & !V2 %in% relapses) 


pv <- pv[,c(3,10,11,12)]
pv$V10 <- pv$V10/1000
pv <- pv[pv$V10<1000,]
pv$daydist <- sqrt((pv$V11 - pv$V12)^2)
colnames(pv) <- c("ibd","geodist","d1","d2","daydist")
pv <- pv[!is.na(pv$geodist),]
Fun <- function(st, width, pv, fun) {fun(subset(pv, geodist >= st & geodist <= st + width - 1)$ibd)}

length(pv$ibd[pv$ibd>0.1155])/length(pv$ibd)

starts <- seq(0, max(pv$geodist), step)
MEAN_ibd <- function(x) mean(x)
pv_plot <- data.frame(starts, ends = starts + width - 1, MEAN_ibd = sapply(starts, Fun, width, pv, MEAN_ibd), N = sapply(starts, Fun, width, pv, length))

pf_plot1 <- pf_plot[pf_plot$starts<=291,]
pv_plot1 <- pv_plot[pv_plot$starts<=291,]

### and finally do the same except calculating standard deviation so we can add error bars

rm(list=ls()[!ls() %in% c("pf_plot1", "pv_plot1")])

width= 31
exclrelapse <- "Yes"
exclPflarge <- "No"
step <- 10

pf <- read.table("Pf_test_for_R_only_p025_mono_no_ge10clones_noninit.txt", header=F)

# pf <- read.table("Pf_test_for_R_only_p025_mono_noninit.txt", header=F)

### optionally exclude infections from Venezuela

# pf <- pf[!grepl("^CEM|Venez|^SPT|^PW",pf$V1) & !grepl("^CEM|Venez|^SPT|^PW",pf$V2),]

pf <- pf[,c(3,10,11,12)]
pf$V10 <- pf$V10/1000
pf <- pf[pf$V10<1000,]
pf$daydist <- sqrt((pf$V11 - pf$V12)^2)
colnames(pf) <- c("ibd","geodist","d1","d2","daydist")
pf <- pf[!is.na(pf$geodist),]

Fun <- function(st, width, pf, fun) {fun(subset(pf, geodist >= st & geodist <= st + width - 1)$ibd)}

starts <- seq(0, max(pf$geodist), step)

SD_ibd <- function(x) sd(x)

pf_plot <- data.frame(starts, ends = starts + width - 1, SD_ibd = sapply(starts, Fun, width, pf, SD_ibd), N = sapply(starts, Fun, width, pf, length))

pf_plot

pf_plot$q05 <- NA
pf_plot$q95 <- NA
for (row in 1:dim(pf_plot)[1]){
  sub <- pf_plot[row,]
  myprob <- sub$SD_ibd
  myN <- sub$N
  boot_dsn <- do(1000) * rflip(myN, prob=myprob)
  pf_plot[row,"q05"] <- quantile(boot_dsn$prop, probs=0.05, na.rm=T)
  pf_plot[row,"q95"] <- quantile(boot_dsn$prop, probs=0.95, na.rm=T)
}

pv <- read.table("Pv_test_for_R_only_p025_mono_noninit.txt", header=F)

pv <- pv[!grepl("^CEM|Venez|^SPT|^PW",pv$V1) & !grepl("^CEM|Venez|^SPT|^PW",pv$V2),]

# relapses <- c("G4G1013",
#              "G4G1019",
#              "G4G1040",
#              "G4G1426",
#              "G4G691",
#              "G4G697",
#              "G4G454",
#              "G7B582",
#              "G4G424",
#              "G7B200",
#              "G7B180")


# pv <- pv %>% filter(!V1 %in% relapses & !V2 %in% relapses) 

pv <- pv[,c(3,10,11,12)]
pv$V10 <- pv$V10/1000
pv <- pv[pv$V10<1000,]
pv$daydist <- sqrt((pv$V11 - pv$V12)^2)
colnames(pv) <- c("ibd","geodist","d1","d2","daydist")
pv <- pv[!is.na(pv$geodist),]
Fun <- function(st, width, pv, fun) {fun(subset(pv, geodist >= st & geodist <= st + width - 1)$ibd)}

length(pv$ibd[pv$ibd>0.1155])/length(pv$ibd)

starts <- seq(0, max(pv$geodist), step)
SD_ibd <- function(x) sd(x)
pv_plot <- data.frame(starts, ends = starts + width - 1, SD_ibd = sapply(starts, Fun, width, pv, SD_ibd), N = sapply(starts, Fun, width, pv, length))

pf_plot2 <- pf_plot[pf_plot$starts<=291,]
pv_plot2 <- pv_plot[pv_plot$starts<=291,]

### add windowed sd to table containing windowed means

pf_plot1$plusstdev <- pf_plot1$MEAN_ibd + pf_plot2$SD_ibd
pv_plot1$plusstdev <- pv_plot1$MEAN_ibd + pv_plot2$SD_ibd
pf_plot1$minusstdev <- pf_plot1$MEAN_ibd - pf_plot2$SD_ibd
pv_plot1$minusstdev <- pv_plot1$MEAN_ibd - pv_plot2$SD_ibd

### do not allow mean minus sd to pass below zero

pv_plot1$minusstdev[pv_plot1$minusstdev<0] <- 0 


### finally plot the last panel showing mean and sd IBD windowed over spatial distance

par(new=F)

plot(pf_plot1$starts,pf_plot1$perc_clon, pch=21, bg=rgb(
  col2rgb("#D41159")[1]/255,
  col2rgb("#D41159")[2]/255,
  col2rgb("#D41159")[3]/255, alpha=1), cex=0, xlim=c(0,300), ylim=c(0,0.45), las=1, ylab="", cex.lab=1.5, xlab="", yaxt='n')

par(new=T)

plot(xaxt='n', pf_plot1$starts,pf_plot1$MEAN_ibd, pch=21, bg=rgb(
  col2rgb("#D41159")[1]/255,
  col2rgb("#D41159")[2]/255,
  col2rgb("#D41159")[3]/255, alpha=1), cex=0, xlim=c(0,300), ylim=c(0,0.45), las=1, ylab="Mean IBD\n", cex.lab=1, xlab="Spatial distance (km)", yaxt='n')

for (row in 1:dim(pf_plot1)[1]){
  sub <- pf_plot1[row,]
  segments(sub$starts, sub$minusstdev, x1 = sub$starts, y1 = sub$plusstdev, lty="dashed", col=rgb(0.5,0.5,0.5, alpha=0.6))
  segments(sub$starts-1, sub$minusstdev, x1 = sub$starts+1, y1 = sub$minusstdev, col=rgb(0.5,0.5,0.5, alpha=0.6))
  segments(sub$starts-1, sub$plusstdev, x1 = sub$starts+1, y1 = sub$plusstdev, col=rgb(0.5,0.5,0.5, alpha=0.6))}


axis(2, at=c(0.00, 0.10, 0.20, 0.30, 0.40), labels=c("0.00", "0.10", "0.20", "0.30", "0.40"), las=1)

par(new=T)

plot(xaxt='n', pf_plot1$starts,pf_plot1$MEAN_ibd, pch=21, bg=rgb(
  col2rgb("#D41159")[1]/255,
  col2rgb("#D41159")[2]/255,
  col2rgb("#D41159")[3]/255, alpha=1), lwd = 0.02, cex=0.8, xlim=c(0,300), ylim=c(0,0.45), las=1, ylab="", cex.lab=2.5, xlab="", xaxt='n', yaxt='n')

for (row in 1:dim(pv_plot1)[1]){
  sub <- pv_plot1[row,]
  #points(x=sub$starts, y=sub$q05, pch="-")
  #points(x=sub$starts, y=sub$q95, pch="-")
  segments(sub$starts, sub$minusstdev, x1 = sub$starts, y1 = sub$plusstdev, lty="dashed", col=rgb(0.5,0.5,0.5, alpha=0.6))
  segments(sub$starts-1, sub$minusstdev, x1 = sub$starts+1, y1 = sub$minusstdev, col=rgb(0.5,0.5,0.5, alpha=0.6))
  segments(sub$starts-1, sub$plusstdev, x1 = sub$starts+1, y1 = sub$plusstdev, col=rgb(0.5,0.5,0.5, alpha=0.6))}
