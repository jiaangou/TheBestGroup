### New GAMs for smoothing over Gall Lake SS isotope data

library(mgcv)

IS <- read.csv(file.choose(), header=TRUE)

## N15

N15<-gam(N15~s(Date, bs="tp"), data = IS, method = "REML")

summary(N15) ##significance is good
plot(N15) ## looks good 
gam.check(N15)## residuals look ok  
qq.gam(N15) ## nice 

newN15 <- with(IS, data.frame(Date = seq(min(Date), max(Date),
                                                length.out = 200)))
newN15 <- cbind(newN15,
                    data.frame(predict(N15, newN15, se.fit = TRUE))) 
newN15 <- transform(newN15, 
                    upper = fit + (2 * se.fit),
                    lower = fit - (2 * se.fit))



plot(IS$Date, IS$N15, col="black", bg="blue4", pch=21,
     xlab="Depth (m)", ylab="N15")
lines(newN15$Date, newN15$fit, lwd=2)
lines(lty=3, newN15$Date, newN15$upper)
lines(lty=3, newN15$Date, newN15$lower)


#######################################
#### C13

C13<-gam(C13~s(Date, bs="tp"), data = IS, method = "REML")

summary(C13) ##significance is good
plot(C13) ## looks good 
gam.check(C13)## residuals look ok  
qq.gam(C13) ## nice 

newC13 <- with(IS, data.frame(Date = seq(min(Date), max(Date),
                                         length.out = 200)))
newC13 <- cbind(newC13,
                data.frame(predict(C13, newC13, se.fit = TRUE))) 
newC13 <- transform(newC13, 
                    upper = fit + (2 * se.fit),
                    lower = fit - (2 * se.fit))



plot(IS$Date, IS$C13, col="black", bg="green4", pch=21,
     xlab="Depth (m)", ylab="C13")
lines(newC13$Date, newC13$fit, lwd=2)
lines(lty=3, newC13$Date, newC13$upper)
lines(lty=3, newC13$Date, newC13$lower)


#######################################
###  CN

CN<-gam(CN~s(Date, bs="tp"), data = IS, method = "REML")

summary(CN) 
plot(CN) 
gam.check(CN) 
qq.gam(CN)

newCN <- with(IS, data.frame(Date = seq(min(Date), max(Date),
                                         length.out = 200)))
newCN <- cbind(newCN,
                data.frame(predict(CN, newCN, se.fit = TRUE))) 
newCN <- transform(newCN, 
                    upper = fit + (2 * se.fit),
                    lower = fit - (2 * se.fit))


###########################################333
############ %C and N
## add in inversion for gamma errors
C<-gam(C~s(Date, bs="tp"), data = IS, family=Gamma(link="log"), method = "REML")

summary(C) 
plot(C) 
gam.check(C) 
qq.gam(C)

fam<-family(C)
fam
str(fam) 

ilink <- fam$linkinv
ilink
Gamma()$linkinv 
ilink <- family(C)$linkinv 

newC<- with(IS, data.frame(Date = seq(min(Date), max(Date),
                                               length.out = 200)))
newC<- cbind(newC,
                    data.frame(predict(C, newC, se.fit = TRUE))) 
newC <- transform(newC, fit = ilink(fit), 
                        upper = ilink(fit + (2 * se.fit)), lower = ilink(fit - (2 * se.fit)))
################################
N<-gam(N~s(Date, bs="tp"), data = IS, family=Gamma(link="log"), method = "REML")

summary(N) 
plot(N) 
gam.check(N) 
qq.gam(N)

fam<-family(N)
fam
str(fam) 

ilink <- fam$linkinv
ilink
Gamma()$linkinv 
ilink <- family(N)$linkinv 

newN<- with(IS, data.frame(Date = seq(min(Date), max(Date),
                                      length.out = 200)))
newN<- cbind(newN,
             data.frame(predict(N, newN, se.fit = TRUE))) 
newN<- transform(newN, fit = ilink(fit), 
                  upper = ilink(fit + (2 * se.fit)), lower = ilink(fit - (2 * se.fit)))
###### plots
plot(IS$Date, IS$CN, col="black", bg="white", pch=21, cex=1.5,
     xlab="Depth (m)", ylab="C:N ratio")
lines(newCN$Date, newCN$fit, lwd=2)
lines(lty=3, newCN$Date, newCN$upper)
lines(lty=3, newCN$Date, newCN$lower)

### %N
plot(IS$Date, IS$N, col="black", bg="white", pch=21, cex=1.5, xlab="Depth (m)", ylab="%N",
     ylim=c(0,2.5))
lines(newN$Date, newN$fit, lwd=2)
lines(lty=3, newN$Date, newN$upper)
lines(lty=3, newN$Date, newN$lower)

### %C
plot(IS$Date, IS$C, col="black", bg="white", pch=21, cex=1.5, xlab="Depth (m)", ylab="%C",
     ylim=c(0,25))
lines(newC$Date, newC$fit, lwd=2)
lines(lty=3, newC$Date, newC$upper)
lines(lty=3, newC$Date, newC$lower)

### N15
plot(IS$Date, IS$N15, col="black", bg="white", pch=21, cex=1.5,
     xlab="Depth (m)", ylab="N15")
lines(newN15$Date, newN15$fit, lwd=2)
lines(lty=3, newN15$Date, newN15$upper)
lines(lty=3, newN15$Date, newN15$lower)

### C13
plot(IS$Date, IS$C13, col="black", bg="white", pch=21, cex=1.5,
     xlab="Depth (m)", ylab="C13")
lines(newC13$Date, newC13$fit, lwd=2)
lines(lty=3, newC13$Date, newC13$upper)
lines(lty=3, newC13$Date, newC13$lower)


### Now make a cluster
library(rioja)

iso <- read.csv(file.choose(), header=TRUE, row.names=1)
diss <- dist((iso))
clust <- chclust(diss, method = "coniss")
plot(clust, hang=-1)


#### Gonna redo PCAs of pigments (SQRT) and diatoms without transect numbers to better show changes
### Also gonna do a PCA of isotope data

### then seeing how they all look, put all three data types together in a 4th PCA

library(vegan)

### SQRT pigments first

gallpigments <- read.csv(file.choose(), header = TRUE, row.names=1)

ppca <- rda(gallpigments)
biplot(ppca)

summary(ppca)

zone <- read.csv(file.choose(), header = TRUE) ##groups

plot(ppca, type = c("points")) 

colvec <- c("black", "black", "black") ##colours
symb <- c(21, 21, 22) ##symbols - circle, square, triangle
back <- c("black", "white", "black") ## fill - black circle, open square, black triangle

with(zone, levels(Zone)) ##groups names
[1] "B1", "B2", "P"

plot(ppca, type = c("points"))
with(zone, 
     points(ppca, display = "sites", col = colvec[Zone], pch = symb[Zone], cex= 1.5, bg = back[Zone]))

with(zone, legend("topleft", legend = levels(Zone), bty = "n", col = colvec, inset = c(0.05,0), pch = symb, pt.bg = back, cex=1.5))

text(ppca, display = "species", cex = 0.8, col = "black")



### Now SQRT 5% diatoms

galldiatoms <- read.csv(file.choose(), header = TRUE, row.names=1)

dpca <- rda(galldiatoms)
biplot(dpca)

summary(dpca)

dsp <- scores(dpca, display ="species")
dsp

zone2 <- read.csv(file.choose(), header = TRUE) ##groups

plot(dpca, type = c("points")) 

colvec <- c("black", "black", "black") ##colours
symb <- c(21, 21, 22) ##symbols - circle, square, triangle
back <- c("black", "white", "black") ## fill - black circle, open square, black triangle

with(zone2, levels(Zone)) ##groups names
[1] "B1", "B2", "P"

plot(dpca, type = c("points"))
with(zone2, 
     points(dpca, display = "sites", col = colvec[Zone], pch = symb[Zone], cex= 1.5, bg = back[Zone]))

with(zone2, legend("topleft", legend = levels(Zone), bty = "n", col = colvec, inset = c(0.05,0), pch = symb, pt.bg = back, cex=1.5))

text(dpca, display = "species", cex = 0.8, col = "black")


### Isotopes

gallisotopes <- read.csv(file.choose(), header = TRUE, row.names=1)

Ipca <- rda(gallisotopes)
biplot(Ipca)

summary(Ipca)
### too pulled by high C values. SQRT transform data? 

zone3 <- read.csv(file.choose(), header = TRUE)

colvec <- c("black", "black", "black") ##colours
symb <- c(21, 21, 22) ##symbols
back <- c("black", "white", "black") ## fill 

with(zone3, levels(Zone)) 
[1]"B1", "B2", "P"


plot(Ipca, type=c("points"), scaling=2)
with(zone3,
     points(Ipca, display = "sites", col = colvec[Zone], pch = symb[Zone], cex= 1.5, bg = back[Zone], 
            scaling=2))
with(zone3, legend("bottom", legend = levels(Zone), 
                   bty = "n", col = colvec, inset = c(0.05,0), pch = symb, pt.bg = back, cex=1.5))

text(Ipca, display = "species", cex =1, col = "black")

########## species scores

dsp <- scores(dpca, display ="species")
dsp

psp <- scores(ppca, display ="species")
psp

isp <- scores(Ipca, display="species")
isp

write.csv(dsp, file = "galldiatspecies.csv")
write.csv(psp, file = "gallpigspecies.csv")
write.csv(isp, file = "gallisospecies.csv")

###### Combined data PCA

gallcombined <- read.csv(file.choose(), header = TRUE, row.names=1)

Cpca <- rda(gallcombined)
biplot(Cpca, sc=3, col="black")
text(Cpca, display="sites", sc=3)

summary(Cpca)

Csp <-scores(Cpca, choices =c(1,2,3), display="species")
Csp

write.csv(Csp, file = "gallcombinedpca.csv")

