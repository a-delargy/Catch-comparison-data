
#### Adam Delargy 

# size structured catch comparison rates

# 17/06/2021

library(selfisher)
library(splines)	# for bs function
library(plyr)	# for data aggregation
library(bbmle)	# for AICtab
library(snow)	# for bootstrapping across cores

#

# function for global legend

legend.fun<-function(...){
  opar<-par(fig=c(0,1,0,1),oma=rep(0,4),
            mar=rep(0,4), new=TRUE)
  on.exit(par(opar))
  plot(0,0,type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#


### load data 

setwd("C:/Users/osp813/OneDrive - Bangor University/Documents/EMFF/Survey data/Vessel comparison 2021")

dat <- read.csv("Comparison data for models.csv", header=TRUE)

dat$Station <- as.factor(dat$Station)

# divide cpuetotal by 10 for better numbers 
# lowest non-zero was 60 before, so 6 after this. 
dat$cpuetotalN <- dat$cpuetotal/10

##### fit models and store predictions

# loop through multiple comparisons

# the loop takes a while, so users may prefer to load the pre-made loop results after instead


newstore <- c()			# storage
plotstore <- c()		# storage 

for(i in 1:length(unique(dat$comp2))){			# loop through each vessel comparison combination

	set_dat <- dat[which(dat$comp2==unique(dat$comp2)[i]),]	# subset to a single vessel comparison
	set_dat <- droplevels(set_dat)

	# scale the scallop length for numerical stability 
	set_dat$sl <- scale(set_dat$Length

	# fit polynomial (Holst and Revill 2009)
	mod1 <- selfisher(prop~(sl+I(sl^2)+I(sl^3)), total=cpuetotalN, qratio=qratio, haul=Station,data=set_dat) 

	# fit splines like Miller (2013)
	mod2 <- selfisher(prop~ bs(sl, df=3), total=cpuetotalN, qratio=qratio, haul=Station,data=set_dat) 
	#mod3 <- selfisher(prop~ bs(sl, df=5), total=cpuetotalN, qratio=qratio, haul=Station,data=set_dat) 

	# compare models 
	fmod <- attr(AICtab(mod1,mod2),"row.names")[1]			# pick model based on lowest AIC
	#fmod <- attr(AICtab(mod1,mod2,mod3),"row.names")[1]		

	newdata <- expand.grid(sl = (unique(set_dat$Length)-mean(set_dat$Length))/sqrt(var(set_dat$Length)), cpuetotalN=1, Station=NA, qratio=1)	# create data.frame to make predictions
	newdata$prop <- predict(get(paste(fmod)), newdata=newdata, type="response")		# make predictions

	# bootstrap for confidence intervals 

	ncpus = 4	# need to adjust this based on how many cores your laptop has 
			# the following spreads the work across your cores
	cl = makeCluster(rep("localhost", ncpus))
	clusterExport(cl, "newdata")
	bs = bootSel(get(paste(fmod)), nsim = 1000, parallel = "snow", cl = cl,FUN = function(mod){predict(mod, newdata = newdata, type = "response")})
	stopCluster(cl)

	quants = apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975))		# get confidence intervals from bootstrap results
	newdata[,c("lo", "mid", "hi")] = t(quants)			# add confidence intervals to prediction data
	
	newdata$Length <- (newdata$sl * sqrt(var(set_dat$Length))) + mean(set_dat$Length)	# convert scaled length back to normal 

	# store newdata

	newdata$comp2 <- rep(set_dat$comp2[1],dim(newdata)[1])		# record comparison name	
	newdata$mod <- rep(paste(fmod), dim(newdata)[1])		# record model that was selected by AIC

	# get residuals and store them (with original data)
	set_plot <- set_dat
	set_plot$res <- residuals(get(paste(fmod)), type="deviance")

	newdata <- newdata[order(newdata$Length),]	# sort by length now for plots later

	plotstore <- rbind(plotstore,set_plot)
	newstore <- rbind(newstore,newdata)
}

#write newdat to avoid having to run loop 

#write.csv(newstore,"Model data 3 knots.csv", row.names=FALSE)
#write.csv(plotstore,"Model residuals 3 knots.csv", row.names=FALSE)

# load data instead of running above loop
newstore <- read.csv("Model data 3 knots.csv", header=TRUE)
plotstore <- read.csv("Model residuals 3 knots.csv", header=TRUE)

# need to raise both total and prop for plotting 
# this is important. Use unraised catches in model, but raised catches after for plotting.

#raising
plotstore$test <- round(plotstore$test*plotstore$samp_test)
plotstore$control <- round(plotstore$control*plotstore$samp_cont)
plotstore$cpuetest <- plotstore$test/plotstore$areatest		# get raised cpues
plotstore$cpuecont <- plotstore$control/plotstore$areacont	# get raised cpues
plotstore$cpuetotal <- plotstore$cpuetest + plotstore$cpuecont	# get raised weighting factor
plotstore$prop <- plotstore$cpuetest/plotstore$cpuetotal	# get raised prop	
plotstore$cpuetotalN <- plotstore$cpuetotal/100			# rescale cpue

# create categorical representation of weighting factors for colouring
plotstore$cpuecat <- cut(plotstore$cpuetotalN, breaks=quantile(plotstore$cpuetotalN,probs=c(0,0.25,0.5,0.75,1)), labels=1:4)

colfunc <- colorRampPalette(c("gray73","black"))
mycols <-  colfunc(length(levels(plotstore$cpuecat)))	# assign colours to each cat level


#### catch comparison rates plot

# again, looped for many comparisons

setwd("C:/Users/osp813/OneDrive - Bangor University/Documents/EMFF/Survey data/Vessel comparison 2021/Figures")
tiff("Figure 4.tiff", width=14, height=14, units="cm", res=1000, pointsize=7)	# plot storage


c.pan <- floor(sqrt(length(unique(plotstore$comp2))))	# automatically determine number of cols to split screen
r.pan <- ceiling(length(unique(plotstore$comp2))/c.pan)	# automatically determine rows
par(mfrow=c(r.pan,c.pan), mar=c(2,2,2,2), oma=c(4,4,4,1))

for(j in 1:length(unique(plotstore$comp2))){		# loop through each comparison
	set_plot <- plotstore[which(plotstore$comp2==unique(plotstore$comp2)[j]),]		# subset to individual comparison
	set_plot <- droplevels(set_plot)
	new2 <- newstore[which(newstore$comp2==unique(plotstore$comp2)[j]),]			# subset to individual comparison
	new2 <- droplevels(new2)

	plot(prop~Length, las=1, xlab="", ylab="", data=set_plot, pch=16, col=mycols[set_plot$cpuecat], main=paste(set_plot$comp2[1]), xlim=range(plotstore$Length))
	lines(prop~Length, col="magenta", lwd=2, data=new2)
	polygon(c(new2$Length,rev(new2$Length)),c(new2$lo,rev(new2$hi)),border=NA,col=adjustcolor("magenta", alpha.f=0.3))
	abline(h=0.5, lty=2, lwd=2,col="grey")
}
mtext("Scallop size (mm)", side=1, line=2, outer=TRUE)
mtext("Catch comparison rate", side=2, line=2, outer=TRUE)
leglab <- c()
for(i in 1:length(levels(plotstore$cpuecat))){
	leglab <- c(leglab, paste0(round(quantile(plotstore$cpuetotalN/100,probs=c(0,0.25,0.5,0.75,1))[i],2)," to ",round(quantile(plotstore$cpuetotalN/100,probs=c(0,0.25,0.5,0.75,1))[i+1],2)))
}
legend.fun("top", legend=leglab, pch=16, col=mycols, ncol=length(leglab), bty='n', title="Total CPUE")

dev.off()

###### residual plots 


tiff("Figure A2.tiff", width=14, height=14, units="cm", res=1000, pointsize=7)

c.pan <- floor(sqrt(length(unique(plotstore$comp2))))	# automatically determine number of cols to split screen
r.pan <- ceiling(length(unique(plotstore$comp2))/c.pan)	# automatically determine rows
par(mfrow=c(r.pan,c.pan), mar=c(2,2,2,2), oma=c(4,4,1,1))


for(j in 1:length(unique(plotstore$comp2))){
	set_plot <- plotstore[which(plotstore$comp2==unique(plotstore$comp2)[j]),]
	set_plot <- droplevels(set_plot)

	plot(res~Length, las=1, xlab="", ylab="", data=set_plot, pch=16, main=paste(set_plot$comp2[1]), xlim=range(plotstore$Length), ylim=range(plotstore$res))
	abline(h=0, lty=2, lwd=2,col="grey")
}
mtext("Scallop size (mm)", side=1, line=2, outer=TRUE)
mtext("Residuals", side=2, line=2, outer=TRUE)

dev.off()

###### END