# Adam Delargy 

# Vessel comparison >= MLS catch analysis 

# see the 'Scallop catch analysis' script for notes

# 29/10/2021


legend.fun<-function(...){
  opar<-par(fig=c(0,1,0,1),oma=rep(0,4),
            mar=rep(0,4), new=TRUE)
  on.exit(par(opar))
  plot(0,0,type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


library(glmmTMB)		# for mixed models 
library(bbmle)			# for AICtab

###### load data 

setwd("C:/Users/osp813/OneDrive - Bangor University/Documents/EMFF/Survey data/Vessel comparison 2021")

dat <- read.csv("Comparison data for models.csv", header=TRUE)
# note comparisons where no scallops were caught by either 
# vessel have already been eliminated

# raise catches 

dat$test <- round(dat$test*dat$samp_test)
dat$control <- round(dat$control*dat$samp_cont)

table(dat$test,dat$comp2,sum)


##### get the data for comparisons 

# subset to groups contains only those over the MLS

dat <- dat[which(dat$Length>=110),]
dat <- droplevels(dat)


storedat <- c()

for(i in 1:length(unique(dat$comp2))){
	set <- dat[which(dat$comp2==unique(dat$comp2)[i]),]
	set <- droplevels(set)
	
	dum1 <- data.frame("Station"=names(tapply(set$test, set$Station, sum)))
	dum1$test <- tapply(set$test, set$Station, sum)
	dum1$areatest <- tapply(set$areatest, set$Station, mean)
	dum1$cont <- tapply(set$control, set$Station, sum)
	dum1$areacont <- tapply(set$areacont, set$Station, mean)
	dum1$cpuetest <- (dum1$test/dum1$areatest)/100
	dum1$cpuecont <- (dum1$cont/dum1$areacont)/100
	dum1$cpuecomb <- dum1$cpuetest + dum1$cpuecont
	dum1$depth <- tapply(set$depth, set$Station, mean)
	dum1$comp2 <- rep(set$comp2[1], dim(dum1)[1])
	dum1$prop <- dum1$cpuetest/(dum1$cpuetest+dum1$cpuecont)

	storedat <- rbind(storedat,dum1)
}


# get cpuecomb as categorical 

storedat$cpuecat <- cut(storedat$cpuecomb, breaks=quantile(storedat$cpuecomb,probs=c(0,0.25,0.5,0.75,1)), labels=1:4)

colfunc <- colorRampPalette(c("gray73","black"))
mycols <-  colfunc(length(levels(storedat$cpuecat)))


### explore potential explanatory variables 

# 1. depth 

c.pan <- floor(sqrt(length(unique(storedat$comp2))))	# automatically determine number of cols to split screen
r.pan <- ceiling(length(unique(storedat$comp2))/c.pan)	# automatically determine rows
par(mfrow=c(r.pan,c.pan), oma=c(1,1,3,1), mar=c(4,4,2,1))

for(i in 1:length(unique(storedat$comp2))){
	set <- storedat[which(storedat$comp2==unique(storedat$comp2)[i]),]
	set <- droplevels(set)

	plot(set$prop ~ set$depth,las=1, pch=16, col=mycols[set$cpuecat], main=set$comp2[1], xlab="Depth (m)", ylab="", ylim=c(0,1), xlim=range(storedat$depth))
	abline(h=0.5, col="grey", lwd=2, lty=2)
}



leglab <- c()
for(i in 1:length(levels(storedat$cpuecat))){
	leglab <- c(leglab, paste0(round(quantile(storedat$cpuecomb/100,probs=c(0,0.25,0.5,0.75,1))[i],1)," to ",round(quantile(storedat$cpuecomb/100,probs=c(0,0.25,0.5,0.75,1))[i+1],1)))
}

legend.fun("top", legend=leglab, pch=16, col=mycols, ncol=length(leglab), bty='n', title="Total CPUE")

# so potentially trends in FV1 vs FV2 and FV2 vs FV3

###



#write loop to fit models 

modres <- matrix(NA, nrow=length(unique(storedat$comp2)), ncol=5)

setwd("C:/Users/osp813/OneDrive - Bangor University/Documents/EMFF/Survey data/Vessel comparison 2021/Figures")
tiff("Figure 3.tiff", width=14, height=14, units="cm", res=1000, pointsize=7)

c.pan <- floor(sqrt(length(unique(storedat$comp2))))	# automatically determine number of cols to split screen
r.pan <- ceiling(length(unique(storedat$comp2))/c.pan)	# automatically determine rows
par(mfrow=c(r.pan,c.pan), oma=c(1,1,3,1), mar=c(4,4,2,1))

stres <- c()

for(i in 1:length(unique(storedat$comp2))){
	set <- storedat[which(storedat$comp2==unique(storedat$comp2)[i]),]
	set <- droplevels(set)

	mod1 <- glmmTMB(prop~1 + (1|Station),weights=cpuecomb, data=set, family='binomial')
	mod2 <- glmmTMB(prop~depth + (1|Station),weights=cpuecomb, data=set, family='binomial')

	# compare models 
	fmod <- attr(AICtab(mod1,mod2),"row.names")[1]
	# save the model with the lowest AIC 

	# however, if model difference is 2 AIC or less then automatically
	# accept the simpler model 

	if(AICtab(mod1,mod2)[1]$dAIC[2] <= 2){fmod <- "mod1"}

	# record preferred model and signifance info 
	modres[i,1] <- set$comp2[1]
	modres[i,2] <- fmod

	# get residuals 
	set_res <- set
	set_res$res <- residuals(get(paste(fmod)), type="pearson")

	sigs <- rep(NA,2)

	sigs[1:length(summary(get(paste(fmod)))$coef$cond[,4])] <- summary(get(paste(fmod)))$coef$cond[,4]
	modres[i,3:4]  <- sigs

	newdat <- expand.grid(1,"Station"=NA,"cpuecomb"=1, "depth"=seq(min(storedat$depth),max(storedat$depth),length.out=25))	
	
	newdat$pred <- predict(get(fmod), newdata=newdat, type="response")

	# get CIs
	
	if(fmod=="mod1"){
		X <- model.matrix(~ 1 , data = newdat)
	}
	if(fmod=="mod2"){
		X <- model.matrix(~ depth , data = newdat)
	}
	#Extract parameters and parameter covariance matrix	
	betas    <- fixef(get(fmod))
	
	#Calculate the fitted values in the predictor scale
	newdat$eta <- as.matrix(X) %*% unlist(betas)
	newdat$Pi  <- exp(newdat$eta) / (1 + exp(newdat$eta))

	#Calculate the SEs on the scale of the predictor function
	newdat$se    <- sqrt(diag(as.matrix(X) %*% matrix(unlist(vcov(get(fmod))),nrow=sqrt(length(unlist(vcov(get(fmod))))),ncol=sqrt(length(unlist(vcov(get(fmod)))))) %*% t(as.matrix(X))))
	newdat$SeUp  <- exp(newdat$eta + 1.96 *newdat$se) / (1 + exp(newdat$eta  + 1.96 *newdat$se))	
	newdat$SeLo  <- exp(newdat$eta - 1.96 *newdat$se) / (1 + exp(newdat$eta  - 1.96 *newdat$se))

	
	# store average 
	modres[i,5] <- mean(newdat$pred)

	if(fmod=="mod1"){
		plot(set$prop,las=1, pch=16, col=mycols[set$cpuecat], main=set$comp2[1], xlab="Haul", ylab="", ylim=c(0,1), xaxt='n')
		abline(h=mean(newdat$pred), col="magenta", lwd=2)
		abline(h=0.5, col="grey", lwd=2, lty=2)
		polygon(c(seq(0,(dim(set)[1])+1,length.out=dim(newdat)[1]),rev(seq(0,dim(set)[1]+1,length.out=dim(newdat)[1]))),c(newdat$SeLo,rev(newdat$SeUp)),border=NA,col=adjustcolor("magenta", alpha.f=0.3))

	}
	if(fmod=="mod2"){
		plot(set$prop ~ set$depth,las=1, pch=16, col=mycols[set$cpuecat],main=set$comp2[1], xlab="Depth (m)", ylab="", ylim=c(0,1), xlim=range(storedat$depth))
		lines(newdat$pred~ newdat$depth, col="magenta", lwd=2)
		abline(h=0.5, col="grey", lwd=2, lty=2)
		polygon(c(newdat$depth,rev(newdat$depth)),c(newdat$SeLo,rev(newdat$SeUp)),border=NA,col=adjustcolor("magenta", alpha.f=0.3))
	}

	stres <- rbind(stres,set_res)
	
	
}
mtext("Catch comparison rate", side=2, line=-1, outer=TRUE)
legend.fun("top", legend=leglab, pch=16, col=mycols, ncol=length(leglab), bty='n', title=expression(Total ~ CPUE ~ (nb ~ per ~ 100~m^2)))

dev.off()


# ignore warnings 

# write results 

modat <- as.data.frame(modres)
colnames(modat) <- c("comp2","model","int_pval","depth_pval","meanest")

# write 
#write.csv(modat, "MLS catch res.csv", row.names=FALSE)

################
