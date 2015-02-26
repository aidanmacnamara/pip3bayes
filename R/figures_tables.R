require(reshape2)
require(plyr)
require(gridExtra)

# FIGURE 1 ----------------------------------------------------------------

fig1 = read.csv("figure_data/data_1.csv", check.names=FALSE)
fig1 = melt(fig1)
names(fig1) = c("Experiment", "Value")
pdf("../figures/figure_1.pdf", height=5, width=8)
ggplot(fig1, aes(x=Value, fill=Experiment)) + geom_histogram(binwidth=0.02) + theme_thesis(20) + ylab("Count")
dev.off()


# FIGURE 2 ----------------------------------------------------------------

# remove messy data
fig2 = read.csv("figure_data/data_2.csv", check.names=TRUE)
toRemove = data.frame(exper=c(2,3,3,3,3,3,3,3), run=c(2,2,3,4,5,6,8,10))
for(i in 1:dim(toRemove)[1]) {
  fig2 = fig2[!(fig2$Experiment==toRemove[i,1] & fig2$Run==toRemove[i,2]),]
}

fig2 = ddply(fig2, .(Experiment,Run,Rep), colwise(min, .(Model.1, User.1, User.3)))
fig2 = melt(fig2, measure=4:6)
names(fig2)[4:5] = c("Model","Likelihood")

fig2Plot = ggplot(fig2, aes(x=Rep, y=Likelihood, color=Model)) + geom_line() + geom_point() + theme_thesis(20)
pdf("../figures/figure_2.pdf", height=5, width=8)
fig2Plot
dev.off()


# AIC/BIC -----------------------------------------------------------------

# calculate AIC/BIC per run
# remove messy data
tab1 = read.csv("figure_data/data_2.csv", check.names=TRUE)
toRemove = data.frame(exper=c(2,3,3,3,3,3,3,3), run=c(2,2,3,4,5,6,8,10))
for(i in 1:dim(toRemove)[1]) {
  tab1 = tab1[!(tab1$Experiment==toRemove[i,1] & tab1$Run==toRemove[i,2]),]
}

tab1 = ddply(tab1, .(Experiment,Run,Rep), colwise(min, .(Model.1, User.1, User.3)))
aic <- function(l, k) {
  out = 2*l + 2*k
  return(out)
}
bic <- function(l, k, n) {
  out = 2*l + (k*log(n))
  return(out)
}

k = c(6,4,6) # parameters for model 1, user 1, user 3
aicRes = matrix(0, nrow=dim(tab1)[1], ncol=3)
bicRes = matrix(0, nrow=dim(tab1)[1], ncol=3)
for(i in 1:3) {
  count = i+3
  aicRes[,i] = aic(tab1[,count],k[i])
  bicRes[,i] = bic(tab1[,count],k[i], 181)
  count = count+1
}

apply(aicRes, 2, mean, na.rm=TRUE)
apply(bicRes, 2, mean, na.rm=TRUE)

write.csv(cbind(aicRes[,1],bicRes[,1],
                aicRes[,2],bicRes[,2],
                aicRes[,3],bicRes[,3]),
          file="figure_data/data_5.csv", row.names=FALSE)

# SIMULATION AND DATA -----------------------------------------------------

sim = read.csv("figure_data/data_7.csv")
n=191
sim$Time = rep(seq(from=0,by=30,length=n))
sim$Name = paste("Replicate ",sim$Experiment,","," Cell ",sim$Run, sep="")
pdf("../figures/model_b.pdf", height=8, width=12)
ggplot(sim) + geom_line(aes(x=Time, y=Sim)) + geom_line(aes(x=Time, y=Data), color="red") + facet_wrap(~Name, ncol=4) + theme_thesis(15) + ylab("Activation")
dev.off()

# for experiment 2 plot
mean_sim = ddply(sim, .(Time), summarize,
               meanData=mean(Data), sdData=sd(Data),
               meanSim=mean(Sim), sdSim=sd(Sim)
               )

pdf("../figures/figure_sim.pdf", height=5, width=8)
ggplot(mean_sim, aes(x=Time)) + geom_line(aes(y=meanSim)) + theme_thesis(25) + ylab("Activation") + geom_ribbon(aes(ymin=meanSim-sdSim, ymax=meanSim+sdSim), alpha=0.1) + geom_line(aes(y=meanData), col="red") + geom_errorbar(aes(ymin=meanData-sdData, ymax=meanData+sdData), col="red", alpha=0.2)
dev.off()

# POSTERIORS --------------------------------------------------------------

post_1 = read.csv("figure_data/data_11.csv")
post_2 = read.csv("figure_data/data_13.csv")
post_1$Group = "Exper 1"
# post_1$Label = paste("Replicate ",post_1$Experiment,","," Cell ",post_1$Run,sep="")
post_2$Group = "Exper 1/2"
postAll = rbind(post_1,post_2)

pdf("../figures/figure_4.pdf", height=8, width=12)
ggplot(post_1, aes(x=Value)) + geom_density(alpha=0.5) + theme_thesis(12) + facet_wrap(~Label, ncol=4, scales="free") + geom_vline(aes(xintercept=Best), color="red") + xlab("PTEN") + ylab("")
dev.off()

# is there a systematic shift in PTEN with the addition of the H2O2 experiment?
mean_1 = ddply(post_1, .(Experiment, Run), summarize, SD=sd(Value), Mean=mean(Value))
mean_2 = ddply(post_2, .(Experiment, Run), summarize, SD=sd(Value), Mean=mean(Value))
mean_1 = mean_1[-12,] # take out replicate 3, run 7
mean_2 = mean_2[-12,]
meanAll = cbind(mean_1[,-3], mean_2[,4])
names(meanAll)[3:4] = c("Exper 1", "Exper 1/2")
meanAll = melt(meanAll, measure=3:4)

pdf("../figures/figure_6.pdf", height=5, width=8)
ggplot(meanAll, aes(variable,value)) + geom_boxplot() + theme_thesis() + xlab("") + ylab(expression(k[f])) + theme(axis.title.y=element_text(angle=0, vjust=0.5))
dev.off()

ratio = mean_1$Mean / mean_2$Mean


# EXPLAIN GRANGER ---------------------------------------------------------

fig = read.csv("figure_data/data_6.csv")
fig$Time = seq(from=0,by=30,length.out=dim(fig)[1])
fig = melt(fig,measure=1:2)
fig[fig$variable=='PH','iSH2.random'] = NA
names(fig)[3:4] = c("Species","Activity")
pdf("../figures/figure_2a.pdf", height=8, width=8)
ggplot(fig) + geom_line(aes(x=Time, y=Activity)) + geom_line(aes(x=Time, y=iSH2.random), col="red") + facet_wrap(~Species, ncol=1) + theme_thesis(25)
dev.off()


# COVARIANCE --------------------------------------------------------------

source("ggplot_themes.R")
name_1 = "PI3K"
name_2 = "PTEN"
post_1 = read.csv("figure_data/data_11.csv")
post_1$Parameter = name_1
post_2 = read.csv("figure_data/data_10.csv")
post_2$Parameter = name_2

# experiment
corrAll = rep(NA,14)
for(i in 1:14) {
  dat = rbind(post_1[post_1$Count==i, c('Parameter','Value')],
              post_2[post_2$Count==i, c('Parameter','Value')]
  )
  samp = sample(1:10000, 1000, replace=FALSE)
  dat = dat[c(samp,samp+10000),]
  
  myScaleX <- scale_x_continuous(limits=c(0.0023,0.0029), expand=c(.00005,.00005))
  myScaleXX <- scale_x_continuous(limits=c(0.0245,0.027), expand=c(.00005,.00005))
  myScaleXY <- scale_y_continuous(limits=c(0.0245,0.027), expand=c(.00005,.00005))
  
  histTop <- ggplot(dat, aes(x=Value[Parameter==name_1])) + geom_histogram(fill="red", alpha=0.2) + theme_bw() + themeHT + ylab(" ") + myScaleX
  
  scatter <- ggplot(dat, aes(x=Value[Parameter==name_1], y=Value[Parameter==name_2],)) + geom_point(shape=1, alpha=0.2) + geom_smooth(method=lm) + theme_bw() + themeS + ylab(name_2) + xlab(name_1) + myScaleX + myScaleXY
  
  histRight <- ggplot(dat, aes(x=Value[Parameter==name_2])) + geom_histogram(fill="blue", alpha=0.2) + coord_flip() + theme_bw() + themeHR + ylab(" ") + myScaleXX
  
  # pdf(file=paste("../figures/posteriors/post_",i,".pdf",sep=""), height=5, width=9)
  pdf(file="../figures/figure_5.pdf", height=5, width=9)
  grid.arrange(histTop, histEmp, scatter, histRight, ncol=2, nrow=2, widths=c(3.5, 1.5), heights=c(1.5, 3.5))
  dev.off()
  corrAll[i] = cor.test(dat[dat$Parameter=="PI3K", 'Value'], dat[dat$Parameter=="PTEN", 'Value'])$estimate
}


# EXPERIMENT 2 DATA -------------------------------------------------------

data = read.csv("../R/figure_data/data_14.csv")
data$Label = paste("Replicate ",data$experiment,","," Cell ",data$run, sep="")
pdf("../figures/figure_s2.pdf", height=10, width=12)
ggplot(data) + geom_line(aes(x=time, y=value_post_inhib_scaled, color=readout)) + facet_wrap(~Label, ncol=4) + theme_thesis(15) + ylab("Activation") + xlab("Time")
dev.off()

