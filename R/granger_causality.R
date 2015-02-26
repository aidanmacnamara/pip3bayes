require(MSBVAR)
require(plyr)
require(reshape2)


# DATA --------------------------------------------------------------------

# check the 2 time series for causality and the time-lag from iSH2 >> PH
data = read.csv("../share/data/single_cell_scaled.csv")

# remove messy data
toRemove = data.frame(exper=c(2,3,3,3,3,3,3,3), run=c(2,2,3,4,5,6,8,10))
for(i in 1:dim(toRemove)[1]) {
  data = data[!(data$experiment==toRemove[i,1] & data$run==toRemove[i,2]),]
}

data$Label = paste("Replicate ",data$experiment,","," Cell ",data$run, sep="")
pdf("../figures/figure_s1.pdf", height=8, width=12)
ggplot(data) + geom_line(aes(x=time, y=value_post_inhib_scaled, color=readout)) + facet_wrap(~Label, ncol=4) + theme_thesis(15) + ylab("Activation") + xlab("Time")
dev.off()

# try using the time course before perturbation
data = data[data$time<2175,]

# use the whole time-series to check causality
ish2 = data[data$readout=='iSH2', 'value_post_inhib_scaled']
ph = data[data$readout=='PH', 'value_post_inhib_scaled']
# ish2 = ts(data[(data$experiment==1 & data$run==1 & data$readout=='iSH2'),'value_post_inhib_scaled'])
# ph = ts(data[(data$experiment==1 & data$run==1 & data$readout=='PH'),'value_post_inhib_scaled'])
allData = ts(cbind(ish2,ph))
plot(allData)


# PER CELL ----------------------------------------------------------------

# question is the departure from the mean across the PH single cell traces noise
# or explained by the 'iSH2' traces?
# apply a sliding window across each trace and sample from noise from the mean

g <- function(ish2, ph, lag=1, plotFig=TRUE) {
  allData = ts(cbind(ish2,ph))
  checkLag = rep(NA,length(lag))
  for(i in 1:length(lag)) {
    checkLag[i] = granger.test(allData, p=lag[i])[2]
  }
  if (plotFig) plot(1:length(lag), checkLag, type="l")
  return(checkLag)
}

# there seems to be a peak in correlation at 7*30 seconds lag, is this shown in each single cell?

perCell = ddply(data, .(run, experiment),
      function(x) g(x[x$readout=='iSH2','value_post_inhib_scaled'],x[x$readout=='PH','value_post_inhib_scaled'], 1:5, FALSE)
)
names(perCell)[3:7] = paste("n=",1:5, sep="")
perCellPlot = melt(perCell, measure=3:7)
names(perCellPlot)[3:4] = c("n","Causality")
fig1Plot = ggplot(perCellPlot, aes(n, Causality)) + geom_boxplot() + theme_thesis()
pdf("../figures/figure_s4.pdf", height=5, width=8)
fig1Plot
dev.off()
aovOut = aov(Causality~n, data=perCellPlot)
TukeyHSD(aovOut)

res = data.frame(matrix(NA,nrow=dim(perCell)[1], ncol=4))
names(res) = c("PH", "iSH2", "p PH", "p iSH2")
for(i in 1:dim(perCell)[1]) {
  ish2_single = data[(data$experiment==perCell$experiment[i] & data$run==perCell$run[i] & data$readout=='iSH2'),'value_post_inhib_scaled']
  ph_single = data[(data$experiment==perCell$experiment[i] & data$run==perCell$run[i] & data$readout=='PH'),'value_post_inhib_scaled']
  res[i,] = granger.test(ts(cbind(ish2_single, ph_single)), 1)
}
res$comp = res$iSH2 > res$PH
wilcox.test(res[res$comp==TRUE,'p iSH2'], res[res$comp==FALSE,'p PH'])

# what does the TRUE/FALSE data look like
for(i in 1:dim(perCell)[1]) {
  if(res$comp[i]==FALSE) {
    ish2_single = data[(data$experiment==perCell$experiment[i] & data$run==perCell$run[i] & data$readout=='iSH2'),'value_post_inhib_scaled']
    ph_single = data[(data$experiment==perCell$experiment[i] & data$run==perCell$run[i] & data$readout=='PH'),'value_post_inhib_scaled']
    plot(ts(cbind(ish2_single,ph_single)))
  }
}

resPlot = cbind(melt(res[,1:2]), melt(res[,3:4])[,2])
names(resPlot) = c("Causal","Strength","P")
resPlot[resPlot$Causal=="iSH2",'Strength'] = -resPlot[resPlot$Causal=="iSH2",'Strength']
ggplot(resPlot) + geom_point(aes(x=Strength,y=-log(P,10),color=Causal)) + theme_thesis()


# RANDOMIZE ---------------------------------------------------------------

# what indices should be excluded from shuffle?
# idxFinder = ddply(data, .(run, experiment), function(x) x$value_post_inhib_scaled[140:160])
# idxFinder = idxFinder[,-c(1,2)]
# names(idxFinder) = 140:160
# idxFinder = melt(idxFinder)
# ggplot(idxFinder) + geom_point(aes(x=variable,y=value))

shuffle <- function(ish2, ph, expers=22, n=20, lag=1, win=3, lastPoint=181) {
  offSet = floor(win/2)
  # shuffleIdx = c(seq(1+offSet,140,by=win), seq(160,lastPoint-offSet,by=win))
  shuffleIdx = seq(1+offSet, lastPoint-offSet, by=win)
  idx = as.numeric(sapply(lastPoint*c(0:(expers-1)), function(x) x + shuffleIdx))
  ish2_copy = ish2
  nRes = matrix(NA,nrow=length(lag),ncol=n)
  plot(ph, type="l")
  points(ish2, type="l", col="blue")
  
  for(j in 1:n) {
    for(i in idx) {
      ish2_copy[(i-offSet):(i+offSet)] = ish2[sample((i-offSet):(i+offSet),win)]
    }
    if(j==1) points(ish2_copy, type="l",col="red")
    nRes[,j] = g(ish2_copy, ph, lag, plotFig=FALSE)
  }
  return(nRes)
}

x = granger.test(allData, 1)
resRandom = data.frame(perCell[,1:2], true=NA, random=NA)

for(i in 1:dim(perCell)[1]) {
  ish2_single = data[(data$experiment==perCell$experiment[i] & data$run==perCell$run[i] & data$readout=='iSH2'),'value_post_inhib_scaled']
  ph_single = data[(data$experiment==perCell$experiment[i] & data$run==perCell$run[i] & data$readout=='PH'),'value_post_inhib_scaled']
  xy = g(ish2_single, ph_single, 1, plotFig=FALSE)
  yy = shuffle(ish2_single, ph_single, 1, 20, lag=1, win=21, lastPoint=145)
  resRandom[i,3:4] = c(xy,mean(yy))
  # boxplot(as.numeric(yy))
  # points(xy,col=2, pch=17)
}

resRandomPlot = melt(resRandom, measure=c(3,4))
names(resRandomPlot)[3:4] = c("Group","Causality")
fig1Plot = ggplot(resRandomPlot, aes(Group, Causality)) + geom_boxplot() + theme_thesis()
pdf("../figures/figure_2b.pdf", height=5, width=8)
fig1Plot
dev.off()

# test causality
x = rnorm(100)
y = c(0,x)
x = c(x,0)
plot(x,y)
plot(x[-length(x)], y[-1])
yy = sample(x,length(x))
plot(x,yy)
granger.test(ts(cbind(x,y)),1)
granger.test(ts(cbind(x,yy)),1)
