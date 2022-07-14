# General Setup ----
# Load Required Packages
library(here)
library(foreach)
library(doParallel)
library(RColorBrewer)

#Source Utility Functions
source('src.R')

# Set Number of Cores 
ncores <- 15

# Experiment 1 ----
nsim = 1000
N  = c(100,500,1000)
L  = c(0.5,1,5,10)
mu.e =mu.d = mu = 0.01
pspace.obj = pspace.wf = expand.grid(nsim=1:nsim,lambda=L,pop=N)
# pspace.wf$N.wf = pspace.wf$lambda * pspace.wf$N
# pspace.wf = unique(pspace.wf[,-c(2,3)])
timesteps <- 5000 #number of timesteps

cl <-  makeCluster(ncores)  
registerDoParallel(cl)

# run models 
# WF:
# res.wf = foreach (i = 1:nrow(pspace.wf), .combine='cbind') %dopar% {
# unlist(wf(N=pspace.wf$N.wf[i],mu=mu,timesteps = 5000,output="sumstat"))
# }
# pspace.wf$div = res.wf[1,]
# pspace.wf$k = res.wf[2,]

# Object Mediated
res.decode = foreach (i = 1:nrow(pspace.obj), .combine='rbind') %dopar% {
unlist(objTr(N=pspace.obj$pop[i],lambda=pspace.obj$lambda[i],mu.e=0,mu.d=mu.d,timesteps = 5000,output="sumstat"))}
pspace.obj.decode = cbind.data.frame(pspace.obj,res.decode)

res.encode = foreach (i = 1:nrow(pspace.obj), .combine='rbind') %dopar% {
unlist(objTr(N=pspace.obj$pop[i],lambda=pspace.obj$lambda[i],mu.e=mu.e,mu.d=0,timesteps = 5000,output="sumstat"))}
pspace.obj.encode = cbind.data.frame(pspace.obj,res.encode)

#stop cluster
stopCluster(cl)

# Plot 
cols = brewer.pal(4, 'Set2')

pdf(file = here("figures","experiment1_diversity.pdf"), width = 10, height = 4)
par(mfrow=c(1,3))
for (i in 1:length(N))
{
	tmp.decode = subset(pspace.obj.decode,pop==N[i])	
	tmp.encode = subset(pspace.obj.encode,pop==N[i])	
	plot(NULL,xlim=c(0.5,4.5),ylim=c(0,1),xlab=expression(lambda),ylab='Diversity',axes=FALSE,main=paste0('N=',N[i]))
	for (j in 1:length(L))
	{
	boxplot(at=j-0.32,x=subset(tmp.decode,lambda==L[j])$div.obj,col=cols[1],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j-0.12,x=subset(tmp.decode,lambda==L[j])$div.mental,col=cols[2],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j+0.12,x=subset(tmp.encode,lambda==L[j])$div.obj,col=cols[3],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j+0.32,x=subset(tmp.encode,lambda==L[j])$div.mental,col=cols[4],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
# 	points(x=jitter(rep(j-0.30,nsim),factor=jf),y=subset(tmp.decode,lambda==L[j])$div.obj,pch=20,col=cols[1])
# 	points(x=jitter(rep(j-0.12,nsim),factor=jf),y=subset(tmp.decode,lambda==L[j])$div.mental,pch=20,col=cols[2])
# 	points(x=jitter(rep(j+0.12,nsim),factor=jf),y=subset(tmp.encode,lambda==L[j])$div.obj,pch=20,col=cols[3])
# 	points(x=jitter(rep(j+0.30,nsim),factor=jf),y=subset(tmp.encode,lambda==L[j])$div.mental,pch=20,col=cols[4])
   	
	}
	axis(1,at=c(1,2,3,4),labels=L)
	axis(2)
	box()
}

legend('bottomright',legend=c('Decoding Error - Objects','Decoding Error - Mental Templates','Encoding Error - Objects','Encoding Error - Mental Templates'),fill=cols,bty='n')
dev.off()


pdf(file = here("figures","experiment1_richness.pdf"), width = 10, height = 4)
par(mfrow=c(1,3))
maxK = max(c(pspace.obj.decode$k.obj,pspace.obj.decode$k.mental,pspace.obj.encode$k.obj,pspace.obj.encode$k.mental))
for (i in 1:length(N))
{
	tmp.decode = subset(pspace.obj.decode,pop==N[i])	
	tmp.encode = subset(pspace.obj.encode,pop==N[i])	
	plot(NULL,xlim=c(0.5,4.5),ylim=c(0,maxK),xlab=expression(lambda),ylab='Richness',axes=FALSE,main=paste0('N=',N[i]))
	for (j in 1:length(L))
	{
	boxplot(at=j-0.32,x=subset(tmp.decode,lambda==L[j])$k.obj,col=cols[1],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j-0.12,x=subset(tmp.decode,lambda==L[j])$k.mental,col=cols[2],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j+0.12,x=subset(tmp.encode,lambda==L[j])$k.obj,col=cols[3],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j+0.32,x=subset(tmp.encode,lambda==L[j])$k.mental,col=cols[4],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
# 	points(x=jitter(rep(j-0.30,nsim),factor=jf),y=subset(tmp.decode,lambda==L[j])$k.obj,pch=20,col=cols[1])
# 	points(x=jitter(rep(j-0.12,nsim),factor=jf),y=subset(tmp.decode,lambda==L[j])$k.mental,pch=20,col=cols[2])
# 	points(x=jitter(rep(j+0.12,nsim),factor=jf),y=subset(tmp.encode,lambda==L[j])$k.obj,pch=20,col=cols[3])
# 	points(x=jitter(rep(j+0.30,nsim),factor=jf),y=subset(tmp.encode,lambda==L[j])$k.mental,pch=20,col=cols[4])
	}
	axis(1,at=c(1,2,3,4),labels=L)
	axis(2)
	box()
	if (i==1)
	{
		legend('topleft',legend=c('Decoding Error - Objects','Decoding Error - Mental Templates','Encoding Error - Objects','Encoding Error - Mental Templates'),fill=cols,bty='n')
	}
}
dev.off()


# Experiment 2 ----

# parameter settings:
nsim = 1000 #number of repetitions for each parameter combination
n.objects = 1000 #target number of objects produced at each timestep 
N <- c(50,500,1000,2000) #number of individuals (producers)
lambda <- n.objects/N # production rate
# create parameter spaces
pspace = data.frame(N=rep(N,each=nsim), lambda=rep(lambda,each=nsim)) 
mu.e = mu.d = mu = 0.01 #error rates
timesteps <- 5000 #number of timesteps


# Setup parallel processing using n-1 threads, with n being maximum possible for the specific machine
cl <-  makeCluster(ncores)  
registerDoParallel(cl)

# run models
res.wf = foreach (i = 1:nsim, .combine='cbind') %dopar% {
unlist(wf(N=n.objects,mu=mu,timesteps = 5000,output="sumstat"))
}
res.wf.div = res.wf[1,]
res.wf.k = res.wf[2,]

res.decode = foreach (i = 1:nrow(pspace), .combine='cbind') %dopar% {
unlist(objTr(N=pspace$N[i],lambda=pspace$lambda[i],mu.e=0,mu.d=mu.d,timesteps = 5000,output="sumstat"))}
res.decode.div=matrix(res.decode[1,],ncol=length(N),nrow=nsim)
res.decode.k=matrix(res.decode[2,],ncol=length(N),nrow=nsim)

res.encode = foreach (i = 1:nrow(pspace), .combine='cbind') %dopar% {
unlist(objTr(N=pspace$N[i],lambda=pspace$lambda[i],mu.e=mu.e,mu.d=0,timesteps = 5000,output="sumstat"))}
res.encode.div=matrix(res.encode[1,],ncol=length(N),nrow=nsim)
res.encode.k=matrix(res.encode[2,],ncol=length(N),nrow=nsim)

#stop cluster
stopCluster(cl)



# Plot Results (pdf)
pdf(file = here("figures","experiment2_richnessdiversity.pdf"), width = 8, height = 9)

par(mfrow=c(2,1))
xs=jitter(rep(c(1,4,7,10,2,5,8,11,14),each=nsim),factor=1.2)

# Diversity
ys.div=(c(res.decode.div,res.encode.div,res.wf.div))

lo.div=quantile(res.wf.div,0.025)
hi.div=quantile(res.wf.div,0.975)
col.div=vector(length=length(ys.div))

col.div[ys.div>=lo.div&ys.div<=hi.div]=rgb(0,0,0,0.05)
col.div[ys.div<=lo.div]=add.alpha("royalblue",0.05)
col.div[ys.div>=hi.div]=add.alpha("indianred",0.05)

plot(xs,ys.div,col=col.div,pch=20,axes=F,xlab="",ylab="Diversity")
abline(h=c(lo.div,hi.div),col=c("royalblue","indianred"),lty=4)
axis(side=1,at=c(1.5,4.5,7.5,10.5,14),labels=c(N,n.objects))
mtext(side=1,line=2.5,"N")
axis(side=2)
axis(side=3,at=c(1.5,4.5,7.5,10.5,14),labels=c(n.objects/N,"NA\n (Wright-Fisher)"))
mtext(side=3,line=2.5,expression(lambda))
box()

text(4,0.2,"Decoding Error",srt=90,cex=0.8)
text(5,0.2,"Enccoding Error",srt=90,cex=0.8)

# Richness
ys.k=(c(res.decode.k,res.encode.k,res.wf.k))

lo.k=quantile(res.wf.k,0.025)
hi.k=quantile(res.wf.k,0.975)
col.k=vector(length=length(ys.k))

col.k[ys.k>=lo.k&ys.k<=hi.k]=rgb(0,0,0,0.05)
col.k[ys.k<=lo.k]=add.alpha("royalblue",0.05)
col.k[ys.k>=hi.k]=add.alpha("indianred",0.05)

plot(xs,ys.k,col=col.k,pch=20,axes=F,xlab="",ylab="Richness (k)")
abline(h=c(lo.k,hi.k),col=c("royalblue","indianred"),lty=4)
axis(side=1,at=c(1.5,4.5,7.5,10.5,14),labels=c(N,n.objects))
mtext(side=1,line=2.5,"N")
axis(side=2)
axis(side=3,at=c(1.5,4.5,7.5,10.5,14),labels=c(n.objects/N,"NA\n (Wright-Fisher)"))
mtext(side=3,line=2.5,expression(lambda))
box()

text(1,44,"Decoding Error",srt=90,cex=0.8)
text(2,44,"Enccoding Error",srt=90,cex=0.8)
dev.off()

# Experiment 3 ----

# Expectation Under Neutrality (Figure)

# changing population size
wf.prog.n300.m01 = wf(N=300,mu=0.01,warmup=10000,timesteps=20000,output="progeny")
wf.prog.n1000.m01 = wf(N=1000,mu=0.01,warmup=10000,timesteps=20000,output="progeny")
wf.prog.n3000.m01 = wf(N=3000,mu=0.01,warmup=10000,timesteps=20000,output="progeny")

# changing innovation rate
wf.prog.n1000.m005 = wf(N=300,mu=0.005,warmup=10000,timesteps=20000,output="progeny")
wf.prog.n3000.m001 = wf(N=300,mu=0.001,warmup=10000,timesteps=20000,output="progeny")                                             

#
pdf(file = here("figures","experiment3_wrightfisher_progeny.pdf"), width = 9, height = 5)

par(mfrow=c(1,2))
plot(wf.prog.n300.m01$d2,pch=20,log="xy",ylab="Probability P(k) of number variants >= k",xlab="k",col="darkgrey",type="b")
lines(wf.prog.n1000.m01$d2,pch=2,col="indianred",type="b")
lines(wf.prog.n3000.m01$d2,pch=3,col="royalblue",type="b")
legend("bottomleft",legend=c(expression(paste("N=300; ",mu,"=0.01")),expression(paste("N=1000; ",mu,"=0.01")),expression(paste("N=3000; ",mu,"=0.01"))),pch=c(20,2,3),col=c("darkgrey","indianred","royalblue"),bty="n")


plot(wf.prog.n300.m01$d2,pch=20,log="xy",ylab="Probability P(k) of number variants >= k",xlab="k",col="grey",type="b")
lines(wf.prog.n1000.m005$d2,pch=2,col="indianred",type="b")
lines(wf.prog.n3000.m001$d2,pch=3,col="royalblue",type="b")
legend("bottomleft",legend=c(expression(paste("N=300; ",mu,"=0.01")),expression(paste("N=300; ",mu,"=0.005")),expression(paste("N=300; ",mu,"=0.001"))),pch=c(20,2,3),col=c("darkgrey","indianred","royalblue"),bty="n")                                             
dev.off()



#Progeny Distribution under object-mediated transmission

# WF for reference
mu <- 0.01
N <- 300
warmup <- 10000
timesteps <- warmup * 2
lambda <- c(0.1,1,5,10)
                 
decode.prog.varlambda <- vector("list",length=length(lambda))
encode.prog.varlambda <- vector("list",length=length(lambda))
legItems.progeny1=vector(length=length(lambda))

for (i in 1:length(lambda))
{
  decode.prog.varlambda[[i]] <-objTr(N=N,mu.e=0,mu.d=mu,lambda=lambda[i],warmup=warmup,timesteps=timesteps,output="progeny")
  encode.prog.varlambda[[i]] <-objTr(N=N,mu.e=mu,mu.d=0,lambda=lambda[i],warmup=warmup,timesteps=timesteps,output="progeny")  
  legItems.progeny1[i] = as.expression(bquote(paste(lambda,"=",.(lambda[i]),"; N=",.(N),sep="")))
}

legItems.progeny1=c(as.expression("Wright-Fisher"),legItems.progeny1)


# actual object mediated transmission model

n.objects <- 1000
N <- c(50,500,1000,2000)
lambda <- n.objects/N                                                                                                             

encode.prog.fixobjects <- vector("list",length=length(N))
decode.prog.fixobjects <- vector("list",length=length(N))
legItems.progeny2=vector(length=length(N))


for (i in 1:length(N))
{
  decode.prog.fixobjects[[i]] <- objTr(N=N[i],mu.e=0,mu.d=mu,lambda=lambda[i],warmup=warmup,timesteps=timesteps,output="progeny")
  encode.prog.fixobjects[[i]] <- objTr(N=N[i],mu.e=mu,mu.d=0,lambda=lambda[i],warmup=warmup,timesteps=timesteps,output="progeny")  
  legItems.progeny2[i] = as.expression(bquote(paste(lambda,"=",.(lambda[i]),"; N=",.(N[i]),sep="")))
}
legItems.progeny2=c(as.expression("Wright-Fisher"),legItems.progeny2)

# Plot Results (pdf)
pdf(file = here("figures","figure5_objectmediated_progeny.pdf"), width = 9, height = 9)

par(mfrow=c(2,2))
# (decoding error)
cc <- brewer.pal(length(lambda),'Set1')
plot(wf.prog.n300.m01$d2.objects,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="a")
for (i in 1:length(decode.prog.varlambda))
{
  lines(decode.prog.varlambda[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2.objects,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2.objects,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="b")
for (i in 1:length(decode.prog.fixobjects))
{
  lines(decode.prog.fixobjects[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2.objects,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")

# (encoding error)
plot(wf.prog.n300.m01$d2.objects,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="c")
for (i in 1:length(encode.prog.varlambda))
{
  lines(encode.prog.varlambda[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2.objects,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2.objects,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="d")
for (i in 1:length(encode.prog.fixobjects))
{
  lines(encode.prog.fixobjects[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2.objects,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")
dev.off()
