# General Setup ----
# Load Required Packages
library(here)
library(RColorBrewer)
options(scipen=9999)

#Source Utility Functions
source('src.R')

#Load Results
load(file=here('simres.RData'))

# Plot Experiment 1a ----
N = c(100,500,1000)
L = c(0.5,1,5,10)
mu.e = mu.d = mu = 0.01
cols = brewer.pal(4, 'Set2')

pdf(file = here("figures","experiment1a_homogeneity.pdf"), width = 10, height = 4)
par(mfrow=c(1,3))
for (i in 1:length(N))
{
	tmp.decode = subset(pspace.obj.decode,pop==N[i])	
	tmp.encode = subset(pspace.obj.encode,pop==N[i])	
	plot(NULL,xlim=c(0.5,4.5),ylim=c(0,1),xlab=expression(lambda),ylab='Homogeneity',axes=FALSE,main=paste0('N=',N[i]))
	for (j in 1:length(L))
	{
	boxplot(at=j-0.32,x=subset(tmp.decode,lambda==L[j])$hom.obj,col=cols[1],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j-0.12,x=subset(tmp.decode,lambda==L[j])$hom.mental,col=cols[2],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j+0.12,x=subset(tmp.encode,lambda==L[j])$hom.obj,col=cols[3],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	boxplot(at=j+0.32,x=subset(tmp.encode,lambda==L[j])$hom.mental,col=cols[4],add=T,axes=F,boxwex=0.3,outline=FALSE,whisklty=1)
	}
	axis(1,at=c(1,2,3,4),labels=L)
	axis(2)
	box()
}

legend('topright',legend=c('Decoding Error - Objects','Decoding Error - Mental Representations','Encoding Error - Objects','Encoding Error - Mental Representations'),fill=cols,bty='n')
dev.off()


pdf(file = here("figures","experiment1a_richness.pdf"), width = 10, height = 4)
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
	}
	axis(1,at=c(1,2,3,4),labels=L)
	axis(2)
	box()
	if (i==1)
	{
		legend('topleft',legend=c('Decoding Error - Objects','Decoding Error - Mental Representations','Encoding Error - Objects','Encoding Error - Mental Representations'),fill=cols,bty='n')
	}
}
dev.off()


# Plot Experiment 1b ----
nsim = 1000
n.objects = 1000
N = c(50,500,100,2000)
lambda = n.objects/N
pspace = data.frame(N=rep(N,each=nsim),lambda=rep(lambda,each=nsim))
timesteps = 5000
pdf(file = here("figures","experiment1b_richnesshomogeneity.pdf"), width = 8, height = 9)

par(mfrow=c(2,1))
xs=jitter(rep(c(1,4,7,10,2,5,8,11,14),each=nsim),factor=1.2)

# Homogeneity
ys.hom=(c(res.decode.hom,res.encode.hom,res.wf.hom))

lo.hom=quantile(res.wf.hom,0.025)
hi.hom=quantile(res.wf.hom,0.975)
col.hom=vector(length=length(ys.hom))

col.hom[ys.hom>=lo.hom&ys.hom<=hi.hom]=rgb(0,0,0,0.05)
col.hom[ys.hom<=lo.hom]=add.alpha("royalblue",0.05)
col.hom[ys.hom>=hi.hom]=add.alpha("indianred",0.05)

plot(xs,ys.hom,col=col.hom,pch=20,axes=F,xlab="",ylab="Homogeneity")
abline(h=c(lo.hom,hi.hom),col=c("royalblue","indianred"),lty=4)
axis(side=1,at=c(1.5,4.5,7.5,10.5,14),labels=c(N,n.objects))
mtext(side=1,line=2.5,"N")
axis(side=2)
axis(side=3,at=c(1.5,4.5,7.5,10.5,14),labels=c(n.objects/N,"NA\n (Wright-Fisher)"))
mtext(side=3,line=2.5,expression(lambda))
box()

text(4,0.8,"Decoding Error",srt=90,cex=0.8)
text(5,0.8,"Enccoding Error",srt=90,cex=0.8)

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

# Plot Experiment 2 ----
lambda = c(0.1,1,5,10)
N = 300
mu = 0.01

## Wright Fisher Progeny
pdf(file = here("figures","experiment2_wrightfisher_progeny.pdf"), width = 9, height = 5)

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

## Object Mediated Progeny
n.objects = 1000
N = c(50,500,1000,2000)
lambda = n.objects/N


pdf(file = here("figures","experiment2_objectmediated_progeny.pdf"), width = 9, height = 9)

par(mfrow=c(2,2))
# (decoding error)
cc <- brewer.pal(length(lambda),'Set1')
plot(wf.prog.n300.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="a")
for (i in 1:length(decode.prog.varlambda))
{
  lines(decode.prog.varlambda[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="b")
for (i in 1:length(decode.prog.fixobjects))
{
  lines(decode.prog.fixobjects[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")

# (encoding error)
plot(wf.prog.n300.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="c")
for (i in 1:length(encode.prog.varlambda))
{
  lines(encode.prog.varlambda[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="d")
for (i in 1:length(encode.prog.fixobjects))
{
  lines(encode.prog.fixobjects[[i]]$d2.objects, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")
dev.off()
