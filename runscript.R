# General Setup ----
# Load Required Packages
library(here)
library(foreach)
library(doParallel)
library(RColorBrewer)

#Source Utility Functions
source('src.R')

# Set Number of Cores 
ncores <- 6


# Experiment 1 ----

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
pdf(file = here("figures","figure2_richnessdiversity.pdf"), width = 8, height = 9)

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


# Plot Results (tiff)
tiff(file = here("figures","figure2_richnessdiversity.tiff"), width = 8, height = 9, units = "in", res = 300)

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




# Experiment 2 ----

# parameter setting
nsim = 100 #number of repetitions for each parameter combination
n.objects = 1000 #target number of objects produced at each timestep 
N <- c(50,500,1000,2000) #number of individuals (producers)
lambda <- n.objects/N # production rate
# create parameter spaces
pspace = data.frame(N=rep(N,each=nsim), lambda=rep(lambda,each=nsim)) 
mu.e = mu.d = mu = 0.01 #error rates
timesteps <- 5500 #number of timesteps
warmup <- 5000
top <- 10

#Setup parallel processing using n-1 threads, with n being maximum possible for the specific machine
cl <-  makeCluster(ncores)  
registerDoParallel(cl)

#Run Model
wf.tr = foreach (i=1:nsim) %dopar% {
tmp=wf(N = n.objects,mu=mu,timesteps=5501,warmup=5000,output="turnover",verbose=FALSE,top=top)
list(tmp$z.frame[,2],tmp$b)
}

encode.tr = foreach(i=1:nrow(pspace)) %dopar% {
tmp=objTr(N = pspace$N[i],lambda=pspace$lambda[i],mu.e=mu.e,mu.d=0,timesteps=5501,warmup=5000,output="turnover",verbose=FALSE,top=top)
list(tmp$z.frame[,2],tmp$bN,tmp$bNp)  
}

decode.tr = foreach(i=1:nrow(pspace)) %dopar% {
tmp=objTr(N = pspace$N[i],lambda=pspace$lambda[i],mu.e=0,mu.d=mu.d,timesteps=5501,warmup=5000,output="turnover",verbose=FALSE,top=top)
list(tmp$z.frame[,2],tmp$bN,tmp$bNp)  
}

#stop cluster
stopCluster(cl)


#Post Processing
wf.rates=matrix(unlist(lapply(wf.tr,'[[',1)),nrow=10)
wf.b=unlist(lapply(wf.tr,'[[',2))

encode.rates.mat=matrix(unlist(lapply(encode.tr,'[[',1)),nrow=top)
encode.rates = vector("list",length=length(N))
encode.b=matrix(unlist(lapply(encode.tr,'[[',2)),nrow=nsim)

decode.rates.mat=matrix(unlist(lapply(decode.tr,'[[',1)),nrow=top)
decode.rates = vector("list",length=length(N))
decode.b=matrix(unlist(lapply(decode.tr,'[[',2)),nrow=nsim)

legItems.turnover=vector(length=length(lambda))


for (i in 1:length(N))
{
  encode.rates[[i]] = encode.rates.mat[,(i*nsim-(nsim-1)):(i*nsim)]
  decode.rates[[i]] = decode.rates.mat[,(i*nsim-(nsim-1)):(i*nsim)]
  legItems.turnover[i] = as.expression(bquote(paste(lambda,"=",.(lambda[i]),"; N=",.(N[i]),sep="")))
}

legItems.turnover = c(legItems.turnover,paste0("Wright-Fisher (N=",n.objects,")"))

                                                                                                                                                                                                                                                               
# Compute Means and 95% percentiles
encode.b.hi=apply(encode.b,2,quantile,0.975)
encode.b.mean=apply(encode.b,2,mean)
encode.b.lo=apply(encode.b,2,quantile,0.025)
decode.b.hi=apply(decode.b,2,quantile,0.975)
decode.b.mean=apply(decode.b,2,mean)
decode.b.lo=apply(decode.b,2,quantile,0.025)
wf.b.hi=quantile(wf.b,0.975)
wf.b.mean=mean(wf.b)
wf.b.lo=quantile(wf.b,0.025)

# Plot Results
## Estimates of b (pdf)
pdf(file = here("figures","figure3_turnover_b_estimate.pdf"), width = 7, height = 5)

plot(0,0,type="n",ylim=range(c(wf.b,encode.b,decode.b)),xlim=c(0.5,14.5),axes=F,xlab="",ylab="")
abline(h=0.86,lty=4)
axis(side=1,at=c(1.5,4.5,7.5,10.5,14),labels=c(N,n.objects))
mtext(side=1,line=2.5,"N")
axis(side=2)
axis(side=3,at=c(1.5,4.5,7.5,10.5,14),labels=c(n.objects/N,"NA\n (Wright-Fisher)"))
mtext(side=3,line=2.5,expression(lambda))
box()     
arrows(x0=c(1,4,7,10),x1=c(1,4,7,10),y0=encode.b.lo,y1=encode.b.hi,lty=1,lwd=2,length=0,col="indianred")
points(c(1,4,7,10),encode.b.mean,pch=20,col="indianred",cex=1.5)
arrows(x0=c(2,5,8,11),x1=c(2,5,8,11),y0=decode.b.lo,y1=decode.b.hi,lty=1,lwd=2,length=0,col="royalblue")
points(c(2,5,8,11),decode.b.mean,pch=20,col="royalblue",cex=1.5)
arrows(x0=14,x1=14,y0=wf.b.lo,y1=wf.b.hi,lwd=2,length=0)
points(14,wf.b.mean,pch=20,cex=1.5)
legend("topright",legend=c("Object-mediated with encoding error","Object-mediated with decoding error","Wright-Fisher"),lty=c(1),pch=20,col=c("indianred","royalblue","black"),cex=0.7)
dev.off()


## Estimates of b (tiff)
tiff(file = here("figures","figure3_turnover_b_estimate.tiff"), width = 7, height = 5,units = "in", res = 300)

plot(0,0,type="n",ylim=range(c(wf.b,encode.b,decode.b)),xlim=c(0.5,14.5),axes=F,xlab="",ylab="")
abline(h=0.86,lty=4)
axis(side=1,at=c(1.5,4.5,7.5,10.5,14),labels=c(N,n.objects))
mtext(side=1,line=2.5,"N")
axis(side=2)
axis(side=3,at=c(1.5,4.5,7.5,10.5,14),labels=c(n.objects/N,"NA\n (Wright-Fisher)"))
mtext(side=3,line=2.5,expression(lambda))
box()     
arrows(x0=c(1,4,7,10),x1=c(1,4,7,10),y0=encode.b.lo,y1=encode.b.hi,lty=1,lwd=2,length=0,col="indianred")
points(c(1,4,7,10),encode.b.mean,pch=20,col="indianred",cex=1.5)
arrows(x0=c(2,5,8,11),x1=c(2,5,8,11),y0=decode.b.lo,y1=decode.b.hi,lty=1,lwd=2,length=0,col="royalblue")
points(c(2,5,8,11),decode.b.mean,pch=20,col="royalblue",cex=1.5)
arrows(x0=14,x1=14,y0=wf.b.lo,y1=wf.b.hi,lwd=2,length=0)
points(14,wf.b.mean,pch=20,cex=1.5)
legend("topright",legend=c("Object-mediated with encoding error","Object-mediated with decoding error","Wright-Fisher"),lty=c(1),pch=20,col=c("indianred","royalblue","black"),cex=0.7)
dev.off()





## Turnover Profile (pdf)
pdf(file = here("figures","figure4_turnover_profile.pdf"), width = 9, height = 5)
par(mfrow=c(1,2))
plot(1:10,apply(wf.rates,1,mean,na.rm=TRUE),type="b",xlab="Top",ylab="Turnover Rate",pch=20,ylim=range(c(wf.rates,unlist(encode.rates)),na.rm=TRUE),main="Encoding Error")
arrows(x0=1:10,x1=1:10,y0=apply(wf.rates,1,quantile,0.025,na.rm=TRUE),y1=apply(wf.rates,1,quantile,0.975,na.rm=TRUE),length=0)
cc <- brewer.pal(length(N),"Set1") 
for (i in 1:length(N))
{
  lines(1:10,apply(encode.rates[[i]],1,mean,na.rm=TRUE),type="b",pch=20,col=cc[i])  
  arrows(x0=1:10,x1=1:10,y0=apply(encode.rates[[i]],1,quantile,0.025,na.rm=TRUE),y1=apply(encode.rates[[i]],1,quantile,0.975,na.rm=TRUE),length=0,col=cc[i])
}
legend("topleft",legend=legItems.turnover,col=c(cc,"black"),pch=20,lty=1,cex=0.8,bty="n")

plot(1:10,apply(wf.rates,1,mean,na.rm=TRUE),type="b",xlab="Top",ylab="Turnover Rate",pch=20,ylim=range(c(wf.rates,unlist(decode.rates)),na.rm=TRUE),main="Decoding Error")
arrows(x0=1:10,x1=1:10,y0=apply(wf.rates,1,quantile,0.025,na.rm=TRUE),y1=apply(wf.rates,1,quantile,0.975,na.rm=TRUE),length=0)

cc <- brewer.pal(length(N),"Set1") 
for (i in 1:length(N))
{
  lines(1:10,apply(decode.rates[[i]],1,mean,na.rm=TRUE),type="b",pch=20,col=cc[i])  
  arrows(x0=1:10,x1=1:10,y0=apply(decode.rates[[i]],1,quantile,0.025,na.rm=TRUE),y1=apply(decode.rates[[i]],1,quantile,0.975,na.rm=TRUE),length=0,col=cc[i])
}
dev.off()



## Turnover Profile (tiff)
tiff(file = here("figure","figure4_turnover_profile.tiff"), width = 9, height = 5,units = "in", res = 300)
par(mfrow=c(1,2))
plot(1:10,apply(wf.rates,1,mean,na.rm=TRUE),type="b",xlab="Top",ylab="Turnover Rate",pch=20,ylim=range(c(wf.rates,unlist(encode.rates)),na.rm=TRUE),main="Encoding Error")
arrows(x0=1:10,x1=1:10,y0=apply(wf.rates,1,quantile,0.025,na.rm=TRUE),y1=apply(wf.rates,1,quantile,0.975,na.rm=TRUE),length=0)
cc <- brewer.pal(length(N),"Set1") 
for (i in 1:length(N))
{
  lines(1:10,apply(encode.rates[[i]],1,mean,na.rm=TRUE),type="b",pch=20,col=cc[i])  
  arrows(x0=1:10,x1=1:10,y0=apply(encode.rates[[i]],1,quantile,0.025,na.rm=TRUE),y1=apply(encode.rates[[i]],1,quantile,0.975,na.rm=TRUE),length=0,col=cc[i])
}
legend("topleft",legend=legItems.turnover,col=c(cc,"black"),pch=20,lty=1,cex=0.8,bty="n")

plot(1:10,apply(wf.rates,1,mean,na.rm=TRUE),type="b",xlab="Top",ylab="Turnover Rate",pch=20,ylim=range(c(wf.rates,unlist(decode.rates)),na.rm=TRUE),main="Decoding Error")
arrows(x0=1:10,x1=1:10,y0=apply(wf.rates,1,quantile,0.025,na.rm=TRUE),y1=apply(wf.rates,1,quantile,0.975,na.rm=TRUE),length=0)

cc <- brewer.pal(length(N),"Set1") 
for (i in 1:length(N))
{
  lines(1:10,apply(decode.rates[[i]],1,mean,na.rm=TRUE),type="b",pch=20,col=cc[i])  
  arrows(x0=1:10,x1=1:10,y0=apply(decode.rates[[i]],1,quantile,0.025,na.rm=TRUE),y1=apply(decode.rates[[i]],1,quantile,0.975,na.rm=TRUE),length=0,col=cc[i])
}
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
pdf(file = here("figures","figure1_wrightfisher_progeny.pdf"), width = 9, height = 5)

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


# WF (tiff)
tiff(file = here("figures","figure1_wrightfisher_progeny.tiff"), width = 9, height = 5, units = "in", res = 300)

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
timesteps <- warmup + 10000
lambda <- c(0.1,1,5,10)
                 
decode.prog.varlambda <- vector("list",length=length(lambda))
encode.prog.varlambda <- vector("list",length=length(lambda))
legItems.progeny1=vector(length=length(lambda))

for (i in 1:length(lambda))
{
  decode.prog.varlambda[[i]] <-objTr(N=N,mu.e=0,mu.d=mu,lambda=lambda[i],warmup=10000,timesteps=20000,output="progeny")
  encode.prog.varlambda[[i]] <-objTr(N=N,mu.e=mu,mu.d=0,lambda=lambda[i],warmup=10000,timesteps=20000,output="progeny")  
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
  decode.prog.fixobjects[[i]] <- objTr(N=N[i],mu.e=0,mu.d=mu,lambda=lambda[i],warmup=10000,timesteps=20000,output="progeny")
  encode.prog.fixobjects[[i]] <- objTr(N=N[i],mu.e=mu,mu.d=0,lambda=lambda[i],warmup=10000,timesteps=20000,output="progeny")  
  legItems.progeny2[i] = as.expression(bquote(paste(lambda,"=",.(lambda[i]),"; N=",.(N[i]),sep="")))
}
legItems.progeny2=c(as.expression("Wright-Fisher"),legItems.progeny2)

# Plot Results (pdf)
pdf(file = here("figures","figure5_objectmediated_progeny.pdf"), width = 9, height = 9)

par(mfrow=c(2,2))
# (decoding error)
cc <- brewer.pal(length(lambda),'Set1')
plot(wf.prog.n300.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="a")
for (i in 1:length(decode.prog.varlambda))
{
  lines(decode.prog.varlambda[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="b")
for (i in 1:length(decode.prog.fixobjects))
{
  lines(decode.prog.fixobjects[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")

# (encoding error)
plot(wf.prog.n300.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="c")
for (i in 1:length(encode.prog.varlambda))
{
  lines(encode.prog.varlambda[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="d")
for (i in 1:length(encode.prog.fixobjects))
{
  lines(encode.prog.fixobjects[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")
dev.off()



# Plot Results (tiff)
tiff(file = here("figures","figure5_objectmediated_progeny.tiff"), width = 9, height = 9, units = "in", res = 300)

par(mfrow=c(2,2))
# (decoding error)
cc <- brewer.pal(length(lambda),'Set1')
plot(wf.prog.n300.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="a")
for (i in 1:length(decode.prog.varlambda))
{
  lines(decode.prog.varlambda[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="b")
for (i in 1:length(decode.prog.fixobjects))
{
  lines(decode.prog.fixobjects[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")

# (encoding error)
plot(wf.prog.n300.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="c")
for (i in 1:length(encode.prog.varlambda))
{
  lines(encode.prog.varlambda[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n300.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny1,col=c(1,cc),pch=c(NA,rep(20,length(lambda))),lty=c(4,rep(NA,length(lambda))),lwd=c(2,rep(NA,length(lambda))),bty="n")

plot(wf.prog.n1000.m01$d2,type="n",lty=4,log="xy",ylab="Probability P(k) of number progeny >= k",xlab="k",col="darkgrey",main="d")
for (i in 1:length(encode.prog.fixobjects))
{
  lines(encode.prog.fixobjects[[i]]$d2, type="b",pch=20,col=cc[i])
}
lines(wf.prog.n1000.m01$d2,lty=4,lwd=2)
legend("bottomleft",legend=legItems.progeny2,col=c(1,cc),pch=c(NA,rep(20,length(N))),lty=c(4,rep(NA,length(N))),lwd=c(2,rep(NA,length(N))),bty="n")
dev.off()

