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

# Experiment 1a ----
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
# pspace.wf$hom = res.wf[1,]
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

# Experiment 1b ----

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
unlist(wf(N=rpois(1,n.objects),mu=mu,timesteps = 5000,output="sumstat"))
}
res.wf.hom = res.wf[1,]
res.wf.k = res.wf[2,]

res.decode = foreach (i = 1:nrow(pspace), .combine='cbind') %dopar% {
unlist(objTr(N=pspace$N[i],lambda=pspace$lambda[i],mu.e=0,mu.d=mu.d,timesteps = 5000,output="sumstat"))}
res.decode.hom=matrix(res.decode[1,],ncol=length(N),nrow=nsim)
res.decode.k=matrix(res.decode[2,],ncol=length(N),nrow=nsim)

res.encode = foreach (i = 1:nrow(pspace), .combine='cbind') %dopar% {
unlist(objTr(N=pspace$N[i],lambda=pspace$lambda[i],mu.e=mu.e,mu.d=0,timesteps = 5000,output="sumstat"))}
res.encode.hom=matrix(res.encode[1,],ncol=length(N),nrow=nsim)
res.encode.k=matrix(res.encode[2,],ncol=length(N),nrow=nsim)

#stop cluster
stopCluster(cl)

# Experiment 2 ----

# Expectation Under Neutrality (Figure)

# changing population size
wf.prog.n300.m01 = wf(N=300,mu=0.01,warmup=10000,timesteps=20000,output="progeny")
wf.prog.n1000.m01 = wf(N=1000,mu=0.01,warmup=10000,timesteps=20000,output="progeny")
wf.prog.n3000.m01 = wf(N=3000,mu=0.01,warmup=10000,timesteps=20000,output="progeny")

# changing innovation rate
wf.prog.n1000.m005 = wf(N=300,mu=0.005,warmup=10000,timesteps=20000,output="progeny")
wf.prog.n3000.m001 = wf(N=300,mu=0.001,warmup=10000,timesteps=20000,output="progeny")                                             

#
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

# Store Output into R Image File
save.image(file=here("simres.RData"))
