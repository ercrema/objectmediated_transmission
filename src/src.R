
### Utility Functions ###

#Utility Function for counting cases
instances <- function(x,cases)
{
  x <- c(x,cases)
  return(table(x) - 1)
}

#Compute turnover rate from frequency matrix
turnover <- function(mat,top)
{
  z<-matrix(NA,nrow=c(nrow(mat)-1),ncol=top)
  for (i in 1:c(nrow(mat)-1))
  {
    tmp.top = min(c(top,sum(mat[i,]>0),sum(mat[i+1,]>0)))
    t1 = names(sort(mat[i,],decreasing=TRUE))
    t2 = names(sort(mat[i+1,],decreasing=TRUE))
    z[i,1:tmp.top]=sapply(1:tmp.top,function(x,t1,t2){return(x-length(intersect(t1[1:x],t2[1:x])))},t1=t1,t2=t2)
  }
return(z)
}

#Compute add alpha values on color (from mages/add.alpha.R)
add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}


### Object mediated transmission ###
#
# mu.e ... encoding error rate
# mu.d ... decoding error rate
# lambda ... production rate
# N ... population size (number of mental templates)
# ta,tb,tc,td ... constants for turn-over rate analysis
# top ... top ranks used for turn-over rate analysis
# warmup ... number of discarded simulation runs
# timesteps ... number of timesteps
# output ... output of analysis:
#	* sumstat... computes diversity and richness for the final timestep
#	* turnover ... carries out turnover rate analysis
#	* progeny ... creates a progeny distribution

objTr <- function(N,mu.e,mu.d,lambda,timesteps,output=c("sumstat","turnover","progeny"),warmup=0,ta=0.55,tb=0.86,tc=0.130,td=1.38,top=NA,verbose=F)
{
	# Initial Setup
	mental.templates = 1:N
	trait.counter = N+1

	# Type specific setup
	if (output=="turnover")
	{
		rawList <- vector("list",length=timesteps-warmup)	
	}	
	if (output=="progeny")
	{
		rawList <- rep(NA,N*(timesteps-warmup)*lambda*2)	
		start = 1
		end = 0
	}
	
	if (verbose){pb <-txtProgressBar(min=1,max=timesteps,style=3)}

	# Initial production
	samplePool = rep(mental.templates,times=rpois(N,lambda))

	for (x in 1:timesteps)
	{

		if (verbose){setTxtProgressBar(pb,x)}

		#transmission
		mental.templates<-sample(size=N,x=samplePool,replace=TRUE)

		#Decoding Error
		index<-which(runif(N)<mu.d)
		if (length(index)>0)
		{
			newTraits<-trait.counter:c(trait.counter+length(index)-1)
			mental.templates[index]=newTraits
			trait.counter=max(newTraits)+1
		}

		#production
		samplePool = rep(mental.templates,times=rpois(N,lambda))

		#Encoding Error
		index<-which(runif(length(samplePool))<mu.e)
		if (length(index)>0)
		{
			newTraits<-trait.counter:c(trait.counter+length(index)-1)
			samplePool[index]=newTraits
			trait.counter=max(newTraits)+1
		}
		
		#store
		if (x>warmup&output=="turnover")
		{
			rawList[[x-warmup]]=samplePool
		}

		if (x>warmup&output=="progeny")
		{
			end = end + length(samplePool)
			rawList[start:end] = samplePool
			start = end + 1
		}

	}
	
	if(verbose){close(pb)}

	# Output calculation
	if (output=="sumstat")
	{
		if (verbose) {print("Computing Diversity and Richness ...")}
		p=table(samplePool)/length(samplePool)
		return(list(div=1-sum(p^2),k=length(p)))
	}

	if (output=="turnover")
	{
		if (verbose) {print("Computing turnover rate ...")}
		cases = unique(unlist(rawList))
		rawMatrix <- t(sapply(rawList,instances,cases=cases))
		#TurnoverRate
		if (is.na(top)) {top = min(apply(rawMatrix,1,function(x){sum(x>0)}))} 
		zMatrix=turnover(rawMatrix,top=top)
		z = apply(zMatrix,2,mean,na.rm=TRUE) #average turn-over per top list of size y
		z[which(is.nan(z))]=NA
		z[which(z==0)]=NA
		
		z.frame=data.frame(y=1:top,obs.z=z)

		aN = td * (mu^ta) * N^tc 
		aNp = aN*lambda 

		z.frame$exp.z.N= (aN*(z.frame$y)^tb)/2 
		z.frame$exp.z.Np= (aNp*(z.frame$y)^tb)/2 #see page 2 on Giometto and Evans on defintion of z

		z.frame$zfitN= z.frame$obs.z/aN
		z.frame$zfitNp= z.frame$obs.z/aNp

		print(z.frame)
		modN=lm(log(zfitN)~log(y),data=z.frame)
		modNp=lm(log(zfitNp)~log(y),data=z.frame)

		bN=as.numeric(coefficients(modN)[2])
		bNp=as.numeric(coefficients(modNp)[2])

		return(list(bN=bN,bNp=bNp,z.frame=z.frame))
	}
	if (output=="progeny")
	{
		if(verbose){print("computing progeny distribution ...")}
		sumTraits=table(rawList)
		#compute the frequencies of each variant
		k=as.numeric(names(table(sumTraits))) #compute the number of instances, e.g. "4" means a variant that had 4 cases in rawMatrix
		qk=table(sumTraits) #compute the matching frequencies, e.g. if "4" is matched to 12, it means that 12 variants had 4 cases in rawMatrix
		d=data.frame(instances=as.numeric(k),nvariants=as.numeric(qk)) 
		cumul<-rev(cumsum(d[nrow(d):1,]$nvariants)/sum(d[,2]))
		d2<-data.frame(k=rev(d[nrow(d):1,]$instances),qkp=cumul)
		return(list(d=d,d2=d2))
	}
}

### Wright-Fisher model ###
#
# mu ... error rate
# lambda ... production rate
# N ... population size (number of mental templates)
# ta,tb,tc,td ... constants for turn-over rate analysis
# top ... top ranks used for turn-over rate analysis
# warmup ... number of discarded simulation runs
# timesteps ... number of timesteps
# output ... output of analysis:
#	* sumstat... computes diversity and richness for the final timestep
#	* turnover ... carries out turnover rate analysis
#	* progeny ... creates a progeny distribution

wf  <- function(N,mu,timesteps,output=c("sumstat","turnover","progeny"),warmup=0,ta=0.55,tb=0.86,tc=0.130,td=1.38,top=NA,verbose=FALSE)
  {
	
	# Initial Setup
	mental.templates = 1:N
	trait.counter = N+1


	# Type specific setup
	if (output=="turnover")
	{
		rawList <- vector("list",length=timesteps-warmup)	
	}	
	if (output=="progeny")
	{
		rawList <- rep(NA,N*(timesteps-warmup))	
		start = 1
		end = 0
	}
	
	if (verbose){pb <-txtProgressBar(min=1,max=timesteps,style=3)}

	for (x in 1:timesteps)
	{

		if (verbose){setTxtProgressBar(pb,x)}
		#transmission
		mental.templates<-sample(size=N,x=mental.templates,replace=TRUE)

		#Transmission Error
		index<-which(runif(N)<mu)
		if (length(index)>0)
		{
			newTraits<-trait.counter:c(trait.counter+length(index)-1)
			mental.templates[index]=newTraits
			trait.counter=max(newTraits)+1
		}


		samplePool = mental.templates
		
		#store
		if (x>warmup&output=="turnover")
		{
			rawList[[x-warmup]]=samplePool
		}

		if (x>warmup&output=="progeny")
		{
			end = end + length(samplePool)
			rawList[start:end] = samplePool
			start = end + 1
		}

	}

	if(verbose){close(pb)}
	
	# Output calculation
	if (output=="sumstat")
	{
		if(verbose){print("computing diversity and richness...")}
		p=table(samplePool)/length(samplePool)
		return(list(div=1-sum(p^2),k=length(p)))
	}

	if (output=="turnover")
	{
		if(verbose){print("computing turnover rate...")}
		cases = unique(unlist(rawList))
		rawMatrix <- t(sapply(rawList,instances,cases=cases))
		#TurnoverRate
		if (is.na(top)) {top = min(apply(rawMatrix,1,function(x){sum(x>0)}))} 
		zMatrix=turnover(rawMatrix,top=top)
		z = apply(zMatrix,2,mean) #average turn-over per top list of size y
		z.frame=data.frame(y=1:top,obs.z=z)

		aN = td * (mu^ta) * N^tc 

		z.frame$exp.z.N= (aN*(z.frame$y)^tb)/2 
		z.frame$zfitN= z.frame$obs.z/aN
		modN=lm(log(zfitN)~log(y),data=z.frame)
		bN=as.numeric(coefficients(modN)[2])

		return(list(b=bN,z.frame=z.frame))
	}
	if (output=="progeny")
	{
		if(verbose){print("computing progeny distribution...")}
		sumTraits=table(rawList)
		#compute the frequencies of each variant
		k=as.numeric(names(table(sumTraits))) #compute the number of instances, e.g. "4" means a variant that had 4 cases in rawMatrix
		qk=table(sumTraits) #compute the matching frequencies, e.g. if "4" is matched to 12, it means that 12 variants had 4 cases in rawMatrix
		d=data.frame(instances=as.numeric(k),nvariants=as.numeric(qk)) 
		cumul<-rev(cumsum(d[nrow(d):1,]$nvariants)/sum(d[,2]))
		d2<-data.frame(k=rev(d[nrow(d):1,]$instances),qkp=cumul)
		return(list(d=d,d2=d2))
	}


}

