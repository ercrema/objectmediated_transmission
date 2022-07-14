### Utility Functions ----

#Counting cases number of defined cases
#x... vector
#cases ... unique cases to be counted
instances <- function(x,cases)
{
  x <- c(x,cases)
  return(table(x) - 1)
}

#Compute add alpha values on color (from mages/add.alpha.R)
add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}


### Object mediated transmission ----
#
# mu.e ... encoding error rate
# mu.d ... decoding error rate
# lambda ... production rate
# N ... population size (number of mental templates)
# warmup ... number of discarded simulation runs
# timesteps ... number of timesteps
# output ... output of analysis:
#	* sumstat... computes diversity and richness for the final timestep
#	* progeny ... creates a progeny distribution

objTr <- function(N,mu.e,mu.d,lambda,timesteps,output=c("sumstat","progeny"),warmup=0,verbose=F)
{
	# Initial Setup
	mental.templates = 1:N
	trait.counter = N+1

	# Type specific setup
	if (output=="progeny")
	{
		rawList.objects <- rep(NA,N*(timesteps-warmup)*lambda*2)	
		rawList.mental <- rep(NA,N*(timesteps-warmup))	
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
		if (x>warmup&output=="progeny")
		{
# 			end = end + length(samplePool)
# 			rawList[start:end] = samplePool
# 			start = end + 1
			rawList.objects = c(rawList.objects,samplePool)
			rawList.mental = c(rawList.mental,mental.templates)
		}

	}
	
	if(verbose){close(pb)}

	# Output calculation
	if (output=="sumstat")
	{
		if (verbose) {print("Computing Diversity and Richness ...")}
		p.objects = table(samplePool)/length(samplePool)
		p.mental = table(mental.templates)/length(mental.templates)
		return(list(div.obj=1-sum(p.objects^2),k.obj=length(p.objects),div.mental=1-sum(p.mental^2),k.mental=length(p.mental)))
	}

	if (output=="progeny")
	{
		if(verbose){print("computing progeny distribution ...")}
		sumTraits.objects=table(rawList.objects)
		sumTraits.mental=table(rawList.mental)
		#compute the frequencies of each variant
		k.objects=as.numeric(names(table(sumTraits.objects))) #compute the number of instances, e.g. "4" means a variant that had 4 cases in rawMatrix
		k.mental=as.numeric(names(table(sumTraits.mental))) #compute the number of instances, e.g. "4" means a variant that had 4 cases in rawMatrix
		qk.objects=table(sumTraits.objects) #compute the matching frequencies, e.g. if "4" is matched to 12, it means that 12 variants had 4 cases in rawMatrix
		qk.mental=table(sumTraits.mental) #compute the matching frequencies, e.g. if "4" is matched to 12, it means that 12 variants had 4 cases in rawMatrix
		d.objects=data.frame(instances=as.numeric(k.objects),nvariants=as.numeric(qk.objects)) 
		cumul.objects<-rev(cumsum(d.objects[nrow(d.objects):1,]$nvariants)/sum(d.objects[,2]))
		d2.objects<-data.frame(k=rev(d.objects[nrow(d.objects):1,]$instances),qkp=cumul.objects)

		d.mental=data.frame(instances=as.numeric(k.mental),nvariants=as.numeric(qk.mental)) 
		cumul.mental<-rev(cumsum(d.mental[nrow(d.mental):1,]$nvariants)/sum(d.mental[,2]))
		d2.mental<-data.frame(k=rev(d.mental[nrow(d.mental):1,]$instances),qkp=cumul.mental)

		return(list(d.objects=d.objects,d2.objects=d2.objects,d.mental=d.mental,d2.mental=d2.mental))
	}
}

### Wright-Fisher model ----
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

wf  <- function(N,mu,timesteps,output=c("sumstat","progeny"),warmup=0,verbose=FALSE)
  {
	
	# Initial Setup
	mental.templates = 1:N
	trait.counter = N+1


	# Type specific setup
	if (output=="progeny")
	{
		rawList  <- numeric()
# 		rawList <- rep(NA,N*(timesteps-warmup))	
# 		start = 1
# 		end = 0
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
		if (x>warmup&output=="progeny")
		{
# 			end = end + length(samplePool)
# 			rawList[start:end] = samplePool
# 			start = end + 1
			rawList = c(rawList,samplePool)
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

