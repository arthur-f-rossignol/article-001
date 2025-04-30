### V5 model: V3 but adding (i) a randomized df per dataset; (ii) an aletrenative tmodel with mu fixed at 0 (=> estimatedtbis) and (iii) set.seed(j+i*N,kind="default",normal.kind="default",sample.kind="default")
#### adapted to be launched laso sequentiually on the cluster (with parallelize<-FALSE)
#### then needs runMCMC_btadjust without clearCompiled on Linux - otherwise DICs are not claculated
#### very strange: on Linux, Rdata are much bigger than on Windows!!!
### loading packages
library(coda)
#library(rjags)
#library(R2WinBUGS)
library(nimble); 
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE); nimbleOptions(MCMCorderPosteriorPredictiveSamplersLast = FALSE); nimbleOptions(enableDerivs = TRUE); library(nimbleHMC)
library(nlme)
library(runMCMCbtadjust)
library(parallel)
library(abind)


## meta data
### number of "repetitions" of cover-type calculations
rep<-1000
rep_to_add<-0

### number of samples for each repetitions
N<-500

### data set size
dataSize<-10
Xresampled<-TRUE
True.slope<-1
sd.X<-1

### whether to parallelize or not
parallelize<-TRUE


### MCMC parameters
Nchains<-4
Nchains.seq<-4 ### number of chains in case parallelize==FALSE must be at most Nchains
	Nchains.seq<-min(Nchains,Nchains.seq)

Rhat.max<-1.05
neff.min<-1000
Npvalues<-5 #must be less than neff.min

name<-paste0("rep",rep,"_N",N,"_dataSize",dataSize,"_",Nchains,"chains","_Rhat",Rhat.max,"_neff",neff.min,"_Xresampled",Xresampled)
name_reduced<-paste0("_N",N,"_dataSize",dataSize,"_",Nchains,"chains","_Rhat",Rhat.max,"_neff",neff.min,"_Xresampled",Xresampled)

name<-gsub("\\.","p",name)




## Nimble code that will be used:

{

### writing different nimble models
 
 
  modelCode_t <- nimbleCode(
  {
    
    #hyperparams
	## Priors
	mu~dnorm(0,sd=5.0)
	#activated to force variances to be finite
	log_df~T(dt(0,1.0,1.0),log(2.0),)
	#log_df~dt(0,1.0,1.0)
	
	
	## transformations of variables:
	df<-exp(log_df)

	#likelihood
	for (i in 1:n)
	{
		y[i]~dt(mu,sigma=SE[i],df=df)
	}
  }
)


##model without mu
  modelCode_tbis <- nimbleCode(
  {
    
    #hyperparams
	## Priors
	#mu~dnorm(0,sd=5.0)
	#activated to force variances to be finite
	log_df~T(dt(0,1.0,1.0),log(2.0),)
	#log_df~dt(0,1.0,1.0)
	
	
	## transformations of variables:
	df<-exp(log_df)

	#likelihood
	for (i in 1:n)
	{
		y[i]~dt(0,sigma=SE[i],df=df)
	}
  }
)


modelCode_Gaussian <- nimbleCode(
  {
    
    #hyperparams
	## Priors
	mu~dnorm(0,sd=5.0)
	
	
	#likelihood
	for (i in 1:n)
	{
		y[i]~dnorm(mu,sd=SE[i])
	}
  }
)


### code for DIC calculation in parallele:
### adapted to cases where some chains are removed because they do not move

calculations.for.modalDIC.parallel<-expression(
  {
     Nactivechains<-length(samplesList.temp)
    ##first, send to each cluster the set of sample parameters it will have to treat in a matrix format
    for (j.EC in 1:Nactivechains)
        {
          samples.to.treat.EC<-as.matrix(samplesList.temp[[j.EC]])
          parallel::clusterExport(cl[j.EC], "samples.to.treat.EC",envir=environment())
        }
    ## second, running calculations within each cluster with the parallel::clusterEvalQ function
    out1 <- parallel::clusterEvalQ(cl[1:Nactivechains], {
      Model1.EC<-Model[[1]]
      names.samples.EC<-dimnames(samples.to.treat.EC)[[2]]
    
      ## third preparing the names of variables to calculate on:
      varNames.EC<-CModel[[1]]$getVarNames()
      DatavarNames.EC<-names(data)
      notDatavarNames.EC<-setdiff(varNames.EC,DatavarNames.EC)
      
      ## fourth writing and compiling the nimbleFunction we will use:
    logProbCalc.EC <- nimbleFunction(
        setup = function(model,names.ref.list,notDatavarNames,DatavarNames) {
        },
    run = function(P = double(1)) {
        values(model,names.ref.list) <<- P
        model$calculate(notDatavarNames)
        return(model$calculate(DatavarNames))
        returnType(double(0))
    })
    logProbCalcPrepared.EC <- logProbCalc.EC(Model1.EC, names.samples.EC, notDatavarNames.EC, DatavarNames.EC)
    ClogProbCalcPrepared.EC <- compileNimble(Model1.EC, logProbCalcPrepared.EC)

      ## fifth, running through all the samples in a sapply function to obtain the logLikelihoods corresponding to each set of parameters:
          logLiks<-sapply(1:(dim(samples.to.treat.EC)[1]),function(toto) 
      {
      ClogProbCalcPrepared.EC$logProbCalcPrepared.EC$run(samples.to.treat.EC[toto,])
      })
    
          return(logLiks)
          gc(verbose = FALSE)
        })
    
    logLiks.EC<-unlist(out1)
    
    
    ## sixth: calculating DICs and estimation of numbers of parameters - outside of clusters for the first type of DIC and having to go back yo one cluster - here taken as the first - to do the calculation of logLikelihood on the mean of parameters as required by the classical formula of DIC:
    
    #mode type DIC; cf. Celeux et al. 2006 Bayesian analysis
    DIC.mode.EC<--4*mean(logLiks.EC)+2*max(logLiks.EC)
    p.DIC.mode.EC<--2*mean(logLiks.EC)+2*max(logLiks.EC)
    
    #calculation of classical DIC; cf. Celeux et al. 2006 Bayesian analysis
    samples.to.treat.EC<-colMeans(as.matrix(samplesList.temp))
    parallel::clusterExport(cl[1], "samples.to.treat.EC",envir=environment())
    out1 <- parallel::clusterEvalQ(cl[1], {
          logLiks.EC<-ClogProbCalcPrepared.EC$logProbCalcPrepared.EC$run(samples.to.treat.EC)
          return(logLiks.EC)
          gc(verbose = FALSE)
        })
    logLiks.meanparams.EC<-unlist(out1)
    DIC.EC<--4*mean(logLiks.EC)+2*logLiks.meanparams.EC
    p.DIC.EC<--2*mean(logLiks.EC)+2*logLiks.meanparams.EC
    
    list(DIC.mode=DIC.mode.EC,p.DIC.mode=p.DIC.mode.EC,DIC=DIC.EC,p.DIC=p.DIC.EC)
    
  }
)




calculations.for.modalDIC<-expression(
  {
    
    Model1.EC<-Model[[1]]
    
    ## first preparing the sampled parameters in a matrix format; uses the as.matrix function specific to mcmc.list objects; also stocking names of the parameters
    samples.List.matrix.EC<-as.matrix(samplesList.temp)
    names.samples.EC<-dimnames(samples.List.matrix.EC)[[2]]
    
    ## second preparing the names of variables to calculate on:
    varNames.EC<-Model1.EC$getVarNames()
    DatavarNames.EC<-names(data)
    notDatavarNames.EC<-setdiff(varNames.EC,DatavarNames.EC)
    
    ## third writing and compiling the nimbleFunction we will use:
    logProbCalc.EC <- nimbleFunction(
        setup = function(model,names.ref.list,notDatavarNames,DatavarNames) {
        },
    run = function(P = double(1)) { ##NB: double(1) means this if of double type and has one dimension
        values(model,names.ref.list) <<- P
        model$calculate(notDatavarNames)
        return(model$calculate(DatavarNames))
        returnType(double(0))
    })
    logProbCalcPrepared.EC <- logProbCalc.EC(Model1.EC, names.samples.EC, notDatavarNames.EC, DatavarNames.EC)
    ClogProbCalcPrepared.EC <- compileNimble(Model1.EC, logProbCalcPrepared.EC)

    
    ## fourth, running through all the samples in a sapply function to obtain the logLikelihoods corresponding to each set of parameters:
    logLiks.EC<-sapply(1:(dim(samples.List.matrix.EC)[1]),function(toto) 
      {
      ClogProbCalcPrepared.EC$logProbCalcPrepared.EC$run(samples.List.matrix.EC[toto,])
      })
    ## fifth: calculating DICs and estimation of numbers of parameters:
    
    #mode type DIC; cf. Celeux et al. 2006 Bayesian analysis
    DIC.mode.EC<--4*mean(logLiks.EC)+2*max(logLiks.EC)
    p.DIC.mode.EC<--2*mean(logLiks.EC)+2*max(logLiks.EC)
    
    #calculation of classical DIC; cf. Celeux et al. 2006 Bayesian analysis
    logLiks.meanparams.EC<-ClogProbCalcPrepared.EC$logProbCalcPrepared.EC$run(colMeans(samples.List.matrix.EC))
    
    DIC.EC<--4*mean(logLiks.EC)+2*logLiks.meanparams.EC
    p.DIC.EC<--2*mean(logLiks.EC)+2*logLiks.meanparams.EC
    
    list(DIC.mode=DIC.mode.EC,p.DIC.mode=p.DIC.mode.EC,DIC=DIC.EC,p.DIC=p.DIC.EC)
    }
)


ModelInits_t <- function()
{
  mu<-rnorm(1)
log_sd<-rnorm(1)
log_df<-rnorm(1,6)

		list (mu=mu, log_df=log_df)
}


ModelInits_tbis <- function()
{
  mu<-rnorm(1)
log_sd<-rnorm(1)
log_df<-rnorm(1,6)

		list (log_df=log_df)
}


ModelInits_Gaus <- function()
{
  mu<-rnorm(1)
log_sd<-rnorm(1)


		list (mu=mu)
}
}

try(setwd("E:\\Dossier Frederic\\Articles\\PASSIFOR2_Spatial_Bias\\DATA_ANALYSIS\\Results_Analysis"))

i<-0
mean_estimates<-array(NA,dim=c(rep,N,2))
se_estimates<-array(NA,dim=c(rep,N,2))
res_Bayest_Intercept<-NULL
res_Bayestbis_Intercept<-NULL
res_BayesGaus_Intercept<-NULL
res_MCMC_dfs_t_Intercept<-NULL
res_MCMC_dfs_tbis_Intercept<-NULL
res_Bayest_Slope<-NULL
res_Bayestbis_Slope<-NULL
res_BayesGaus_Slope<-NULL
res_MCMC_dfs_t_Slope<-NULL
res_MCMC_dfs_tbis_Slope<-NULL

data_cover_rate_gauss_semult<-array(NA,dim=c(rep,N,2))
data_cover_rate_theort1_semult<-array(NA,dim=c(rep,N,2))
data_cover_rate_theort1_sesdt<-array(NA,dim=c(rep,N,2))
data_cover_rate_theort2_semult<-array(NA,dim=c(rep,N,2))
data_cover_rate_theort2_sesdt<-array(NA,dim=c(rep,N,2))
data_cover_rate_estimatedt_semult<-array(NA,dim=c(rep,N,Npvalues,2))
data_cover_rate_estimatedt_sesdt<-array(NA,dim=c(rep,N,Npvalues,2))
data_cover_rate_estimatedtbis_semult<-array(NA,dim=c(rep,N,Npvalues,2))
data_cover_rate_estimatedtbis_sesdt<-array(NA,dim=c(rep,N,Npvalues,2))
data_cover_rate_bestmodel_semult<-array(NA,dim=c(rep,N,Npvalues,2))
data_cover_rate_bestmodel_sesdt<-array(NA,dim=c(rep,N,Npvalues,2))



file.to.load<-list.files()[regexpr("Results_",list.files())==1&regexpr("_InterceptandSlopemodel_V5.RData",list.files())>0&regexpr("_simplifiedanalysiscoverrate",list.files())>0  & regexpr(name_reduced,list.files())>0]
if (length(file.to.load)>0){
load(file.to.load)
}

rep<-500
rep_to_add<-500

if (i<(rep+rep_to_add))
{
		added_elem<-array(NA,dim=c(rep+rep_to_add-i,N,2))
		mean_estimates<-abind(mean_estimates,added_elem,along=1)
		se_estimates<-abind(se_estimates,added_elem,along=1)
		data_cover_rate_gauss_semult<-abind(data_cover_rate_gauss_semult,added_elem,along=1)
		data_cover_rate_theort1_semult<-abind(data_cover_rate_theort1_semult,added_elem,along=1)
		data_cover_rate_theort1_sesdt<-abind(data_cover_rate_theort1_sesdt,added_elem,along=1)
		data_cover_rate_theort2_semult<-abind(data_cover_rate_theort2_semult,added_elem,along=1)
		data_cover_rate_theort2_sesdt<-abind(data_cover_rate_theort2_sesdt,added_elem,along=1)
		added_elem<-array(NA,dim=c(rep+rep_to_add-i,N,Npvalues,2))
		data_cover_rate_estimatedt_semult<-abind(data_cover_rate_estimatedt_semult,added_elem,along=1)
		data_cover_rate_estimatedt_sesdt<-abind(data_cover_rate_estimatedt_sesdt,added_elem,along=1)
		data_cover_rate_estimatedtbis_semult<-abind(data_cover_rate_estimatedtbis_semult,added_elem,along=1)
		data_cover_rate_estimatedtbis_sesdt<-abind(data_cover_rate_estimatedtbis_sesdt,added_elem,along=1)
		data_cover_rate_bestmodel_semult<-abind(data_cover_rate_bestmodel_semult,added_elem,along=1)
		data_cover_rate_bestmodel_sesdt<-abind(data_cover_rate_bestmodel_sesdt,added_elem,along=1)
}

iInit<-i+1
for (i in iInit:(rep+rep_to_add))

{ ended_i<-FALSE
set.seed(1+i*N,kind="default",normal.kind="default",sample.kind="default")
X<-rnorm(dataSize,mean=0,sd=sd.X)

for (j in 1:N)

	{ 
	
	if (Xresampled) {### original version which is not so problematic in case Xresampled
		set.seed(j+i*N,kind="default",normal.kind="default",sample.kind="default")
		} else {### to avoid colinearity between data (below) and X (above)
		set.seed(j+1+i*N,kind="default",normal.kind="default",sample.kind="default")
		}
	if (Xresampled) {X<-rnorm(dataSize,mean=0,sd=sd.X)}
	data<-rnorm(dataSize,0+True.slope*X)
	lModel<-lm(data~X)
	mean_estimates[i,j,]<-summary(lModel)$coef[,1]
	se_estimates[i,j,]<-summary(lModel)$coef[,2]
	}


## first calculus for Intercept	


### estimation of Bayesian models (t and Gaus) on these N mean results
	
	params <- c("mu","log_df")

	Inits<-lapply(1:Nchains,function(x){ModelInits_t()})

	ModelConsts=list(n=N, SE=se_estimates[i,,1])
	ModelData=list(y=mean_estimates[i,,1])

if (parallelize)
{
temp<-runMCMC_btadjust(code=modelCode_t, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains, inits = {if(Nchains==1){Inits[3]} else {Inits}}, conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC.parallel, parallelize=TRUE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
}	 else {




if (Nchains.seq==1){
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},length(params)),0.95),digits=3,scientific=FALSE))
temp<-runMCMC_btadjust(code=modelCode_t, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=1, inits = Inits[1], conv.max=Gew.Max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					}
					else
					{
					
					temp<-runMCMC_btadjust(code=modelCode_t, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains.seq, inits = Inits[1:Nchains.seq], conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					
					}


}
	
	res_Bayest_Intercept[[i]]<-c(attributes(temp)[c("final.params","final.diags")],summary=list(summary(temp)))
	tempt<-temp
	
	



	params <- c("log_df")

	Inits<-lapply(1:Nchains,function(x){ModelInits_tbis()})

	ModelConsts=list(n=N, SE=se_estimates[i,,1])
	ModelData=list(y=mean_estimates[i,,1]-0.0)


if (parallelize)
{
temp<-runMCMC_btadjust(code=modelCode_tbis, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains, inits = {if(Nchains==1){Inits[3]} else {Inits}}, conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(parallelize=TRUE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
}	 else {

if (Nchains.seq==1){
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},length(params)),0.95),digits=3,scientific=FALSE))
temp<-runMCMC_btadjust(code=modelCode_tbis, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=1, inits = Inits[1], conv.max=Gew.Max,neff.min=neff.min,
					control.MCMC=list(parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					}
					else
					{
					
					temp<-runMCMC_btadjust(code=modelCode_tbis, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains.seq, inits = Inits[1:Nchains.seq], conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					
					}



}
	res_Bayestbis_Intercept[[i]]<-c(attributes(temp)[c("final.params","final.diags")],summary=list(summary(temp)))
	temptbis<-temp
	
	

	params <- c("mu")

	Inits<-lapply(1:Nchains,function(x){ModelInits_Gaus()})

	ModelConsts=list(n=N, SE=se_estimates[i,,1])
	ModelData=list(y=mean_estimates[i,,1])



if (parallelize)
{
temp<-runMCMC_btadjust(code=modelCode_Gaussian, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains, inits = {if(Nchains==1){Inits[3]} else {Inits}}, conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC.parallel, parallelize=TRUE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
}	 else {

if (Nchains.seq==1){
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},length(params)),0.95),digits=3,scientific=FALSE))
temp<-runMCMC_btadjust(code=modelCode_Gaussian, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=1, inits = Inits[1], conv.max=Gew.Max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					}
					else
					{
					
					temp<-runMCMC_btadjust(code=modelCode_Gaussian, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains.seq, inits = Inits[1:Nchains.seq], conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					
					}



}
	res_BayesGaus_Intercept[[i]]<-c(attributes(temp)[c("final.params","final.diags")],summary=list(summary(temp)))

	
	
	
	
	
	#### command to remove temporary files & directories associated to nimble generated code:
	path<-"C:\\Users\\FGOSSE~1\\AppData\\Local\\Temp\\"
		for(file in list.files("C:\\Users\\FGOSSE~1\\AppData\\Local\\Temp\\")) {
		if (!is.na(file.info(file.path(path,file))$mtime))
			{if (length(list.files(file.path(path,file)))==1)
				{if (list.files(file.path(path,file))=="nimble_generatedCode")
					{
					if (abs(difftime(file.info(file.path(path,file))$mtime,Sys.time(),units="secs"))>600& abs(difftime(file.info(file.path(path,file))$mtime,Sys.time(),units="secs"))<6000) {
					 try(file.remove(file.path(path, file)))
					 ### added the following due to: https://stackoverflow.com/questions/22122968/how-to-remove-the-file-when-permisson-denied-in-r
					 unlink(file.path(path, file),recursive=TRUE)

						}
					}
				}
			}
		}

	
	set.seed(j+i*N,kind="default",normal.kind="default",sample.kind="default")
	MCMC_dfs_t<-exp(unlist(tempt[,"log_df"]))
	if (length(MCMC_dfs_t)>neff.min)
	{MCMC_dfs_t<-sample(exp(unlist(tempt[,"log_df"])),neff.min, replace=FALSE)}
	res_MCMC_dfs_t_Intercept[[i]]<-MCMC_dfs_t



	set.seed(j+i*N,kind="default",normal.kind="default",sample.kind="default")
	MCMC_dfs_tbis<-exp(unlist(temptbis))
	if (length(MCMC_dfs_tbis)>neff.min)
	{MCMC_dfs_tbis<-sample(exp(unlist(temptbis)),neff.min, replace=FALSE)}
	res_MCMC_dfs_tbis_Intercept[[i]]<-MCMC_dfs_tbis

	#### calculation of cover rates with the different methods



 for (j in 1:N)
	{ 
	data_cover_rate_gauss_semult[i,j,1]<-pnorm(0,mean=mean_estimates[i,j,1],sd=se_estimates[i,j,1])
	data_cover_rate_theort2_semult[i,j,1]<-pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=dataSize-length(coef(lModel))-1)
	data_cover_rate_theort1_semult[i,j,1]<-pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=dataSize-length(coef(lModel)))
	data_cover_rate_theort2_sesdt[i,j,1]<-pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt((dataSize-length(coef(lModel))-1-2)/(dataSize-length(coef(lModel))-1))),df=dataSize-length(coef(lModel))-1)
	data_cover_rate_theort1_sesdt[i,j,1]<-pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt((dataSize-length(coef(lModel))-2)/(dataSize-length(coef(lModel))))),df=dataSize-length(coef(lModel)))
	
	
	temp<-sapply(MCMC_dfs_t,function(x){pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedt_semult[i,j,,1]<-c(temp[1],temprand,mean(temp),pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=mean(MCMC_dfs_t)),pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=median(MCMC_dfs_t)))
	if((res_BayesGaus_Intercept[[i]])$final.params$extra$DIC.mode>(res_Bayest_Intercept[[i]])$final.params$extra$DIC.mode)
		{
		data_cover_rate_bestmodel_semult[i,j,,1]<-c(temp[1],temprand,mean(temp),pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=mean(MCMC_dfs_t)),pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=median(MCMC_dfs_t)))
		} else {
		data_cover_rate_bestmodel_semult[i,j,,1]<-rep(pnorm(0,mean_estimates[i,j,1],se_estimates[i,j,1]),5)
		
		}
	
	temp<-sapply(MCMC_dfs_tbis,function(x){pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedtbis_semult[i,j,,1]<-c(temp[1],temprand,mean(temp),pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=mean(MCMC_dfs_tbis)),pt(0-mean_estimates[i,j,1]/se_estimates[i,j,1],df=median(MCMC_dfs_tbis)))
	
	
	temp<-sapply(MCMC_dfs_t,function(x){pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(x/(x-2))),df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedt_sesdt[i,j,,1]<-c(temp[1],temprand,mean(temp),pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(mean(MCMC_dfs_t)/(mean(MCMC_dfs_t)-2))),df=mean(MCMC_dfs_t)),pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(median(MCMC_dfs_t)/(median(MCMC_dfs_t)-2))),df=median(MCMC_dfs_t)))
	if((res_BayesGaus_Intercept[[i]])$final.params$extra$DIC.mode>(res_Bayest_Intercept[[i]])$final.params$extra$DIC.mode)
		{
		data_cover_rate_bestmodel_sesdt[i,j,,1]<-c(temp[1],temprand,mean(temp),pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(mean(MCMC_dfs_t)/(mean(MCMC_dfs_t)-2))),df=mean(MCMC_dfs_t)),pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(median(MCMC_dfs_t)/(median(MCMC_dfs_t)-2))),df=median(MCMC_dfs_t)))
		} else {
		data_cover_rate_bestmodel_sesdt[i,j,,1]<-rep(pnorm(0,mean_estimates[i,j,1],se_estimates[i,j,1]),5)
		}
		
	temp<-sapply(MCMC_dfs_tbis,function(x){pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(x/(x-2))),df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedtbis_sesdt[i,j,,1]<-c(temp[1],temprand,mean(temp),pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(mean(MCMC_dfs_tbis)/(mean(MCMC_dfs_tbis)-2))),df=mean(MCMC_dfs_tbis)),pt(0-mean_estimates[i,j,1]/(se_estimates[i,j,1]*sqrt(median(MCMC_dfs_tbis)/(median(MCMC_dfs_tbis)-2))),df=median(MCMC_dfs_tbis)))
	
	
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
## first calculus for Slope


### estimation of Bayesian models (t and Gaus) on these N mean results
	
	params <- c("mu","log_df")

	Inits<-lapply(1:Nchains,function(x){ModelInits_t()})

	ModelConsts=list(n=N, SE=se_estimates[i,,2])
	ModelData=list(y=mean_estimates[i,,2])

if (parallelize)
{
temp<-runMCMC_btadjust(code=modelCode_t, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains, inits = {if(Nchains==1){Inits[3]} else {Inits}}, conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC.parallel, parallelize=TRUE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
}	 else {




if (Nchains.seq==1){
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},length(params)),0.95),digits=3,scientific=FALSE))
temp<-runMCMC_btadjust(code=modelCode_t, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=1, inits = Inits[1], conv.max=Gew.Max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					}
					else
					{
					
					temp<-runMCMC_btadjust(code=modelCode_t, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains.seq, inits = Inits[1:Nchains.seq], conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					
					}


}
	
	res_Bayest_Slope[[i]]<-c(attributes(temp)[c("final.params","final.diags")],summary=list(summary(temp)))
	tempt<-temp
	
	



	params <- c("log_df")

	Inits<-lapply(1:Nchains,function(x){ModelInits_tbis()})

	ModelConsts=list(n=N, SE=se_estimates[i,,2])
	ModelData=list(y=mean_estimates[i,,2]-True.slope)


if (parallelize)
{
temp<-runMCMC_btadjust(code=modelCode_tbis, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains, inits = {if(Nchains==1){Inits[3]} else {Inits}}, conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(parallelize=TRUE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
}	 else {

if (Nchains.seq==1){
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},length(params)),0.95),digits=3,scientific=FALSE))
temp<-runMCMC_btadjust(code=modelCode_tbis, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=1, inits = Inits[1], conv.max=Gew.Max,neff.min=neff.min,
					control.MCMC=list(parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					}
					else
					{
					
					temp<-runMCMC_btadjust(code=modelCode_tbis, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains.seq, inits = Inits[1:Nchains.seq], conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					
					}



}
	res_Bayestbis_Slope[[i]]<-c(attributes(temp)[c("final.params","final.diags")],summary=list(summary(temp)))
	temptbis<-temp
	
	

	params <- c("mu")

	Inits<-lapply(1:Nchains,function(x){ModelInits_Gaus()})

	ModelConsts=list(n=N, SE=se_estimates[i,,2])
	ModelData=list(y=mean_estimates[i,,2])



if (parallelize)
{
temp<-runMCMC_btadjust(code=modelCode_Gaussian, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains, inits = {if(Nchains==1){Inits[3]} else {Inits}}, conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC.parallel, parallelize=TRUE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
}	 else {

if (Nchains.seq==1){
Gew.Max<-as.double(format(quantile(sapply(1:100000,function(x,N){max(abs(rnorm(N)))},length(params)),0.95),digits=3,scientific=FALSE))
temp<-runMCMC_btadjust(code=modelCode_Gaussian, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=1, inits = Inits[1], conv.max=Gew.Max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					}
					else
					{
					
					temp<-runMCMC_btadjust(code=modelCode_Gaussian, constants = ModelConsts, data = ModelData,params.save=params,params.conv=params,
					niter.min=20000,niter.max=Inf,nburnin.min=10000,nburnin.max=Inf,thin.min=10,thin.max=Inf,
					Nchains=Nchains.seq, inits = Inits[1:Nchains.seq], conv.max=Rhat.max,neff.min=neff.min,
					control.MCMC=list(extraCalculations=calculations.for.modalDIC, parallelize=FALSE),
					control=list(time.max=3000,round.thinmult=TRUE,print.diagnostics=TRUE,Ncycles.target=3,check.convergence.firstrun=TRUE))
					
					}



}
	res_BayesGaus_Slope[[i]]<-c(attributes(temp)[c("final.params","final.diags")],summary=list(summary(temp)))

	
	
	set.seed(j+i*N,kind="default",normal.kind="default",sample.kind="default")
	MCMC_dfs_t<-exp(unlist(tempt[,"log_df"]))
	if (length(MCMC_dfs_t)>neff.min)
	{MCMC_dfs_t<-sample(exp(unlist(tempt[,"log_df"])),neff.min, replace=FALSE)}
	res_MCMC_dfs_t_Slope[[i]]<-MCMC_dfs_t



	set.seed(j+i*N,kind="default",normal.kind="default",sample.kind="default")
	MCMC_dfs_tbis<-exp(unlist(temptbis))
	if (length(MCMC_dfs_tbis)>neff.min)
	{MCMC_dfs_tbis<-sample(exp(unlist(temptbis)),neff.min, replace=FALSE)}
	res_MCMC_dfs_tbis_Slope[[i]]<-MCMC_dfs_tbis

	#### calculation of cover rates with the different methods



 for (j in 1:N)
	{ 
	data_cover_rate_gauss_semult[i,j,2]<-pnorm(True.slope,mean=mean_estimates[i,j,2],sd=se_estimates[i,j,2])
	data_cover_rate_theort2_semult[i,j,2]<-pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=dataSize-length(coef(lModel))-1)
	data_cover_rate_theort1_semult[i,j,2]<-pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=dataSize-length(coef(lModel)))
	data_cover_rate_theort2_sesdt[i,j,2]<-pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt((dataSize-length(coef(lModel))-1-2)/(dataSize-length(coef(lModel))-1))),df=dataSize-length(coef(lModel))-1)
	data_cover_rate_theort1_sesdt[i,j,2]<-pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt((dataSize-length(coef(lModel))-2)/(dataSize-length(coef(lModel))))),df=dataSize-length(coef(lModel)))
	
	
	temp<-sapply(MCMC_dfs_t,function(x){pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedt_semult[i,j,,2]<-c(temp[1],temprand,mean(temp),pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=mean(MCMC_dfs_t)),pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=median(MCMC_dfs_t)))
	if((res_BayesGaus_Slope[[i]])$final.params$extra$DIC.mode>(res_Bayest_Slope[[i]])$final.params$extra$DIC.mode)
		{
		data_cover_rate_bestmodel_semult[i,j,,2]<-c(temp[1],temprand,mean(temp),pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=mean(MCMC_dfs_t)),pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=median(MCMC_dfs_t)))
		} else {
		data_cover_rate_bestmodel_semult[i,j,,2]<-rep(pnorm(True.slope,mean_estimates[i,j,2],se_estimates[i,j,2]),5)
		
		}
	
	temp<-sapply(MCMC_dfs_tbis,function(x){pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedtbis_semult[i,j,,2]<-c(temp[1],temprand,mean(temp),pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=mean(MCMC_dfs_tbis)),pt((True.slope-mean_estimates[i,j,2])/se_estimates[i,j,2],df=median(MCMC_dfs_tbis)))
	
	
	temp<-sapply(MCMC_dfs_t,function(x){pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(x/(x-2))),df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedt_sesdt[i,j,,2]<-c(temp[1],temprand,mean(temp),pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(mean(MCMC_dfs_t)/(mean(MCMC_dfs_t)-2))),df=mean(MCMC_dfs_t)),pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(median(MCMC_dfs_t)/(median(MCMC_dfs_t)-2))),df=median(MCMC_dfs_t)))
	if((res_BayesGaus_Slope[[i]])$final.params$extra$DIC.mode>(res_Bayest_Slope[[i]])$final.params$extra$DIC.mode)
		{
		data_cover_rate_bestmodel_sesdt[i,j,,2]<-c(temp[1],temprand,mean(temp),pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(mean(MCMC_dfs_t)/(mean(MCMC_dfs_t)-2))),df=mean(MCMC_dfs_t)),pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(median(MCMC_dfs_t)/(median(MCMC_dfs_t)-2))),df=median(MCMC_dfs_t)))
		} else {
		data_cover_rate_bestmodel_sesdt[i,j,,2]<-rep(pnorm(True.slope,mean_estimates[i,j,2],se_estimates[i,j,2]),5)
		}
		
	temp<-sapply(MCMC_dfs_tbis,function(x){pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(x/(x-2))),df=x)}    )
	temprand<-sample(temp,1)
	data_cover_rate_estimatedtbis_sesdt[i,j,,2]<-c(temp[1],temprand,mean(temp),pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(mean(MCMC_dfs_tbis)/(mean(MCMC_dfs_tbis)-2))),df=mean(MCMC_dfs_tbis)),pt((True.slope-mean_estimates[i,j,2])/(se_estimates[i,j,2]*sqrt(median(MCMC_dfs_tbis)/(median(MCMC_dfs_tbis)-2))),df=median(MCMC_dfs_tbis)))
	
	
	}

if (length(res_BayesGaus_Intercept)==i&length(res_Bayest_Intercept)==i&length(res_Bayestbis_Intercept)==i&length(res_BayesGaus_Slope)==i&length(res_Bayest_Slope)==i&length(res_Bayestbis_Slope)==i&sum(!is.na(data_cover_rate_bestmodel_sesdt[i,,,]))>0) {ended_i=TRUE}


if (ended_i)
{
	name<-paste0("rep",i-1,"_N",N,"_dataSize",dataSize,"_",Nchains,"chains","_Rhat",Rhat.max,"_neff",neff.min,"_Xresampled",Xresampled)
	name<-gsub("\\.","p",name)
	if (file.exists(paste0("Results_simplifiedanalysiscoverrate_",name,"_InterceptandSlopemodel_V5.RData"))) {
	  #Delete file if it exists
	  file.remove(paste0("Results_simplifiedanalysiscoverrate_",name,"_InterceptandSlopemodel_V5.RData"))
	}
	name<-paste0("rep",i,"_N",N,"_dataSize",dataSize,"_",Nchains,"chains","_Rhat",Rhat.max,"_neff",neff.min,"_Xresampled",Xresampled)
	name<-gsub("\\.","p",name)
	save(list=ls(),file=paste0("Results_simplifiedanalysiscoverrate_",name,"_InterceptandSlopemodel_V5.RData"))
}
				
}

### saving results

if (ended_i)
{
	name<-paste0("rep",i-1,"_N",N,"_dataSize",dataSize,"_",Nchains,"chains","_Rhat",Rhat.max,"_neff",neff.min,"_Xresampled",Xresampled)
	name<-gsub("\\.","p",name)
	if (file.exists(paste0("Results_simplifiedanalysiscoverrate_",name,"_InterceptandSlopemodel_V5.RData"))) {
	  #Delete file if it exists
	  file.remove(paste0("Results_simplifiedanalysiscoverrate_",name,"_InterceptandSlopemodel_V5.RData"))
	}
	name<-paste0("rep",i,"_N",N,"_dataSize",dataSize,"_",Nchains,"chains","_Rhat",Rhat.max,"_neff",neff.min,"_Xresampled",Xresampled)
	name<-gsub("\\.","p",name)
	save(list=ls(),file=paste0("Results_simplifiedanalysiscoverrate_",name,"_InterceptandSlopemodel_V5.RData"))
}
		