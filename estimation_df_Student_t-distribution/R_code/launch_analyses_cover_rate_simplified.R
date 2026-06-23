##############################################################################
#### Program to perform automatic analyses of Rdata of Simulations techniques cover rate.....r
######################### will be put in dfmethod_xlsfile that should be created at the god format in terms of sheets and columns independently from this file
##############################################################################

#### program to be pasted/executed in R
######### !!!!!!!!!!!!!!!!!!! only when Rdata are not being updated - otherwise files that existed at one point in time may have disappeared

#####################Notes:
######### core results in Rdata will all be in the following arrays: cover_rate produced by the different methods:
# [5] "data_cover_rate_bestmodel_semult"         "data_cover_rate_bestmodel_sesdt"          "data_cover_rate_estimatedt_semult"        "data_cover_rate_estimatedt_sesdt"        
# [9] "data_cover_rate_estimatedtbis_semult"     "data_cover_rate_estimatedtbis_sesdt"      "data_cover_rate_gauss_semult"             "data_cover_rate_theort1_semult"          
#[13] "data_cover_rate_theort1_sesdt"            "data_cover_rate_theort2_semult"           "data_cover_rate_theort2_sesdt"                                           



#### preliminary parameters that will guide all the sequence: a priori nothing to parametrise

	### should we recalculate Model Comparisons (DICs) even if they exist in the current RData? Default to FALSE
	#ReCalculate_MC<-FALSE
	### should we rewrite the results of model comparisons (and potentially other info) in the xlsx file? Default to TRUE
	#ReWrite_MC<-TRUE
	### should we accept to change the RData files if necessary (if new DICs_result)? Default to TRUE
	#Update_RDatafiles<-TRUE
	### should we recalculate graphics even if they already exist? Default to FALSE
	#ReDo_Graphs<-FALSE
	

#names of files
	dfmethod_xlsfile_void<-"Results_estimation_df_methods_spatialanal_void.xlsx"
	dfmethod_xlsfile_final<-"Results_estimation_df_methods_spatialanal_V5.xlsx"
	
###############################################################################



#### DEFAULT WORKING DIRECTORIES
WD.results<-"E:\\Dossier Frederic\\Articles\\PASSIFOR2_Spatial_Bias\\DATA_ANALYSIS\\Results_Analysis"



### loading libraries:
library(openxlsx)


### loading xlsx workbook:
setwd(WD.results)
dfmethod_workbook<-loadWorkbook(file=dfmethod_xlsfile_void)

#######################################################################"""""
### common functions to do calculations
transform_pvalue<-function(p){pmin(p,1-p)*2}


### df_calculations: results will be a list with cover_rate and p_values

mean_feo<-function(x) {
	x<-as.vector(x)
	x<-x[is.finite(x)]
	mean(x)
}

df_calculations<-function(x){
rep<-get("i",2)
iref<-get("i",2)
prop0p05_cover_rate<-NULL
pvaluesvsunif_cover_rate<-NULL
transformedpvaluesvsunif_cover_rate<-NULL
for (i in 1:iref)
{
#print(i)
if (sum(is.finite(x[i,]))>5)
	{
	temp<-mean(x[i,]<0.025|x[i,]>0.975)
	prop0p05_cover_rate<-c(prop0p05_cover_rate,ifelse(class(temp)=="try-error",NA,temp))
	pvaluesvsunif_cover_rate<-c(ks.test(x[i,],"punif")$p.value,pvaluesvsunif_cover_rate)
	transformedpvaluesvsunif_cover_rate<-c(ks.test(transform_pvalue(x[i,]),"punif")$p.value,transformedpvaluesvsunif_cover_rate)
	}
}
list(cover_rate=mean(prop0p05_cover_rate),
		p_value_cover_rate={p<-pbinom(mean(prop0p05_cover_rate)*N*iref,size=N*iref,prob=0.05); min(p,1-p)*2},
		p_value_test_uniformity=ks.test(as.vector(x),"punif")$p.value,
		p_value_test_uniformity_twice=ks.test(pvaluesvsunif_cover_rate,"punif")$p.value,
		p_value_test_uniformity_twice_transformed_p=ks.test(transform_pvalue(pvaluesvsunif_cover_rate),"punif")$p.value,
		p_value_test_uniformity_twice_transformedfirst_p=ks.test(transformedpvaluesvsunif_cover_rate,"punif")$p.value)
}

add_colours<-function(x,sheet,active_df) {
### defining columns to be analyzed
selected_columns<-which(regexpr("p_value",names(active_df))==1)
## changing background colour for very very significant results
selected_cells<-which(active_df[,selected_columns]<pvalue_very_very_signif,arr.ind=TRUE)
selected_cells<-cbind(selected_cells[,1],selected_columns[selected_cells[,2]])
if (length(selected_cells)>0)
{
	if (!is.na(selected_cells[1,1]))
	{
		for (i in 1: dim(selected_cells))
		{
		addStyle(x, sheet = sheet, style_very_very_signif, rows = selected_cells[i,1]+1, cols = selected_cells[i,2], gridExpand = TRUE)
		}
	}
}

## changing background colour for very significant results
selected_cells<-which(active_df[,selected_columns]<pvalue_very_signif&active_df[,selected_columns]>=pvalue_very_very_signif,arr.ind=TRUE)
selected_cells<-cbind(selected_cells[,1],selected_columns[selected_cells[,2]])
if (length(selected_cells)>0)
{
	if (!is.na(selected_cells[1,1]))
	{
		for (i in 1: dim(selected_cells))
		{
		addStyle(x, sheet = sheet, style_very_signif, rows = selected_cells[i,1]+1, cols = selected_cells[i,2], gridExpand = TRUE)
		}
	}
}


## changing background colour for significant results
selected_cells<-which(active_df[,selected_columns]<pvalue_signif&active_df[,selected_columns]>=pvalue_very_signif,arr.ind=TRUE)
selected_cells<-cbind(selected_cells[,1],selected_columns[selected_cells[,2]])
if (length(selected_cells)>0)
{
	if (!is.na(selected_cells[1,1]))
	{
		for (i in 1: dim(selected_cells))
		{
		addStyle(x, sheet = sheet, style_signif, rows = selected_cells[i,1]+1, cols = selected_cells[i,2], gridExpand = TRUE)
		}
	}
}

x
}
##############################################################################"

###############################################################################"
### common parameters
style_very_very_signif<-createStyle(fgFill="red")
style_very_signif<-createStyle(fgFill="orange")
style_signif<-createStyle(fgFill="yellow")

pvalue_very_very_signif<-0.001
pvalue_very_signif<-0.01
pvalue_signif<-0.05
########################################################################################""


### determining the files to analyze
names.files<-list.files()
names.files.to.analyze<-names.files[regexpr(".RData",names.files)>0&regexpr("_V5",names.files)>0&regexpr("Results_simplifiedanalysiscoverrate",names.files)==1&regexpr("_onlyTheort1",names.files)<0]

########################################################################################""
### going through the files to analyze:
for (nf in names.files.to.analyze)
{	print(paste0("Treating: ",nf))

	### attaching file nf in 2nd position:
	attach(nf,pos=2)
	
	### specifyning the parameters of the file:
	le_innf<-length(gregexpr("_",nf)[[1]])
	regressionContext<-substring(nf,first=gregexpr("_",nf)[[1]][le_innf-1]+1,last=gregexpr("_",nf)[[1]][le_innf]-1)
	Xresampled.local<-NA
	if (exists("Xresampled",2)) {Xresampled.local<-get("Xresampled",2)}
	start_result0<-c(dataSize=get("dataSize",2),regressionContext=regressionContext,rep=get("i",2),N=get("N",2),Nchains=get("Nchains",2),
					RhatMax=get("Rhat.max",2),neffMin=get("neff.min",2),Xresampled=Xresampled.local)
					
	regressionParameters<-{if (length(dim(data_cover_rate_bestmodel_semult))==4){c(1,2)} else {c(1)}}

	#### going through the parameters
	for (p in regressionParameters)
	{
		if (p==1) {regressionParameter<-"Intercept"}
		if (p==2) {regressionParameter<-"Slope"}
		start_result<-c(start_result0,regressionParameter=regressionParameter)
		
		### treating data_cover_rate_gauss_semult
			sheetname="Gaussian_semult"
			treateddata=data_cover_rate_gauss_semult
			active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
			#### specifying working data
			if (regressionContext=="Interceptmodel")
				{treateddata<-treateddata} else {treateddata<-treateddata[,,p]}
			#### making calculations on working data
			res<-df_calculations(treateddata)
			#### putting data in the active_df
			active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)))[names(active_df)])
			writeData(dfmethod_workbook,sheet=sheetname,active_df)
			dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
			print(paste0("end ",sheetname))
			
		
		### treating data_cover_rate_theort1_semult
			sheetname="Theort1_semult"
			treateddata=data_cover_rate_theort1_semult
			active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
			#### specifying working data
			if (regressionContext=="Interceptmodel")
				{treateddata<-treateddata} else {treateddata<-treateddata[,,p]}
			#### making calculations on working data
			res<-df_calculations(treateddata)
			#### putting data in the active_df
			active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)))[names(active_df)])
			writeData(dfmethod_workbook,sheet=sheetname,active_df)
			dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
			print(paste0("end ",sheetname))
		
		### treating data_cover_rate_theort2_semult
			sheetname="Theort2_semult"
			treateddata=data_cover_rate_theort2_semult
			active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
			#### specifying working data
			if (regressionContext=="Interceptmodel")
				{treateddata<-treateddata} else {treateddata<-treateddata[,,p]}
			#### making calculations on working data
			res<-df_calculations(treateddata)
			#### putting data in the active_df
			active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)))[names(active_df)])
			writeData(dfmethod_workbook,sheet=sheetname,active_df)
			dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
			print(paste0("end ",sheetname))
		
		### treating data_cover_rate_theort1_sesdt
			sheetname="Theort1_sesdt"
			treateddata=data_cover_rate_theort2_sesdt
			active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
			#### specifying working data
			if (regressionContext=="Interceptmodel")
				{treateddata<-treateddata} else {treateddata<-treateddata[,,p]}
			#### making calculations on working data
			res<-df_calculations(treateddata)
			#### putting data in the active_df
			active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)))[names(active_df)])
			writeData(dfmethod_workbook,sheet=sheetname,active_df)
			dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
			print(paste0("end ",sheetname))
		
		### treating data_cover_rate_theort2_sesdt
			sheetname="Theort2_sesdt"
			treateddata=data_cover_rate_theort2_sesdt
			active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
			#### specifying working data
			if (regressionContext=="Interceptmodel")
				{treateddata<-treateddata} else {treateddata<-treateddata[,,p]}
			#### making calculations on working data
			res<-df_calculations(treateddata)
			#### putting data in the active_df
			active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)))[names(active_df)])
			writeData(dfmethod_workbook,sheet=sheetname,active_df)
			dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
			print(paste0("end ",sheetname))
		
		
		### treating data_cover_rate_bestmodel_semult
			for (m in  1:dim(data_cover_rate_bestmodel_semult)[3])
				{sheetname=paste0("Bestmodel_semult_method",m)
				treateddata=data_cover_rate_bestmodel_semult
				active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
				#### specifying working data
				if (regressionContext=="Interceptmodel")
					{treateddata<-treateddata[,,m]} else {treateddata<-treateddata[,,m,p]}
				#### making calculations on working data
				res<-df_calculations(treateddata)
				converged<-NULL
				neffsreached<-NULL
				#### getting information about convergence and neffmin reached
				for (ii in 1:dim(data_cover_rate_bestmodel_semult)[1])
					{if (length(dim(data_cover_rate_bestmodel_semult))==3)
						{
						if (sum(is.finite(data_cover_rate_bestmodel_semult[ii,,]))>0) {
							if (mean_feo(data_cover_rate_bestmodel_semult[ii,,]==data_cover_rate_estimatedt_semult[ii,,])>0.5)
								{converged<-c(converged,res_Bayest[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_Bayest[[ii]]$final.params$neffs.reached)
								} else {
								converged<-c(converged,res_BayesGaus[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_BayesGaus[[ii]]$final.params$neffs.reached)
								}
							}
						} else {
						
						if (sum(is.finite(data_cover_rate_bestmodel_semult[ii,,,]))>0) {
							if (mean_feo(data_cover_rate_bestmodel_semult[ii,,,]==data_cover_rate_estimatedt_semult[ii,,,])>0.5)
								{if (regressionParameter=="Intercept")
									{
										converged<-c(converged,res_Bayest_Intercept[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_Bayest_Intercept[[ii]]$final.params$neffs.reached)
									} else {
										converged<-c(converged,res_Bayest_Slope[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_Bayest_Slope[[ii]]$final.params$neffs.reached)
									} 
								} else {
								if (regressionParameter=="Intercept")
									{
										converged<-c(converged,res_BayesGaus_Intercept[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_BayesGaus_Intercept[[ii]]$final.params$neffs.reached)
									} else {
										converged<-c(converged,res_BayesGaus_Slope[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_BayesGaus_Slope[[ii]]$final.params$neffs.reached)
									} 
								}
							}
						}
					}
				convneff<-c(percent_non_convergent=1-mean(converged),percent_neffmin_unreached=1-mean(neffsreached))
				#### putting data in the active_df
				active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)),t(convneff))[names(active_df)])
				writeData(dfmethod_workbook,sheet=sheetname,active_df)
				dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
				print(paste0("end ",sheetname))
				}
		
		### treating data_cover_rate_bestmodel_sesdt
			for (m in  1:dim(data_cover_rate_bestmodel_sesdt)[3])
				{sheetname=paste0("Bestmodel_sesdt_method",m)
				treateddata=data_cover_rate_bestmodel_sesdt
				active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
				#### specifying working data
				if (regressionContext=="Interceptmodel")
					{treateddata<-treateddata[,,m]} else {treateddata<-treateddata[,,m,p]}
				#### making calculations on working data
				res<-df_calculations(treateddata)
				converged<-NULL
				neffsreached<-NULL
				#### getting information about convergence and neffmin reached
				for (ii in 1:dim(data_cover_rate_bestmodel_sesdt)[1])
					{if (length(dim(data_cover_rate_bestmodel_sesdt))==3)
						{
						if (sum(is.finite(data_cover_rate_bestmodel_sesdt[ii,,]))>0) {
							if (mean_feo(data_cover_rate_bestmodel_sesdt[ii,,]==data_cover_rate_estimatedt_sesdt[ii,,])>0.5)
								{converged<-c(converged,res_Bayest[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_Bayest[[ii]]$final.params$neffs.reached)
								} else {
								converged<-c(converged,res_BayesGaus[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_BayesGaus[[ii]]$final.params$neffs.reached)
								}
							}
						} else {
						
						if (sum(is.finite(data_cover_rate_bestmodel_sesdt[ii,,,]))>0) {
							if (mean_feo(data_cover_rate_bestmodel_sesdt[ii,,,]==data_cover_rate_estimatedt_sesdt[ii,,,])>0.5)
								{if (regressionParameter=="Intercept")
									{
										converged<-c(converged,res_Bayest_Intercept[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_Bayest_Intercept[[ii]]$final.params$neffs.reached)
									} else {
										converged<-c(converged,res_Bayest_Slope[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_Bayest_Slope[[ii]]$final.params$neffs.reached)
									} 
								} else {
								if (regressionParameter=="Intercept")
									{
										converged<-c(converged,res_BayesGaus_Intercept[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_BayesGaus_Intercept[[ii]]$final.params$neffs.reached)
									} else {
										converged<-c(converged,res_BayesGaus_Slope[[ii]]$final.params$converged)
										neffsreached<-c(neffsreached,res_BayesGaus_Slope[[ii]]$final.params$neffs.reached)
									} 
								}
							}
						}
					}
				convneff<-c(percent_non_convergent=1-mean(converged),percent_neffmin_unreached=1-mean(neffsreached))
				#### putting data in the active_df
				active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)),t(convneff))[names(active_df)])
				writeData(dfmethod_workbook,sheet=sheetname,active_df)
				dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
				print(paste0("end ",sheetname))
				}
		
		
		### treating data_cover_rate_estimatedt_semult
			for (m in  1:dim(data_cover_rate_estimatedt_semult)[3])
				{sheetname=paste0("Estimatedt_semult_method",m)
				treateddata=data_cover_rate_estimatedt_semult
				active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
				#### specifying working data
				if (regressionContext=="Interceptmodel")
					{treateddata<-treateddata[,,m]} else {treateddata<-treateddata[,,m,p]}
				#### making calculations on working data
				res<-df_calculations(treateddata)
				converged<-NULL
				neffsreached<-NULL
				#### getting information about convergence and neffmin reached
				for (ii in 1:dim(data_cover_rate_estimatedt_semult)[1])
					{
						if (length(dim(data_cover_rate_estimatedt_semult))==3)
						{
							if (sum(is.finite(data_cover_rate_estimatedt_semult[ii,,]))>0) {
								converged<-c(converged,res_Bayest[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_Bayest[[ii]]$final.params$neffs.reached)
							}
						} else {
							if (regressionParameter=="Intercept")
									{
										if (sum(is.finite(data_cover_rate_estimatedt_semult[ii,,,]))>0) {
											converged<-c(converged,res_Bayest_Intercept[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayest_Intercept[[ii]]$final.params$neffs.reached)
										}
									} else {
										if (sum(is.finite(data_cover_rate_estimatedt_semult[ii,,,]))>0) {
											converged<-c(converged,res_Bayest_Slope[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayest_Slope[[ii]]$final.params$neffs.reached)
											}
									} 
						}
					}
				convneff<-c(percent_non_convergent=1-mean(converged),percent_neffmin_unreached=1-mean(neffsreached))
				#### putting data in the active_df
				active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)),t(convneff))[names(active_df)])
				writeData(dfmethod_workbook,sheet=sheetname,active_df)
				dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
				print(paste0("end ",sheetname))
				}
		
		### treating data_cover_rate_estimatedt_sesdt
			for (m in  1:dim(data_cover_rate_estimatedt_sesdt)[3])
				{sheetname=paste0("Estimatedt_sesdt_method",m)
				treateddata=data_cover_rate_estimatedt_sesdt
				active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
				#### specifying working data
				if (regressionContext=="Interceptmodel")
					{treateddata<-treateddata[,,m]} else {treateddata<-treateddata[,,m,p]}
				#### making calculations on working data
				res<-df_calculations(treateddata)
				converged<-NULL
				neffsreached<-NULL
				#### getting information about convergence and neffmin reached
				for (ii in 1:dim(data_cover_rate_estimatedt_sesdt)[1])
					{
						if (length(dim(data_cover_rate_estimatedt_sesdt))==3)
						{
							if (sum(is.finite(data_cover_rate_estimatedt_sesdt[ii,,]))>0) {
								converged<-c(converged,res_Bayest[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_Bayest[[ii]]$final.params$neffs.reached)
							}
						} else {
							if (regressionParameter=="Intercept")
									{
										if (sum(is.finite(data_cover_rate_estimatedt_sesdt[ii,,,]))>0) {
											converged<-c(converged,res_Bayest_Intercept[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayest_Intercept[[ii]]$final.params$neffs.reached)
										}
									} else {
										if (sum(is.finite(data_cover_rate_estimatedt_sesdt[ii,,,]))>0) {
											converged<-c(converged,res_Bayest_Slope[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayest_Slope[[ii]]$final.params$neffs.reached)
											}
									} 
						}
					}
				convneff<-c(percent_non_convergent=1-mean(converged),percent_neffmin_unreached=1-mean(neffsreached))
				#### putting data in the active_df
				active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)),t(convneff))[names(active_df)])
				writeData(dfmethod_workbook,sheet=sheetname,active_df)
				dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
				print(paste0("end ",sheetname))
				}
				
				
		
		### treating data_cover_rate_estimatedtbis_semult
			for (m in  1:dim(data_cover_rate_estimatedtbis_semult)[3])
				{sheetname=paste0("Estimatedtbis_semult_method",m)
				treateddata=data_cover_rate_estimatedtbis_semult
				active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
				#### specifying working data
				if (regressionContext=="Interceptmodel")
					{treateddata<-treateddata[,,m]} else {treateddata<-treateddata[,,m,p]}
				#### making calculations on working data
				res<-df_calculations(treateddata)
				converged<-NULL
				neffsreached<-NULL
				#### getting information about convergence and neffmin reached
				for (ii in 1:dim(data_cover_rate_estimatedtbis_semult)[1])
					{
						if (length(dim(data_cover_rate_estimatedtbis_semult))==3)
						{
							if (sum(is.finite(data_cover_rate_estimatedtbis_semult[ii,,]))>0) {
								converged<-c(converged,res_Bayestbis[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_Bayestbis[[ii]]$final.params$neffs.reached)
							}
						} else {
							if (regressionParameter=="Intercept")
									{
										if (sum(is.finite(data_cover_rate_estimatedtbis_semult[ii,,,]))>0) {
											converged<-c(converged,res_Bayestbis_Intercept[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayestbis_Intercept[[ii]]$final.params$neffs.reached)
										}
									} else {
										if (sum(is.finite(data_cover_rate_estimatedtbis_semult[ii,,,]))>0) {
											converged<-c(converged,res_Bayestbis_Slope[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayestbis_Slope[[ii]]$final.params$neffs.reached)
											}
									} 
						}
					}
				convneff<-c(percent_non_convergent=1-mean(converged),percent_neffmin_unreached=1-mean(neffsreached))
				#### putting data in the active_df
				active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)),t(convneff))[names(active_df)])
				writeData(dfmethod_workbook,sheet=sheetname,active_df)
				dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
				print(paste0("end ",sheetname))
				}
		
		### treating data_cover_rate_estimatedtbis_sesdt
			for (m in  1:dim(data_cover_rate_estimatedtbis_sesdt)[3])
				{sheetname=paste0("Estimatedtbis_sesdt_method",m)
				treateddata=data_cover_rate_estimatedtbis_sesdt
				active_df<-read.xlsx(dfmethod_workbook,sheet=sheetname,skipEmptyRows = FALSE, skipEmptyCols = FALSE)
				#### specifying working data
				if (regressionContext=="Interceptmodel")
					{treateddata<-treateddata[,,m]} else {treateddata<-treateddata[,,m,p]}
				#### making calculations on working data
				res<-df_calculations(treateddata)
				converged<-NULL
				neffsreached<-NULL
				#### getting information about convergence and neffmin reached
				for (ii in 1:dim(data_cover_rate_estimatedtbis_sesdt)[1])
					{
						if (length(dim(data_cover_rate_estimatedtbis_sesdt))==3)
						{
							if (sum(is.finite(data_cover_rate_estimatedtbis_sesdt[ii,,]))>0) {
								converged<-c(converged,res_Bayestbis[[ii]]$final.params$converged)
								neffsreached<-c(neffsreached,res_Bayestbis[[ii]]$final.params$neffs.reached)
							}
						} else {
							if (regressionParameter=="Intercept")
									{
										if (sum(is.finite(data_cover_rate_estimatedtbis_sesdt[ii,,,]))>0) {
											converged<-c(converged,res_Bayestbis_Intercept[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayestbis_Intercept[[ii]]$final.params$neffs.reached)
										}
									} else {
										if (sum(is.finite(data_cover_rate_estimatedtbis_sesdt[ii,,,]))>0) {
											converged<-c(converged,res_Bayestbis_Slope[[ii]]$final.params$converged)
											neffsreached<-c(neffsreached,res_Bayestbis_Slope[[ii]]$final.params$neffs.reached)
											}
									} 
						}
					}
				convneff<-c(percent_non_convergent=1-mean(converged),percent_neffmin_unreached=1-mean(neffsreached))
				#### putting data in the active_df
				active_df<-rbind.data.frame(active_df,cbind.data.frame(t(start_result),t(unlist(res)),t(convneff))[names(active_df)])
				writeData(dfmethod_workbook,sheet=sheetname,active_df)
				dfmethod_workbook<-add_colours(dfmethod_workbook,sheet=sheetname,active_df)
				print(paste0("end ",sheetname))
				}
				
				
	}
	### detaching file in 2nd position: (nf)
	detach(pos=2)
}





### saving xlsx workbook:
setwd(WD.results)
saveWorkbook(dfmethod_workbook,file=dfmethod_xlsfile_final,overwrite=TRUE)
