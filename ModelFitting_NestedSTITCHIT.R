args <- commandArgs(TRUE)
require("methods")
require("glmnet")
require("doMC")
require("parallel")
if(length(args) < 1) {
	args <- c("--help")
}
 
## Help 
if("--help" %in% args) {
	cat("
	INVOKE offers linear regression with Lasso, Ridge, and Elastic Net regularisation.
	Arguments:
	--outDir Output directory (will be created if it does not exist)
	--dataDir Directory containing the data
	--response Name of the response variable
	--cores Number of cores to be use (default 1)
	--fixedAlpha Use a fixed value for the alpha parameter in elastic net regulatisation, do not perform a grid search
	--alpha Stepsize to optimise the alpha parameter in elastic net regularisation (default 0.05)
	--testsize Size of test data[%] (default 0.2)
	--regularisation L for Lasso, R for Ridge, and E for Elastic net (default E)
	--innerCV Number of folds for inner cross-validation (default 6)
	--constraint Specifies a constraint on the coefficent sign, enter N for negative and P for positive constraint
	--seed Random seed used for random number generation (default random)
	--asRData Store feature coefficients as RData files (default FALSE)
	--logResponse Flag indicating whether the response variable should be log transformed (default TRUE)
	--coefP p-value threshold for model coefficient (default 1, all OLS coefs will be returned)
	--help=print this text
")
	q(save="no")
}

# Process command arguments
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL$outDir)) {
	cat("No output directory specified. Use the --outDir option to specify an output directory.")
	q(save="no")
}
argsL$outDir<-paste0(argsL$outDir,"/")

if(is.null(argsL$dataDir)) {
	cat("No data directory specified. Use the --dataDir option to specify a data directory.")
	q(save="no")
}
argsL$dataDir<-paste0(argsL$dataDir,"/")
Data_Directory <- argsL$dataDir

if(is.null(argsL$response)) {
	cat("No response variable name specified. Use the --response option to specify a response variable.")
	q(save="no")
}

if(is.null(argsL$testsize)){
	argsL$testsize <- 0.2
}

if(is.null(argsL$innerCV)){
	argsL$innerCV <- 6
}

if(is.null(argsL$alpha)) {
	argsL$alpha <- 0.05
}

if(is.null(argsL$cores)) {
	argsL$cores <- 1
}

if(is.null(argsL$coefP)) {
	argsL$coefP <- 1
}

if(is.null(argsL$regularisation)){
	argsL$regularisation<-c("E")
}

if(is.null(argsL$constraint)){
	lower_bound <- NULL
	upper_bound <- NULL
}else if(argsL$constraint=="P"){
	lower_bound <- 0
}else if(argsL$constraint=="N"){
	upper_bound <- 0
}

if(is.null(argsL$fixedAlpha)){
	argsL$fixedAlpha <- -1
}

if (is.null(argsL$performance)){
	argsL$performance <- TRUE
}

if (! is.null(argsL$seed)){
	set.seed(as.numeric(argsL$seed))
}


if (is.null(argsL$asRData)){
	argsL$asRData <- FALSE
}

if (is.null(argsL$logResponse)){
	argsL$logResponse <- TRUE
}


registerDoMC(cores = argsL$cores)

permute<-function(x,resPos){
s<-sample(length(x))
s<-s[which(s != resPos)]
c(x[s],x[resPos])
}

#Check output directory, create it if necessary
dir.create(argsL$outDir,showWarning=FALSE)

#Initilaise lists for storage of intermediate results
FileList<-list.files(path=Data_Directory)
numFiles=length(FileList)
pearson_correlation<-vector("list",numFiles)
spearman_correlation<-vector("list",numFiles)
test_error<-vector("list",numFiles)
rss_error<-vector("list",numFiles)
ftest_result<-vector("list",numFiles)
coefficients<-vector("list",numFiles)
coefficientsF<-vector("list",numFiles)
Sample_View<-vector("list",numFiles)
validSamples<-vector("logical",numFiles)
spearmanPassed<-vector("numeric",numFiles)

#Print sample names
if (length(FileList)==0){
	exit()
}
counter<-0
for(Sample in FileList){
	counter<-counter+1
	}

#Loop through sample files
i<-0
for(Sample in FileList){
	i<-i+1
	#Loading and preprocessing data
	M<-read.table(paste(Data_Directory,Sample,sep=""),header=TRUE,sep="",row.names=1)
	
	M<-unique(M)
	M<-data.frame(M)
	FeatureNames_temp<-colnames(M)
	Response_Variable_location_temp <- grep(argsL$response,FeatureNames_temp)

	vectorLength<-nrow(M)    
	pearson_correlation[[i]]<-cbind(vector("list",vectorLength),vector("list",vectorLength))
	spearman_correlation[[i]]<-cbind(vector("list",vectorLength),vector("list",vectorLength))
	test_error[[i]]<-vector("list",vectorLength)
	coefficients[[i]]<-vector("list",vectorLength)

	if (min(M[,Response_Variable_location_temp]) >= 0){
		if (argsL$logResponse == TRUE){
			M<-log2(M+1)
		}
	}else{
		Response_Variable_location_temp <- grep(argsL$response,FeatureNames_temp)
		M[,-Response_Variable_location_temp]<-log2(M[,-Response_Variable_location_temp]+1)
	}
	SD<-apply(M,2,sd)
	Feature_zero_SD<-as.vector(which(SD==0))
	if(length(Feature_zero_SD)>0){
		if (Response_Variable_location_temp %in% Feature_zero_SD){
			validSamples[i]=FALSE
			next;
			}
		M<-M[,-c(Feature_zero_SD)]
	}
	if (is.null(dim(M))){
		validSamples[i]=FALSE
		spearmanPassed[i]=1
		next;
	}
	if (dim(M)[2] < 2){
		validSamples[i]=FALSE
		spearmanPassed[i]=1
		next;
	}
	if (length(which(M==0))>(dim(M)[1]*dim(M)[2]*0.5)){
		validSamples[i]=FALSE
		spearmanPassed[i]=1
		next;
	}

	FeatureNames<-colnames(M)
	M<-data.frame(scale(M,center=TRUE, scale=TRUE))
     if (dim(M)[1] < 30){
          validSamples[i]=FALSE
		spearmanPassed[i]=1
          next;
     }else{
          validSamples[i]=TRUE
     }

	name<-unlist(unlist(strsplit(Sample, ".txt")))
	Response_Variable_location<- grep(argsL$response,FeatureNames)

	predictedAll<-c()
	measuredAll<-c()
	Train_Data<-M

	# Split the features from response
	x_train<-as.matrix(Train_Data[,-Response_Variable_location])
	y_train<-as.vector(unlist(Train_Data[,Response_Variable_location,drop=FALSE]))

	#Creating alpha vector
 	A<-c()
	if(argsL$regularisation=="L"){
		alphaslist <- c(1.0)
	}else{
		if(argsL$regularisation=="R"){
		alphaslist <- c(0.0)
		}else{
			alphaslist<-seq(0,1,by=as.numeric(argsL$alpha))
		}
	}
	#Learning model on training data
	if(argsL$regularisation=="E"){   
		if(argsL$fixedAlpha==-1){
			if(is.null(argsL$constraint)){
				elasticnet<-mclapply(alphaslist, function(x){;cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
			}else{ 
				if(argsL$constraint=="P"){
					elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
					}else{
						if(argsL$constraint=="N"){
						elasticnet<-mclapply(alphaslist, function(x){cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV))}, mc.cores=argsL$cores)
						}
					}
	      		}
		}else{
		x=argsL$fixedAlpha
		if(is.null(argsL$constraint)){
			elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
		}else{ 
			if(argsL$constraint=="P"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
			}else{
				if(argsL$constraint=="N"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
						}
				}
      	}
		}   
	}else{
	x=alphaslist[1]

	if(is.null(argsL$constraint)){
			elasticnet<-cv.glmnet(x_train, y_train,alpha=x,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
			}else{ 
			if(argsL$constraint=="P"){
				elasticnet<-cv.glmnet(x_train, y_train,alpha=x,lower=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
			}else{
				if(argsL$constraint=="N"){
					elasticnet<-cv.glmnet(x_train, y_train,alpha=x,upper=0,nfolds=as.numeric(argsL$innerCV),parallel=TRUE)
					}
				}
    			}
		}
	if(length(elasticnet[[1]]) > 1){
		if (argsL$regularisation=="E"){
				if(argsL$fixedAlpha==-1){
					for (j in 1:length(alphaslist)) {
						A[j]<-min(elasticnet[[j]]$cvm)
					}
					#Determine best alpha value from training data
					index<-which(A==min(A), arr.ind=TRUE)
					model<-elasticnet[[index]]
			}
		}else{
		model<-elasticnet
			}
	saveRDS(model,file=paste0(argsL$outDir,"ElasticNet_model_",unlist(unlist(strsplit(Sample,".txt"))),".RDS"))
	write.table(rownames(model[[8]][2]$beta),file=paste0(argsL$outDir,"Selected_Regions_ElasticNet_",unlist(unlist(strsplit(Sample,".txt"))),".bed"),quote=F,row.names=T,col.names=F,sep="\t")
	}

	k=1
	if (length(elasticnet[[1]]) > 1){ 
		#Determine error of the best alpha model on hold out data and on training data
		predict_fit_train<-predict(model, x_train, s="lambda.min")
		coefficients[[i]][[k]]<-coef(model, s = "lambda.min")
		pearson_correlation[[i]][k]<-cor(predict_fit_train,y_train)
		spearman_correlation[[i]][k]<-cor(predict_fit_train,y_train,method='spearman')
		predictedAll<-c(predictedAll,predict_fit_train)
		measuredAll<-c(measuredAll,y_train)
		test_error[[i]][k]<-sum((y_train-predict_fit_train)^2)/length(y_train)
		rss_error[k]<-sum((y_train-predict_fit_train)^2)
	}

	#Learning the model once on the full data set
	if (! is.null(argsL$seed)){
		set.seed(as.numeric(argsL$seed))
		}
	#Determine nonzero model coefficients
	modelCoefMatrix<-c()
	for (j in 1:length(coefficients[[i]])){
		if (length(coefficients[[i]][[j]]>1)){
			modelCoefMatrix<-rbind(modelCoefMatrix,coefficients[[i]][[j]][,1])
			}
		}
	if (length(modelCoefMatrix) != 0){
	medianModelCoefMatrix<-apply(modelCoefMatrix,2,median)[-1]
	nObs<-dim(M)[1]
	# Partition data into test and training data sets
	if (length(which(medianModelCoefMatrix!=0))){
		if (length(which(medianModelCoefMatrix!=0))>=nObs){
			ols_Data<-M[,c(order(abs(medianModelCoefMatrix),decreasing=T)[1:(nObs-2)],Response_Variable_location)]
		}else{
			ols_Data<-M[,c(which(medianModelCoefMatrix!=0),Response_Variable_location)]
			} 
		model<-lm(Expression~.,ols_Data)	
		model.coefs<-summary(model)$coefficients[,c(1,4)]
		signif.coefs<-which(model.coefs[,2]<=as.numeric(argsL$coefP))
		model.coefs.signif<-model.coefs[signif.coefs,]
		if (length(signif.coefs > 0)){
			for (j in  1:length(row.names(model.coefs.signif))){
				row.names(model.coefs.signif)[j]<-gsub(".","\t",row.names(model.coefs.signif)[j],fixed=T)
			}
			saveRDS(model,file=paste0(argsL$outDir,"OLS_model_",unlist(unlist(strsplit(Sample,".txt"))),".RDS"))
			if (length(signif.coefs)>1){
				if (row.names(model.coefs.signif)[1]=="(Intercept)"){
					write.table(model.coefs.signif[-1,],file=paste0(argsL$outDir,"Selected_Regions_",unlist(unlist(strsplit(Sample,".txt"))),".bed"),quote=F,row.names=T,col.names=F,sep="\t")
				}else{
				write.table(model.coefs.signif,file=paste0(argsL$outDir,"Selected_Regions_",unlist(unlist(strsplit(Sample,".txt"))),".bed"),quote=F,row.names=T,col.names=F,sep="\t")
				}
			}else{
			cat(paste0(gsub(".","\t",row.names(model.coefs)[signif.coefs],fixed=T),"\t",model.coefs.signif[1],"\t",model.coefs.signif[2],"\n"),file=paste0(argsL$outDir,"Selected_Regions_",unlist(unlist(strsplit(Sample,".txt"))),".bed"))
				}
			}
		}
	}
}
	
###############################
###Writing model performance### #Change to print out
###############################
if (argsL$performance == TRUE){
	for (i in 1:length(FileList)){
		if (validSamples[i]==FALSE){
     	          next;
	     }
	     cm<-mean(unlist(pearson_correlation[[i]]),na.rm=TRUE)
	     csd<-var(unlist(pearson_correlation[[i]]),na.rm=TRUE)
	     cms<-mean(unlist(spearman_correlation[[i]]),na.rm=TRUE)
     	csds<-var(unlist(spearman_correlation[[i]]),na.rm=TRUE)
	     erm<-mean(unlist(test_error[[i]]),na.rm=TRUE)
	     ersd<-var(unlist(test_error[[i]]),na.rm=TRUE)
	     Sample_View[[i]]<-data.frame(Sample_Name=FileList[i],Pearson=cm,Searman=cms,MSE=erm)
	}		
	Sample_ViewF<-do.call("rbind",Sample_View)
	cat(paste(gsub("Segmentation_","",gsub("_Spearman_10.txt","", as.character(Sample_ViewF[,1]))),Sample_ViewF[,2],Sample_ViewF[,3],Sample_ViewF[,4],"\n",sep=" "))
}
