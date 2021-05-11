args <- commandArgs(TRUE)
library("methods")
library("glmnet")
if(length(args) < 1) {
	args <- c("--help")
}
 
## Help 
if("--help" %in% args) {
	cat("
	INVOKE offers linear regression with Lasso, Ridge, and Elastic Net regularisation.
	Arguments:
	--outDir Output directory (will be created if it does not exist)
	--testdata Path to the file containing the test data
	--response Name of the response variable
	--targetGeneID ID of the gene of interest
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

if(is.null(argsL$testdata)) {
	cat("No data directory specified. Use the --testData option to specify a data directory.")
	q(save="no")
}
Data_Path <- argsL$testdata

if(is.null(argsL$targetGeneID)) {
	cat("No target geneID specified. Use the --targetGeneID option to specify a target gene ID.")
	q(save="no")
}
geneID <- argsL$targetGeneID


if(is.null(argsL$response)) {
	cat("No response variable name specified. Use the --response option to specify a response variable.")
	q(save="no")
}
response <- argsL$response


#Check output directory, create it if necessary
dir.create(argsL$outDir,showWarning=FALSE)

#Read the model
model <- readRDS(paste0(Data_Path,"_ROut/OLS_model_Segmentation_",geneID,"_Spearman_10.RDS",sep=""))

#Read the test data
tmp <- read.table(paste0(Data_Path,"/Segmentation_",geneID,"_TestData.txt"),header=T)
M <- log2(tmp+1)
M <- scale(M,center=T,scale=T)
test.data <- M[,c(1:dim(M)[2]-1)]
expression<- M[,dim(M)[2]]

result<-predict.lm(model,data.frame(test.data))
cat(paste(geneID,cor(result,expression,method="pearson"),cor(result,expression,method="spearman"),"\n",sep=" "))

