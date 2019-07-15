#### IMPORTS ############################################
library(impute)
library(missForest)
library(Biobase)
library(pcaMethods)
library(ggfortify)
library(imputeLCMD)
#########################################################

#' Metabolomic Dataset Imputation
#'
#' Impute the given dataset with different method options.
#' Produces <filename>_imputed.txt, containing the imputed dataset;
#' See Details for available Imputation Methods
#'
#' Available imputation methods:
#' \itemize{
#' \item "knn": From the \code{impute} package, use the k nearest neighboors to impute the values;
#' \item "RF": From the \code{missForest} package, use RandomForest algorithm to impute the values;
#' \item "QRILC": From the \code{imputeLCMD} package, use Quantile regression to impute the values;
#' \item "SVD": From the \code{pcaMethods} package, use SVDimpute algorithm as proposed by Troyanskaya et al, 2001. to impute the values;
#' \item "mean","median", ""median", "0", "HM": simple value replacement, either by the mean, median, 0 of Half minimum of the row;
#' }
#'
#' @param file file containing the dataset to impute;
#' \itemize{
#' \item a column "Compound";
#' \item a column "Metabolite";
#' \item and all the columns sample from <sampleStart> to the end of the file;
#' }
#'       the rest doesn't matter and the names are optional, as long as the column position is entered.
#' @param method default "knn", the chosen method for replacing the missing values. Can be "knn", "RF",
#'  "QRILC", "SVD", "mean", "median", "HM" or "0". See Details.
#'
#' @param k default 2, the k used for the knn imputation;
#' @param npcs default 3, npcs for SVD method;
#' @param sigma default 0.1, tune sigma parameter for QRILC method;
#' @param nTree default 30, number of tree for the RF method;
#' @param transformation default "None", can be "scale" or "log";
#' @param na.string default "NA", string to consider as NA in the dataset;
#' @param compound default NULL, position of the compound column if named otherwise;
#' @param metabolite default NULL, position of the metabolite column if named otherwise;
#' @param sampleStart default 3, 1st column of the actual data;
#'
#' @return df, the imputed dataset as a dataframe.
#'
#' @examples
#' for a dataset with the following header ; Compound, m/z, Metabolite, RT, Sample #1, ...
#' imputation("dummySet.tsv", method="knn", transformation="log", sampleStart=5)
#'
#' for a dataset with the following header ; compound, m/z, metabolite, RT, Sample #1, ...
#' imputation("dummySet.tsv", method="knn", transformation="log", metabolite = 3, compound = 1, sampleStart=5)
#'
#'
#' @seealso \itemize{
#' \item \code{impute} package\url{https://www.rdocumentation.org/packages/impute}
#' \item \code{missForest} package \url{https://www.rdocumentation.org/packages/missForest}
#' \item \code{imputeLCMD} package \url{https://www.rdocumentation.org/packages/imputeLCMD}
#' \item \code{pcaMethods} package \url{https://www.rdocumentation.org/packages/pcaMethods}
#' }
#'
#' @import impute
#' @import missForest
#' @import Biobase
#' @import pcaMethods
#' @import ggfortify
#' @import imputeLCMD
#'
#' @export
imputation <- function(file,
                       k=2, method="knn",
                       npcs=3,
                       sigma=0.1,
                       nTree=30,
                       na.string="NA",
                       transformation="None",
                       compound=NULL,
                       metabolite=NULL,
                       sampleStart=3){

  #### Specify Data files and parameters ###################################################
  inputFile <- file
  outputFile <- paste0(strsplit(file, ".", fixed=TRUE)[[1]][1],"_imputed.txt")
  knn_k = k

  myData1 <- read.table(inputFile,header=TRUE,sep="\t")

  if (!is.null(compound)){
    colnames(myData1)[compound] <- "Compound"
  }

  if (!is.null(metabolite)){
    colnames(myData1)[metabolite] <- "Metabolite"
  }
  if(is.null(myData1$Compound)){
    stop("no Compound column or column number entered")

  }
  if(is.null(myData1$Metabolite)){
    stop("no Metabolite column or column number entered")

  }

  myData1_noZero <- as.matrix(myData1[,sampleStart:NCOL(myData1)])

  myData1 <- cbind(Compound=myData1$Compound, Metabolite=myData1$Metabolite,
                   myData1[,sampleStart:NCOL(myData1)])

  for(i in 1:(dim(myData1)[2]-2)){
    if((class(myData1[,i+2]) == "factor")|(class(myData1[,i+2]) == "character")){
      warning(paste(paste0("Warning: sampleStart might be wrong: column #", i), " are factors or characters"))
    }
  }
  #### Transformation ####################################################################
  if(transformation=="log"){

    #Add 1 in each cell to overcome issue of 0
    myData1_noZero_temp <- myData1_noZero + 1

    #Transform intensity data to the log scale
    myDataMatrix <- log(myData1_noZero_temp)

    #Retrieve row log-mean values
    rowMean <- rowMeans(myDataMatrix,na.rm=TRUE)

    #Substract row log-mean from each value (recenter each row to 0)
    for ( i in 1:NROW(myDataMatrix)){
      myDataMatrix[i,] = myDataMatrix[i,] - rowMean[i]
    }

  }else if(transformation=="scale"){
    myDataMatrixT <- scale(t(myData1_noZero))
    myDataMatrix <- t(myDataMatrixT)
  }else{
    myDataMatrix <- as.matrix(myData1_noZero)
  }




  #### Imputation method #############################################################
  if (method=="knn"){
    myDataMatrixImputed <- impute.knn(myDataMatrix,
                                      k=knn_k,rowmax=0.8)$data

  }else if(method== "RF"){
    myDataMatrixImputed <- missForest(myDataMatrix, decreasing = T, ntree=nTree)$ximp

  }else if(method=="mean"){
    k <- which(is.na(myDataMatrix), arr.ind=TRUE)
    myDataMatrix[k] <- rowMeans(myDataMatrix, na.rm=TRUE)[k[,1]]
    myDataMatrixImputed <- myDataMatrix

  }else if(method=="median"){
    k <- which(is.na(myDataMatrix), arr.ind=TRUE)
    myDataMatrix[k] <- rowMedians(myDataMatrix, na.rm=TRUE)[k[,1]]
    myDataMatrixImputed <- myDataMatrix

  }else if(method=="0"){
    myDataMatrix[is.na(myDataMatrix)] <- 0
    myDataMatrixImputed <- myDataMatrix

  }else if(method=="HM"){
    k <- which(is.na(myDataMatrix), arr.ind=TRUE)
    myDataMatrix[k] <- min(myDataMatrix[k[,1],], na.rm=TRUE)/2
    myDataMatrixImputed <- myDataMatrix

  }else if(method=="SVD"){
    myDataMatrixImputed <- pca(t(myDataMatrix), method="svdImpute", nPcs=npcs, center = F)
    myDataMatrixImputed <- t(myDataMatrixImputed@completeObs)

  }else if(method=="QRILC"){
    if(transformation!="log"){
      warning("Log transformation is necessary for this method")
    }
    myDataMatrixImputed <- impute.QRILC(myDataMatrix, tune.sigma= sigma)[[1]]

  }else{
    stop(paste0("Not a valid method: ", method))
  }


  #### Transform back ################################################################
  if(transformation=="log"){

    for ( i in 1:NROW(myDataMatrixImputed)){
      myDataMatrixImputed[i,] = exp(myDataMatrixImputed[i,] + rowMean[i])
    }

  }else if(transformation=="scale"){
    myDataMatrixImputed <-( myDataMatrixImputed *attr(myDataMatrixT, 'scaled:scale'))+ attr(myDataMatrixT, 'scaled:center')
  }

  #### Round ##########################################################################
  myDataMatrixImputed <- round(myDataMatrixImputed)


  #### Write the imputed data to the output files #####################################
  print("saving imputated data")

  write.table(data.frame(Compound=myData1$Compound, Metabolite=myData1$Metabolite, myDataMatrixImputed),
              file = outputFile, row.names=FALSE,sep="\t")


  return(df)
}

#' Test for the various methods
#'
#' Test imputation accuracy for the given dataset with different method options. Compute the NRMSE values for the desired number of test.
#' If asked, produces an output <filename>_Accuracy.txt with the columns
#' \itemize{
#' \item "Method";
#' \item "missing_proportion";
#' \item "transformation";
#' \item  and	"NRMSE"
#' }
#'
#'
#' Will compute de NRMSE (Normalized Root Mean Squared Error) for an imputation test.
#' Available methods:
#' \itemize{
#' \item "knn": From the \code{impute} package, use the k nearest neighboors to impute the values;
#' \item "RF": From the \code{missForest} package, use RandomForest algorithm to impute the values;
#' \item "QRILC": From the \code{imputeLCMD} package, use Quantile regression to impute the values;
#' \item "SVD": From the \code{pcaMethods} package, use SVDimpute algorithm as proposed by Troyanskaya et al, 2001. to impute the values;
#' \item "mean","median", ""median", "0", "HM": simple value replacement, either by the mean, median, 0 of Half minimum of the row;
#' }
#'
#' @param input file containing the test dataset; Must contain:
#' \itemize{
#' \item a column "Compound";
#' \item a column "Metabolite";
#' \item and all the columns sample from <sampleStart> to the end of the file;
#' }
#' the rest doesn't matter and the names are optional, as long as the column position is entered.
#' *MUST BE COMPLETE (no NA values) for accurate results
#' @param output default NULL, name of the test results file if not NULL
#' @param method default "knn", the chosen method for replacing the missing values. Can be "knn", "RF",
#'  "QRILC", "SVD", "mean", "median", "HM" or "0". See Details.
#' @param nTest, default 10, number of test to loop;
#' @param k default 2, the k used for the knn imputation;
#' @param npcs default 3, npcs for SVD method;
#' @param sigma default 0.1, tune sigma parameter for QRILC method;
#' @param nTree default 30, number of tree for the RF method;
#' @param transformation default "None", can be "scale" or "log";
#' @param missing default 0.05, proportion of missing values;
#' @param missingType, default "MCAR" (missing completely at random), can be "MNAR" (not at ramdom), will target the values under the median;
#' @param na.string default "NA", string to consider as NA in the dataset;
#' @param compound default NULL, position of the compound column if named otherwise;
#' @param metabolite default NULL, position of the metabolite column if named otherwise;
#' @param sampleStart default 3, 1st column of the actual data;
#'
#' @return resDf, the result NRMSE dataframe.
#'
#' @examples
#' for a dataset with the following header ; Compound, m/z, Metabolite, RT, Sample #1, ...
#' imputationTest("dummySet.tsv", method="knn", transformation="log", sampleStart=5)
#'
#' for a dataset with the following header ; compound, m/z, metabolite, RT, Sample #1, ...
#' imputationTest("dummySet.tsv", method="knn", transformation="log", metabolite = 3, compound = 1, sampleStart=5)
#'
#'
#' @seealso \itemize{
#' \item \code{impute} package\url{https://www.rdocumentation.org/packages/impute}
#' \item \code{missForest} package \url{https://www.rdocumentation.org/packages/missForest}
#' \item \code{imputeLCMD} package \url{https://www.rdocumentation.org/packages/imputeLCMD}
#' \item \code{pcaMethods} package \url{https://www.rdocumentation.org/packages/pcaMethods}
#' }
#'
#' @import impute
#' @import missForest
#' @import Biobase
#' @import pcaMethods
#' @import ggfortify
#' @import imputeLCMD
#'
#' @export
imputationTest <- function(output=NULL, input,
                           k=2, method="knn",
                           npcs=3,
                           sigma=0.1,
                           nbTest=10,
                           nTree=30,
                           na.string="NA",
                           missing=0.05,
                           missingType="MCAR",
                           transformation="None",
                           sampleStart=3){


  outputAccuracyFile <- output

  #### Load complete dataset ##################################################
  data <-read.table(input, header=TRUE,sep="\t")



  #### Parameter for tests ####################################################
  knn_k <-k

  # Number of imputation test
  nTest = nbTest


  #### Matrix of imputation error #############################################
  total <- dim(data)[1]*dim(data)[2]

  # Number of fake missing value to impute in each test
  nMissingPerTest = round(missing*total)

  res <- c()


  #### Fix the random seed for the random number generator ####################
  rng.seed=345

  #### For each test ##########################################################
  for( i in 1:nTest)
  {
    print(paste(" Test run #", toString(i)))

    # Change the random seed
    rng.seed = rng.seed +1

    # Make a copy of the original data matrix, only the intensity column
    testMatrix <- data[, sampleStart:dim(data)[2]]

    for(i in 1:(dim(testMatrix)[2])){
      if((class(testMatrix[,i]) == "factor")|(class(testMatrix[,i]) == "character")){
        warning(paste(paste0("Warning: sampleStart might be wrong: column #", i), " are factors or characters"))
      }
    }
    #### Randomly set n value to missing  ##################################

    if(missingType=="MNAR"){
      if(missing > 0.1){
        warning("number of missing too great, adjusted at 0.05")
        missing <- 0.05
        round(missing*total)
      }
      for(j in 1:nMissingPerTest)
      {
        ok = FALSE
        while(!ok)
        {

          colIdx = sample(1:NCOL(testMatrix), 1, replace = FALSE, prob = NULL)
          rowIdx = sample(1:NROW(testMatrix), 1, replace = FALSE, prob = NULL)

          val = testMatrix[rowIdx,colIdx]
          condition <- (val <= (range(testMatrix[rowIdx,], na.rm=T)[2]-range(testMatrix[rowIdx,], na.rm=T)[1])/2)

           if(!is.na(condition)){
            if((!is.na(val)) & condition)

            {
              ok = TRUE
            }
          }

        }

        testMatrix[rowIdx,colIdx] <- NA

      }




    }else{

      for(j in 1:nMissingPerTest)
      {
        ok = FALSE
        while(!ok)
        {

          colIdx = sample(1:NCOL(testMatrix), 1, replace = FALSE, prob = NULL)
          rowIdx = sample(1:NROW(testMatrix), 1, replace = FALSE, prob = NULL)
          val = testMatrix[rowIdx,colIdx]
          if(!is.na(val))
          {
            ok = TRUE
          }

        }

        testMatrix[rowIdx,colIdx] <- NA

      }
    }



    #### Transformation ####################################################################
    if(transformation=="log"){
      #Add 1 in each cell to overcome issue of 0
      myData1_noZero_temp <- testMatrix + 1
      #Transform intensity data to the log scale
      myDataMatrix <- log(myData1_noZero_temp)
      #Retrieve row log-mean values
      rowMean <- rowMeans(myDataMatrix,na.rm=TRUE)
      #Substract row log-mean from each value (recenter each row to 0)
      for ( i in 1:NROW(myDataMatrix)){
        myDataMatrix[i,] = myDataMatrix[i,] - rowMean[i]
      }
      testMatrix <- myDataMatrix


    }else if(transformation=="scale"){
      testMatrixT <- scale(t(testMatrix))
      testMatrix <- t(testMatrixT)

    }


    testMatrix <- as.matrix(testMatrix)


    #### Do the imputation for this test and retrieve the imputed data
    if (method=="knn"){
      myRes <- impute.knn(testMatrix,k=knn_k,rowmax=0.8)$data
      myDataMatrixImputed <- myRes
      imputedData <- myDataMatrixImputed

    }else if(method=="RF"){
      imputedData <- missForest(testMatrix, decreasing = T, ntree=nTree)$ximp

    }else if(method=="mean"){
      myDataMatrix <- testMatrix
      k <- which(is.na(myDataMatrix), arr.ind=TRUE)
      myDataMatrix[k] <- rowMeans(myDataMatrix, na.rm=TRUE)[k[,1]]
      imputedData <- myDataMatrix

    }else if(method=="median"){
      myDataMatrix <- testMatrix
      k <- which(is.na(myDataMatrix), arr.ind=TRUE)
      myDataMatrix[k] <- rowMedians(myDataMatrix, na.rm=TRUE)[k[,1]]
      imputedData <- myDataMatrix

    }else if(method=="0"){
      myDataMatrix <- testMatrix
      myDataMatrix[is.na(myDataMatrix)] <- 0
      imputedData <- myDataMatrix
    }else if(method=="HM"){
      k <- which(is.na(testMatrix), arr.ind=TRUE)
      testMatrix[k] <- min(testMatrix[k[,1],], na.rm=TRUE)/2
      imputedData <-testMatrix

    }else if(method=="SVD"){

      imputedData<- pca(t(testMatrix), method="svdImpute", nPcs=npcs, center = F)
      imputedData<- t(imputedData@completeObs)

    }else if(method=="QRILC"){
      if(transformation!="log"){
        warning("Log transformation is necessary for this method")
      }
      imputedData<- impute.QRILC(testMatrix, tune.sigma=sigma)[[1]]


    }else{
      stop(paste0("Not a valid method: ", method))
    }


    #### Transform back ################################################################
    if(transformation=="log"){
      for ( i in 1:NROW(imputedData)){
        imputedData[i,] = exp(imputedData[i,] + rowMean[i])

      }

    }else if(transformation=="scale"){

      imputedData <- (imputedData * attr(testMatrixT, 'scaled:scale') +attr(testMatrixT, 'scaled:center'))


    }


    #### compute the NRMSE ########################################################
    dat.impute<- imputedData
    dat.true <- data[,sampleStart:dim(data)[2]]

    NRMSE<- NRMSE(dat.impute, dat.true)
    res <- c(res, NRMSE)


  }

  #### Generates the result matrix #############################################

  resDf <- cbind(Method=method, missing_proportion=missing, transformation=transformation,NRMSE=mean(res))

  if(!is.null(output)){
    #Write the results to file
    write.table(resDf, file =outputAccuracyFile, row.names=FALSE,sep="\t", quote = F)
  }


  return(resDf)

}


#' NRMSE computing
#'
#' Compute NRMSE for accuracy mesure. Compare tow matrices, one being imputed from the original one with holes added.
#'
#' NRMSE is sqrt((sum((dat.impute-dat.true)^2) / totalNbValues ) / variance(dat.true))
#'
#' @param dat.impute imputed dataframe
#' @param dat.true original dataframe
#'
#' @return NRMSE value
#'
#' @examples
#' NRMSE(ImputedDataset, OriginalCompleteDataset)
#'
#' @export
NRMSE <- function(dat.impute, dat.true){
  x <- dat.impute-dat.true
  x2 <- x^2
  somme <- sum(colSums(x2))
  variance <- sd(unlist(dat.true))^2
  total <- dim(dat.true)[1]*dim(dat.true)[2]
  NRMSE <- sqrt((somme/(total))/variance)
  return(NRMSE)
}
