#### Imports ####
library(reshape2)
library(Hmisc)
library(gplots)
#################

#' Correlation Matrix and tree
#'
#' Creates the correlation matrix and tree for a metabolomic dataset.
#' Produces:
#' - <filename>.pdf, containing the correlation matrix
#' - <filename>_pairs.txt, containing the correlation pairs above <upperLimit> or less than <lowerLimit>
#'
#' The significance of the correlation is noted with a "*" in the tile.
#' Conditions for the Pearson's Correlation Test:
#'   - Independent samples;
#'   - normal distribution of the data
#'       *Since this condition is not true for most metabolite,
#'        the test is set as a Spearman's rank Correlation test by default.*
#'
#' @param filename name of the datafile
#'     Must contain:
#'       -a column named "Compound";
#'       -a column named "Metabolite";
#'       -and all the columns sample from <sampleStart> to the end of the file;
#'       the rest doesn't matter and the names are optional, as long as the column position is
#'       entered.;
#' @param output default <dataset filename>.pdf, name of the pdf file;
#' @param na default FALSE, FALSE removes the untargeted meatabolite;
#' @param landscape default FALSE, orientation of the pdf file;
#' @param compound default NULL, position of the compound column if named otherwise;
#' @param metabolite default NULL, position of the metabolite column if named otherwise;
#' @param sampleStart default 3, 1st column of the actual data.
#' @param corDigit default 3,number of digit for the correlation coefficient;
#' @param testType default spearman, type of correlation test. Can be "spearman"or "pearson";
#' @param plotWidth default 0.5, width of the top space above the heatmap;
#' @param plotHeigth default 0.2, width of the side space beside the heatmap;
#' @param textSize default 2, size of the text in the heatmap tiles and the dendrogram;
#' @param upperLimit default 0.6, positive correlation limit to output the metabolite pair;
#' @param lowerLimit default -0.6, negative correlation limit to output the metabolite pair.
#'
#' @return the rcorr objet, with:
#' \itemize{
#' \item \code{r} The correlation matrix
#' \item \code{P} The P-value matrix
#' }
#'
#' @examples
#' for a dataset with the following header ; Compound, m/z, Metabolite, RT, Sample #1, ...
#' corMatrix("dummySet.tsv", sampleStart=5, testType="spearman")
#'
#' for a dataset with the following header ; compound, m/z, metabolite, RT, Sample #1, ...
#' corMatrix("dummySet.tsv", compound=1, metabolite=3, sampleStart=5, testType="spearman")
#'
#' @import reshape2
#' @import Hmisc
#' @import gplots
#'
#' @seealso \itemize{
#' \item \code{reshape2} package\url{https://www.rdocumentation.org/packages/reshape2}
#' \item \code{Hmisc} package \url{https://www.rdocumentation.org/packages/Hmisc}
#' \item \code{gplot} package \url{https://www.rdocumentation.org/packages/gplot}
#' }
#' @export
corMatrix <- function(filename, output=paste0(strsplit(filename, ".", fixed=TRUE)[[1]][1], "_Correlation.pdf"),
                      na=F,
                      landscape=F,
                      dendrogram=T,
                      heatmap=T,
                      compound=NULL,
                      metabolite=NULL,
                      sampleStart=3,
                      corDigit=3,
                      testType="spearman",
                      plotWidth=0.5,
                      plotHeigth=0.2,
                      textSize=0.3,
                      upperLimit=0.6,
                      lowerLimit=-0.6){
  #### Load Data #####################################################
  DATA <- read.table(filename, header=T, sep="\t", na.strings="NA", fill=TRUE)
  if (!is.null(compound)){
    colnames(DATA)[compound] <- "Compound"

  }

  if (!is.null(metabolite)){
    colnames(DATA)[metabolite] <- "Metabolite"
  }

  if(is.null(DATA$Compound)){
    stop("no Compound column or column number entered")

  }
  if(is.null(DATA$Metabolite)){
    stop("no Metabolite column or column number entered")

  }
  #### Remove untargeted metabolite if asked ###################################
  if(!na){
    DATA <- DATA[!is.na(DATA$Metabolite),]
  }

  #hypothesis
  # H0: the metabolites are not linearly correlated rho = 0
  # H1: the metabolites are linearly correlated rho != 0


  # Application conditions
  # - Independent samples -> assumed
  # - normal Distribution -> assumed that they are not, so we use the spearman test by default



  data <- DATA[, sampleStart:dim(DATA)[2]]

  for(i in 1:(dim(data)[2])){

    if((class(data[,i]) == "factor")|(class(data[,i]) == "character")){
      warning(paste(paste0("Warning: sampleStart might be wrong: column #", i), " are factors or characters"))
    }
  }
  names <- paste(DATA$Compound, DATA$Metabolite, sep="_")

  #### Correlation matrix ################################################
  rownames(data) <- names


  data <- as.matrix(t((data)))

  if(testType =="pearson"){
    warning("Application conditions for Pearson's Correlation Test
   - Independent samples -> assumed
   - normal Distribution -> assumed ")
  }else{
    warning("Application conditions for Spearman's Correlation test
       - Independent samples -> assumed")

  }
  data.rcorr <- rcorr(data, type=testType)
  data.coeff <- data.rcorr$r
  data.pvalue <- data.rcorr$P

  #### Heatmap building ########################################################
  coeff <- data.coeff
  pvalue<- data.pvalue

  data.coeff[is.na(data.coeff)]<-0
  # Correlation Clustering
  v <- hclust(dist(data.coeff), method = "complete")

  # Formatiing cell label
  coeff <- data.coeff
  pvalue<- data.pvalue

  coeff <- round(coeff, digit=corDigit)
  coeff[coeff == 1] <- NA
  coeff[(coeff < upperLimit)&(coeff > lowerLimit)] <- ""
  coeff[is.na(coeff)] <- ""

  pvalue[pvalue > 0.05] <- NA
  pvalue[is.nan(pvalue)] <- NA
  pvalue[pvalue <= 0.05] <- "*"
  pvalue[is.na(pvalue)] <- ""


  label= matrix(paste( pvalue,coeff, sep="\n"),dim(data.coeff)[1], dim(data.coeff)[1])



  if(dendrogram|heatmap){
    #### Option landscape or portrait ############################################
    if(landscape){
      pdf(file=paste(getwd(), output, sep = "/"), 20, 20,onefile = T,
          paper = "a4r")
    }else{
      pdf(file=paste(getwd(), output, sep = "/"), 20, 20,onefile = T)
    }

    if(heatmap){
      heatmap.2(data.coeff,
                main = paste("Correlation matrix for",
                             strsplit(filename, ".", fixed=TRUE)[[1]][1]), # heat map title
                col= bluered,
                cellnote = label,
                notecol="black",      # change font color of cell labels to black
                notecex = textSize,
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(12,9),     # widens margins around plot
                dendrogram="row",     # only draw a row dendrogram
                as.dendrogram(v),
                key.title = NA,
                key.xlab = "Correlation Coefficient",
                srtCol=315, adjCol = c(0,1),
                srtRow=315, adjRow = c(0,1),
                lhei=c(plotHeigth,4),
                lwid=c(plotWidth,4),
                key.par = list(cex=0.5),
                Colv="Rowv")
    }


    if(dendrogram){
      par(cex=textSize, cex.main= 1/textSize)
      plot(v,
           xlab="",
           ylab="",
           main=paste("Clustering for",
                      strsplit(filename, ".", fixed=TRUE)[[1]][1]),
           sub="", axes=F)

    }
    dev.off()

  }


  melted_cor<- melt(data.coeff)
  melted_p_value <- melt(data.pvalue)


  dat <- data.frame(metabolite1=melted_cor$Var1, metabolite2=melted_cor$Var2,
                    cor=melted_cor$value, pvalue=melted_p_value$value)

  dat <- dat[dat$pvalue <= 0.05,]
  dat[dat==0] <- "<2.22507e-308"
  dat <- dat[((dat$cor >= upperLimit)|(dat$cor <= lowerLimit)),]

  write.table(((dat[order(dat$cor, na.last=NA), ])[c(TRUE, FALSE),]), file=paste0(strsplit(filename, ".", fixed=TRUE)[[1]][1], "_pairs.txt"),
                quote=F, row.names=F,sep="\t")
  return(data.rcorr)
}


