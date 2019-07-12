options(stringsAsFactors=FALSE)
#' Customization of QC output
#'
#' Enable the user to choose which parameter to output from a dataframe containing the columns compound, metabolite, nbna,	sd,	mu,	CV,	and remove.
#'
#' @param data the input QC dataframe;
#' @param COMPOUND default=TRUE, output the column named "compound";
#' @param NBNA default=TRUE, output the column named "nbna";
#' @param SD default=TRUE, output the column named "sd";
#' @param MU default=TRUE, output the column named "mu";
#' @param CV default=TRUE, output the column named "CV";
#' @param REMOVE default=TRUE, output the column named "remove".
#'
#' @return subdata the desired dataframe,
#'
#'
#' @examples
#' QCcustomization(QCresult$QC, REMOVE=FALSE)
#'
#' @export
QCcustomization <- function(data, COMPOUND=TRUE, COHORTE=TRUE, METABOLITE=TRUE, NBNA=TRUE, SD=TRUE, MU=TRUE, CV=TRUE, REMOVE=TRUE){

  cols <- c(COMPOUND, COHORTE, METABOLITE, NBNA, SD, MU, CV, REMOVE)
  subdata <- cbind(data[,cols])
  return(subdata)

}

#' Data Quality Control
#'
#' Quality control for metabolomic data. Will remove metabolites with more than a determined thershold of missing values for all patient.
#' Produces the following files:
#' \itemize{
#' \item "QC_ <data> .txt": Summary of the QC;
#' \item "REMOVED_QC_<data>.txt": Summary of the removed metabolite;
#' \item output, cleaned dataset file(optional);
#' \item "REMOVED_<output>.txt", the removed set of metabolite (optional)
#'
#' }
#'
#' Can be used on any kind of file containing rows with x (gene, metabolite, etc.) and column with y (Sample, patients, etc.).
#'
#' @param data name of the initial metabolomic dataset file.  Must contain:
#' \itemize{
#' \item a column "Compound";
#' \item a column "Metabolite";
#' \item and all the columns sample from <sampleStart> to the end of the file;
#' }
#' the rest doesn't matter and the names "Compound" and "Metabolite" are optional, as long as the column position is specified in params. ;
#'
#' @param output default NULL, name of the output cleaned dataset to produce if not NULL;
#' @param cohort default "sample", name of the sample's cohort;
#' @param compound default NULL, position of the compound column if named otherwise;
#' @param metabolite default NULL, position of the metabolite column if named otherwise;
#' @param sampleStart default 6, the index of the 1st sample column.
#'
#' @return a list, containing the cleaned dataset and the QC dataframe
#' \itemize{
#' \item dataset cleaned dataset
#' \item QC the QC dataframe containing  the columns compound, metabolite, nbna,	sd,	mu,	CV,	 and remove.
#'
#' }
#'
#' @examples
#' for a dataset with the following header ; Compound, m/z, Metabolite, RT, Sample #1, ...
#' qualityControl("dummySet.txt", sampleStart = 5)
#'
#' for a dataset with the following header ; compound, m/z, metabolite, RT, Sample #1, ...
#' qualityControl("dummySet.txt", sampleStart = 5, compound = 1, metabolite = 3)
#'
#'
#' @export
qualityControl <- function(data, output=NULL, cohort="",
               missing=0.2,
               compound=NULL,
               metabolite=NULL,
               sampleStart=6){

  ###### Read data  #################################################
  DAT <- read.table(data, header=T, sep="\t", na.strings="", fill=TRUE)

  dat <- DAT

  if (!is.null(compound)){
    colnames(dat)[compound] <- "Compound"

  }

  if (!is.null(metabolite)){
    colnames(dat)[metabolite] <- "Metabolite"
  }

  if(is.null(dat$Compound)){
    stop("no Compound column or column number entered")

  }
  if(is.null(dat$Metabolite)){
    stop("no Metabolite column or column number entered")

  }


  dat$Cohort= cohort
  GENMOD <- dat$Cohort

  # les 5 premieres colonnes donnent les infos sur le metabolite
  IDS <- data.frame(Compound=dat$Compound, Metabolite=dat$Metabolite)

  METABO <- dat[,((sampleStart+1):dim(dat)[2]-1)]
  METABO <- suppressWarnings(data.matrix(METABO))

  ######### QC data  et sommaire ##############

  # remplacer les 0 par NA
  METABO[METABO == 0] <- NA

  # compte des NA par groupe
  nbna <- apply(is.na(METABO),1, sum)

  # standard deviation par groupe
  sd <- apply(METABO,1, sd, na.rm=T)

  # moyenne par groupe

  mu <- apply(METABO,1 , mean, na.rm=T)


  # coefficient de variation
  CV <- sd/mu

  # on identifie les metabolites avec plus de 20% de valeurs manquantes par groupe
  rm  <- (nbna > (missing * dim(METABO)[2]))

  # Build the name of the output QC_ file
  name<- paste0("QC_", strsplit(data, ".", fixed=TRUE)[[1]][1])

  QC <- data.frame(compound =IDS$Compound , cohort=GENMOD, metabolite=IDS$Metabolite, nbna, sd, mu, CV, remove=rm)
  print("Saving QC data")

  write.table(QCcustomization(QC[QC$remove == FALSE,], REMOVE=FALSE), paste0(name,".txt"), quote=F, row.names=F,sep="\t")
  write.table(QCcustomization(QC[QC$remove==TRUE,], REMOVE=FALSE), paste0("REMOVED_",paste0(name,".txt")), quote=F, row.names=F,sep="\t")

  # Retirer les métabolites concernés
  clean_dataset  <- as.data.frame(METABO)

  clean_dataset$remove=QC$remove

  remove_set <- clean_dataset[clean_dataset$remove==TRUE,]
  clean_dataset <- clean_dataset[clean_dataset$remove==FALSE,]

  # Ré-insérer les IDS dans le fichier et générer le nouveau dataset
  IDS$remove=QC$remove


  clean_dataset  <- cbind(IDS[IDS$remove==FALSE, c(1,2)], clean_dataset[,1:(dim(clean_dataset)[2]-1)])
  removed_set <- cbind(IDS[IDS$remove==TRUE, c(1,2)], remove_set[,1:(dim(remove_set)[2]-1)])
  # Output le cleaned_dataset si demandé
  if(!(is.null(output))){
    print("Saving cleaned data")
    write.table(clean_dataset, output, quote=F, row.names=F,sep="\t")
    write.table(removed_set, paste0("REMOVED_",output), quote=F, row.names=F,sep="\t")
  }



  label<- paste(QC[QC$remove==FALSE,"compound"],
       QC[QC$remove==FALSE,"metabolite"], sep = "_")
  dup_label<- duplicated(label)

  i <- 0
  for (dup in dup_label){
    i <- i+1

    if(dup){
      warning(paste(paste(paste("row", i), "is a duplicated compound_metabolite:"), label[i]))
    }

  }


  j <- 0
  for(i in unlist(QC[QC$remove==FALSE,]$sd)){
    j<- j+1
    if(is.na(i)){
      warning(paste0("Sd == NA for row: ", label[j]))
    }else if(i==0){
      warning(paste0("Sd == 0 for row: ", label[j]))

    }

  }

  return(list("dataset"= clean_dataset, "QC"=QC))
}
