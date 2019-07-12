#### Library imports #########################################################
library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggpubr)
library(scales)
###############################################################################


#' Simple Violin Plot Constructor
#'
#' Construct Violin plots for a metabolomic dataset cleaned by QC.r (not required).
#' Produces a pdf File containing the plots.
#'
#' Can be used on any kind of file containing rows with x (gene, metabolite, etc.) and column with y (Sample, patients, etc.).
#'
#' @param input name of the cleaned datafile;
#'     Must contain:
#'      -a column named "Compound";
#'      -a column named "Metabolite";
#'      -and all the columns sample from <sampleStart> to the end of the file;
#'       the rest doesn't matter and the names are optional, as long as the column position is entered.
#' @param output default <input>.pdf, name of the final pdf file;
#' @param na default FALSE, FALSE removes the untargeted meatabolite from the plots;
#' @param numberByGraph default 1, number of metabolite violin plots by graph window;
#' @param colNumberByPage default 5, number of column in the output pdf file;
#' @param rowNumberByPage default 2, number of row in the output pdf file;
#' @param landscape default FALSE, orientation of the pdf file;
#' @param compound default NULL, position of the compound column if named otherwise;
#' @param metabolite default NULL, position of the metabolite column if named otherwise;
#' @param sampleStart default 3, 1st column of the compounds intensity.
#'
#' @return None
#'
#' @examples
#' for a dataset with the following header ; Compound, m/z, Metabolite, RT, Sample #1, ...
#' violinPlotQC<- function("dummySet.txt", sampleStart = 5)
#'
#' for a dataset with the following header ; compound, m/z, metabolite, RT, Sample #1, ...
#' violinPlotQC<- function("dummySet.txt", compound = 1, metabolite = 3, sampleStart = 5)
#'
#' @import ggplot2
#' @import tidyr
#' @import gridExtra
#' @import ggpubr
#' @import scales
#'
#'@seealso \itemize{
#' \item \code{ggplot2} package\url{https://www.rdocumentation.org/packages/ggplot2}
#' \item \code{tidyr} package \url{https://www.rdocumentation.org/packages/tidyr}
#' \item \code{gridExtra} package \url{https://www.rdocumentation.org/packages/gridExtra}
#' \item \code{ggpubr} package \url{https://www.rdocumentation.org/packages/ggpubr}
#' \item \code{scales} package \url{https://www.rdocumentation.org/packages/scales}
#' }
#' @export
violinPlotQC<- function(input,
                           output=paste0(strsplit(input, ".",
                                                  fixed=TRUE)[[1]][1], ".pdf"),
                           na=F,
                           numberByGraph=1,
                           colNumberByPage=5,
                           rowNumberByPage=2,
                           landscape=F,
                           compound=NULL,
                           metabolite=NULL,
                           sampleStart=3){
  #### Data Loading ############################################################
  df<- read.table(input, header=T, sep="\t")

  if (!is.null(compound)){
    colnames(df)[compound] <- "Compound"

  }

  if (!is.null(metabolite)){
    colnames(df)[metabolite] <- "Metabolite"
  }

  if(is.null(df$Compound)){
    stop("no Compound column or column number entered")

  }
  if(is.null(df$Metabolite)){
    stop("no Metabolite column or column number entered")

  }
  for(i in sampleStart:(dim(df)[2])){
    if((class(df[,i]) == "factor")|(class(df[,i]) == "character")){
      warning(paste(paste0("Warning: sampleStart might be wrong: column #", i), " are factors or characters"))
    }
  }
  #### Remove untargeted metabolite if asked ###################################
  if(!na){
    df <- df[!is.na(df$Metabolite),]
  }

  #### Create the ID  ##########################################################
  df<- unite(df, "ID", c("Compound", "Metabolite"), sep="\n")
  rownames(df) <- df$ID
  df$ID <- (as.factor(df$ID)) # For the Violin plots


  print("Computing stats")
  # standard deviation
  sd <- (apply(df[,sampleStart:dim(df)[2]],1, sd, na.rm=T))


  # mean
  mu <- (apply(df[,sampleStart:dim(df)[2]],1 , mean, na.rm=T))
  df$mu <-mu



  # variation coefficient
  CV <- sd/mu
  df$CV <- CV
  NAwarning <- FALSE
  j <- 0
  for(i in CV){
    j<- j+1
    if(is.na(i)){
      NAwarning <- TRUE
      warning(paste0("CV == NA for row: ", df[j, "ID"]))
    }else if(i==0){
      warning(paste0("CV == 0 for row: ", df[j, "ID"]))

    }
  }

  #### Option landscape or portrait ############################################
  if(landscape){
    pdf(file=paste(getwd(), output, sep = "/"), 20, 20,onefile = T,
        paper = "a4r")
  }else{
    pdf(file=paste(getwd(), output,sep = "/"), 20, 20,onefile = T)
  }


  print("Plotting")
  #### Option for all metabolite on a single plot ###############################
  if (class(numberByGraph) != class(0)){

    #Format the data for the plot

    METABO <- gather(df, Sample, Value, sampleStart:(dim(df)[2]-2))
    METABO$Value <- sapply(METABO$Value,as.numeric)
    METABO$CV <- round(METABO$CV, 4)

    # Plotting
    p <- ggplot(METABO, aes(x=ID, y=Value))

    # Get the y axis minimum for the format of the CV
    gb = ggplot_build(p)
    ymin = gb$layout$panel_params[[1]]$y.range[1]

    # Add the Violins
    p <- (p + geom_violin(trim=T) +
      # Boxplot
      geom_boxplot(width=0.05) +
      # Stat summary
      stat_summary(fun.y=mean, geom="point", color="red") +
      # x and y axis setting
        # x and y axis settings
        theme(axis.text.x = element_text(angle = -90,
                                         hjust=0, vjust = 0.8, size = 6)) +
        scale_y_continuous(labels = comma) +
      # CV text settings
      annotate("text", x=METABO$ID, y=ymin, label=paste("CV:" , METABO$CV), size=1.5, angle = -90) +

      # Main title
      ggtitle(paste("Distribution of metabolites for ",
                    strsplit(input, ".", fixed=TRUE)[[1]][1], sep="\n")) +
      theme(plot.title = element_text(face="bold", hjust = 0.5, size=10))+
      # General theme
      theme(panel.border =
              element_blank(),
            panel.grid.major =
              element_line(size = 0.5, linetype = 'solid', colour = "grey"),
            panel.grid.minor =
              element_line(size = 0.5, linetype = 'solid', colour = "grey"),
            panel.background = element_blank(),
            axis.line = element_line(colour = "grey")) +
      # Remove the x and y axis titles
      labs(x= "", y="")
    )
    print(p)

  #### Violin plot will be divided by <numberByGraph> and #######################
  ####<colNumberByPage> x <rowNumberByPage>/page ################################
  }else{

    # Verify for the numberByGraph to adjust the Titles and make sure there's no 0
    if(numberByGraph<=0){
      print("enter a valid numberByGraph (must be > 0)")

    }else{
      # Counter for the list index
      k <-0
      liste= list()

      # Create graph with correct number of plot by graph
      for(i in seq(from=1, to=dim(df)[1], by= numberByGraph)){
        k <- (k+1)

        if((i+numberByGraph-1) > (dim(df)[1])){
          j <- (dim(df)[1])
        }else{
          j <- (i+numberByGraph-1)
        }

        # Select the data
        dat <- df[i:j,]

        #Format the data for the plot
        METABO <- gather(dat, Sample, Value, sampleStart:(dim(df)[2]-2))
        METABO$Value <- sapply(METABO$Value,as.numeric)
        METABO$CV <- round(METABO$CV, 4)

        #plotting
        p <- ggplot(METABO, aes(x=ID, y=Value))

        # Get the y axis minimum for the format of the CV
        gb = ggplot_build(p)
        ymin = gb$layout$panel_params[[1]]$y.range[1]

        # Add the Violins
        p <- p + geom_violin(trim=T) +
          # Boxplot
          geom_boxplot(width=0.05) +
          # Stats summary
          stat_summary(fun.y=mean, geom="point", color="red") +
          # x and y axis settings
          theme(axis.text.x = element_text(angle = -90,
                                           hjust=0, vjust = 0.8, size = 6)) +
          scale_y_continuous(labels = comma) +
          # General theme
          theme(panel.border = element_blank(),
                panel.grid.major =
                  element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                panel.grid.minor =
                  element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                panel.background = element_blank(),
                axis.line = element_line(colour = "grey")) +

          # Remove the x and y labels
          labs(x= "", y="")

        # Adjust the main title
        if(numberByGraph==1){
          p <- p + ggtitle(METABO$ID)+
                theme(plot.title = element_text(face="bold", hjust = 0.5, size=20)) +
                # CV text settings
                annotate("text", x=METABO$ID, y=ymin,
                         label=paste("CV:", METABO$CV), size=4)


        }else{
          p <- p + ggtitle(paste("Distribution of metabolites for ",
                        strsplit(input, ".", fixed=TRUE)[[1]][1], sep="\n")) +
                theme(plot.title = element_text(face="bold", hjust = 0.5, size=20)) +
                # CV text settings
                annotate("text", x=METABO$ID, y=ymin,
                         label=paste("CV:", METABO$CV), size=3, angle= -90)

              }


        # Add to plot list
        liste[[k]] <- p

      }

      # Print in PDF with right number of plots by page
      print("Saving PDF file")
      print(ggarrange(plotlist = liste, ncol = colNumberByPage,
                nrow= rowNumberByPage, align = "hv"))
    }
  }
  dev.off()
  if (NAwarning){
    print("It's possible the sampleStart causes CV == NA for all metabolite, or there's a line with only NAs")
  }
}


#' Comparative Violin Plot Constructor
#'
#' Construct Violin plots from imputed dataset and cleaned by QC.r.
#' Produces a pdf File containing the plots.
#'
#' it is possible to put any 2 corresponding files, as long as they have the same row and column names. (The imputation part is not necessary, but it is the principal input type of this function. The legend and labels does contain the term "Imputed" by default.)
#'
#' @param input name of the cleaned datafile;
#'     Must contain:
#'      -a column named "Compound";
#'      -a column named "Metabolite";
#'      -and all the columns sample from <sampleStart> to the end of the file;
#'       the rest doesn't matter and the names are optional, as long as the column position is entered..
#' @param inputImputed name of the imputed datafile;
#'     Must contain:
#'      -a column named "Compound";
#'      -a column named "Metabolite";
#'      -and all the columns sample from <sampleStart> to the end of the file;
#'       the rest doesn't matter and the names are optional, as long as the column position is entered.
#' @param output default <input>.pdf, name of the final pdf file;
#' @param na default FALSE, FALSE removes the untargeted meatabolite from the plots;
#' @param numberByGraph default 1, number of metabolite violin plots by graph window;
#' @param colNumberByPage default 5, number of column in the output pdf file;
#' @param rowNumberByPage default 2, number of row in the output pdf file;
#' @param landscape default FALSE, orientation of the pdf file;
#' @param compound default NULL, position of the compound column if named otherwise;
#' @param metabolite default NULL, position of the metabolite column if named otherwise;
#' @param sampleStart default 3, 1st column of the compounds intensity.
#' @param compoundImp default NULL, position of the compound column for Imputed dataset if named otherwise;
#' @param metaboliteImp default NULL, position of the metabolite column for Imputed dataset if named otherwise;
#' @param sampleStartImp default 3, 1st column of the compounds intensity for Imputed dataset.
#' @param onlyImputed default TRUE, TRUE will display only plots for metabolites with imputed values. (nb_imputed > 0)
#'
#' @return None
#'
#' @import ggplot2
#' @import tidyr
#' @import gridExtra
#' @import ggpubr
#' @import scales
#'
#' @examples
#' for 2 dataset with the following headers ; Compound, Metabolite, Sample #1, ...
#' violinPlotQC<- function("dummySet.txt")
#'
#' for an imputed dataset with the following header ; compound, metabolite,  Sample #1, ...
#' violinPlotQC<- function("dummySet.txt", compoundImp = 1, metaboliteImp = 2)
#'
#'@seealso \itemize{
#' \item \code{ggplot2} package\url{https://www.rdocumentation.org/packages/ggplot2}
#' \item \code{tidyr} package \url{https://www.rdocumentation.org/packages/tidyr}
#' \item \code{gridExtra} package \url{https://www.rdocumentation.org/packages/gridExtra}
#' \item \code{ggpubr} package \url{https://www.rdocumentation.org/packages/ggpubr}
#' \item \code{scales} package \url{https://www.rdocumentation.org/packages/scales}
#' }
#' @export
violinPlotImp<- function(input,
                    inputImputed=paste0(strsplit(input, ".", fixed=TRUE)[[1]][1],"_imputed.txt"),
                    output=paste0(strsplit(input, ".", fixed=TRUE)[[1]][1], "_Imputed.pdf"),
                    na=F,
                    numberByGraph=1,
                    colNumberByPage=5,
                    rowNumberByPage=2,
                    landscape=F,
                    compound=NULL,
                    metabolite=NULL,
                    sampleStart=3,
                    compoundImp=NULL,
                    metaboliteImp=NULL,
                    sampleStartImp=3,
                    onlyImputed=TRUE){

  #### Data Loading and formating ############################################################
  df_avant <- read.table(input, header=T, sep="\t")

  if (!is.null(compound)){
    colnames(df_avant)[compound] <- "Compound"

  }

  if (!is.null(metabolite)){
    colnames(df_avant)[metabolite] <- "Metabolite"
  }

  if(is.null(df_avant$Compound)){
    stop("no Compound column or column number entered for input set")

  }
  if(is.null(df_avant$Metabolite)){
    stop("no Metabolite column or column number entered for input set")

  }

  for(i in sampleStart:(dim(df_avant)[2])){
    if((class(df_avant[,i]) == "factor")|(class(df_avant[,i]) == "character")){
      warning(paste(paste0("Warning: sampleStart might be wrong: column #", i), " are factors or characters"))
    }
  }
  df_avant$Imp <- "Not imputed"
  df_avant$nb_na <- apply(df_avant[ ,(sampleStart:dim(df_avant)[2])], 1, function(x) sum(is.na(x)))


  df_imputed <- read.table(inputImputed, header=T, sep="\t")

  if (!is.null(compoundImp)){
    colnames(df_imputed)[compoundImp] <- "Compound"
  }

  if (!is.null(metaboliteImp)){
    colnames(df_imputed)[metaboliteImp] <- "Metabolite"
  }

  if(is.null(df_imputed$Compound)){
    stop("no Compound column or column number entered for imputed set")

  }
  if(is.null(df_imputed$Metabolite)){
    stop("no Metabolite column or column number entered for imputed set")

  }

  for(i in sampleStartImp:(dim(df_imputed)[2])){
    if((class(df_imputed[,i]) == "factor")|(class(df_imputed[,i]) == "character")){
      warning(paste(paste0("Warning: sampleStartImp might be wrong: column #", i), " are factors or characters"))
    }
  }

  df_imputed$nb_na <- apply(df_avant[ ,(sampleStartImp:dim(df_avant)[2])], 1, function(x) sum(is.na(x)))



  df_imputed$Compound <- paste(df_imputed$Compound, "imputed")
  df_imputed$Imp <- "Imputed"

  if((dim(df_imputed)[2] - sampleStartImp +1)!=((dim(df_avant)[2] - sampleStart +1))){

    stop("Number of sample columns is unequal for both dataset")
  }
  #### Create the combined df #############################################################
  df <- rbind(df_avant[, c(which(colnames(df_avant)=="Compound"), which(colnames(df_avant)=="Metabolite"), sampleStart:dim(df_avant)[2])],
                df_imputed[, c(which(colnames(df_imputed)=="Compound"), which(colnames(df_imputed)=="Metabolite"),
                               sampleStartImp:dim(df_imputed)[2])])


  df$nb_na <- paste("nb_imp:", df$nb_na)
  df$nb_na[df$nb_na == "nb_imp: 0"] <- ""

  if (onlyImputed){
    df <- df[!(df$nb_na == ""),]
  }
  df$nb_na[df$Imp == "Not imputed"] <- ""

  sampleStart <- 3


  #### Remove untargeted metabolite if asked ###################################
  if(!na){
    df <- df[!is.na(df$Metabolite),]
  }

  #### Create the ID  ##########################################################
  df<- unite(df, "ID", c("Compound", "Metabolite"), sep="\n")
  df$ID <- as.factor(df$ID) # For the Violin plots


  print("Computing stats")
  # standard deviation
  sd <- apply(df[,sampleStart:(dim(df)[2]-2)],1, sd, na.rm=T)

  # mean
  mu <- apply(df[,sampleStart:(dim(df)[2]-2)],1 , mean, na.rm=T)
  df$mu <-mu

  # variation coefficient
  CV <- sd/mu
  df$CV <- CV

  NAwarning <- FALSE
  j <- 0
  for(i in CV){
    j<- j+1
    if(is.na(i)){
      NAwarning <- T
      warning(paste0("CV == NA for row: ", df[j, "ID"]))
    }else if(i==0){
      warning(paste0("CV == 0 for row: ",df[j, "ID"]))

    }
  }
  #### Nombre de metabolite ####################################################
  rowLength <- (dim(df)[1])/2

  #### Option landscape or portrait ############################################
  if(landscape){
    pdf(file=paste(getwd(), output, sep = "/"), 20, 20,onefile = T,
        paper = "a4r")
  }else{
    pdf(file=paste(getwd(), output, sep = "/"), 20, 20,onefile = T)
  }


  print("Plotting")
  #### Option for all metabolite on a single plot ###############################
  if (class(numberByGraph) != class(0)){

    #Format the data for the plot

    METABO <- gather(df, Sample, Value, sampleStart:(dim(df)[2]-4))
    METABO$Value <- sapply(METABO$Value,as.numeric)
    METABO$CV <- round(METABO$CV, 4)

    # Plotting
    p <- ggplot(METABO, aes(x=ID, y=Value, fill=Imp))

    # Get the y axis minimum for the format of the CV
    gb = ggplot_build(p)
    ymin = gb$layout$panel_params[[1]]$y.range[1]

    # Add the Violins
    p <- (p + geom_violin(trim=T) +
            # Boxplot
            geom_boxplot(width=0.05) +
            # Stat summary
            stat_summary(fun.y=mean, geom="point", color="red") +
            # x and y axis setting
            # x and y axis settings
            theme(axis.text.x = element_text(angle = -90,
                                             hjust=0, vjust = 0.8, size = 6)) +
            scale_y_continuous(labels = comma) +
            # CV text settings
            annotate("text", x=METABO$ID, y=ymin, label=paste(METABO$CV, METABO$nb_na, sep="\n"), size=1.5, angle = -90) +
            # Main title
            ggtitle(paste("Distribution of metabolites for ",
                          strsplit(input, ".", fixed=TRUE)[[1]][1], sep="\n")) +
            theme(plot.title = element_text(face="bold", hjust = 0.5, size=10))+
            # General theme
            theme(panel.border =
                    element_blank(),
                  panel.grid.major =
                    element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                  panel.grid.minor =
                    element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "grey")) +
            # Remove the x and y axis titles
            labs(x= "", y="")
    )
    print("Saving PDF file")
    print(p)

    #### Violin plot will be divided by <numberByGraph> and #######################
    ####<colNumberByPage> x <rowNumberByPage>/page ################################
  }else{

    # Verify for the numberByGraph to adjust the Titles and make sure there's no 0
    if(numberByGraph<=0){
      print("enter a valid numberByGraph (must be > 0)")

    }else{
      # Counter for the list index
      k <-0
      liste= list()


      # Create graph with correct number of plot by graph
      for(i in seq(from=1, to=rowLength, by= numberByGraph)){
        k <- (k+1)

        if((i+numberByGraph-1) > (dim(df)[1])){
          j <- rowLength
        }else{
          j <- (i+numberByGraph-1)
        }

        # Select the data
        dat <- df[c((i:j), ((i+rowLength):(j+rowLength))),]

        #Format the data for the plot
        METABO <- gather(dat, Sample, Value, sampleStart:(dim(df)[2]-4))
        METABO$Value <- sapply(METABO$Value,as.numeric)
        METABO$CV <- round(METABO$CV, 4)


        #plotting
        p <- ggplot(METABO, aes(x=ID, y=Value, fill=Imp))

        # Get the y axis minimum for the format of the CV
        gb = ggplot_build(p)
        ymin = gb$layout$panel_params[[1]]$y.range[1]

        # Add the Violins
        p <- p + geom_violin(trim=T) +
          # Boxplot
          geom_boxplot(width=0.05) +
          # Stats summary
          stat_summary(fun.y=mean, geom="point", color="red") +
          # x and y axis settings
          theme(axis.text.x = element_text(angle = -90,
                                           hjust=0, vjust = 0.8, size = 6)) +
          scale_y_continuous(labels = comma) +
          # General theme
          theme(panel.border = element_blank(),
                panel.grid.major =
                  element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                panel.grid.minor =
                  element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                panel.background = element_blank(),
                axis.line = element_line(colour = "grey")) +

          # Remove the x and y labels
          labs(x= "", y="")

        # Adjust the main title
        if(numberByGraph==1){
          p <- p + ggtitle(METABO$ID)+
            theme(plot.title = element_text(face="bold", hjust = 0.5, size=20)) +
            # CV text settings
            annotate("text", x=METABO$ID, y=ymin,
                     label=paste(paste("CV", METABO$CV), METABO$nb_na, sep="\n"),
                     size=3, angle=-90)


        }else{
          p <- p + ggtitle(paste("Distribution of metabolites for ",
                                 strsplit(input, ".", fixed=TRUE)[[1]][1], sep="\n")) +
            theme(plot.title = element_text(face="bold", hjust = 0.5, size=20)) +
            # CV text settings
            annotate("text", x=METABO$ID, y=ymin,
                     label=paste(paste("CV:", METABO$CV), METABO$nb_na, sep="\n"),
                     size=2, angle= -90)

        }


        # Add to plot list
        liste[[k]] <- p

      }

      # Print in PDF with right number of plots by page
      print("Saving PDF file")
      print(ggarrange(plotlist = liste, ncol = colNumberByPage,
                      nrow= rowNumberByPage, align = "hv"))
    }
  }
  dev.off()
  if (NAwarning){
    warning("It's possible the sampleStart causes CV == NA for all metabolite, or there's a line with only NAs")
  }
}


