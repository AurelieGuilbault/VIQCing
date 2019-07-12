Aurelie Guilbault
12 July 2019
Last edited: 12 July 2019
VIQCing package
Programming Language: R

The VIQCing (Visualization, Imputation, Quality Control) package has been made to help processing Metabolomic Data. It contains a robust pipeline that can be adapted depending on the input type and it needs. 

Tools: 
- QC(data, output=NULL, cohort="", 
               missing=0.2,
               compound=NULL,
               metabolite=NULL,
               sampleStart=6)
 
- violinPlotQC(input, 
                           output=paste0(strsplit(input, ".", 
                                                  fixed=TRUE)[[1]][1], ".pdf"), 
                           na=F,
                           numberByGraph=1, 
                           colNumberByPage=5, 
                           rowNumberByPage=2, 
                           landscape=F,
                           compound=NULL,
                           metabolite=NULL,
                           sampleStart=3)

- violinPlotImp(input, 
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
                    onlyImputed=TRUE)

- imputation(file, 
                       k=2, method="knn",
                       npcs=3,
                       sigma=0.1,
                       nTree=30,
                       na.string="NA", 
                       transformation="None",
                       compound=NULL,
                       metabolite=NULL,
                       sampleStart=3)

 - imputationTest(output=NULL, input,
                           k=2, method="knn",
                           npcs=3,
                           sigma=0.1,
                           nbTest=10, 
                           nTree=30,
                           na.string="NA", 
                           missing=0.05, 
                           missingType="MCAR",
                           transformation="None",
                           sampleStart=3)

- corMatrix (filename, output=paste0(strsplit(filename, ".", fixed=TRUE)[[1]][1], "_Correlation.pdf"),
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
                      lowerLimit=-0.6)





