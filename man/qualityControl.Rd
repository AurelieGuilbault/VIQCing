% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.r
\name{qualityControl}
\alias{qualityControl}
\title{Data Quality Control}
\usage{
qualityControl(data, output = NULL, cohort = "", missing = 0.2,
  compound = NULL, metabolite = NULL, sampleStart = 6)
}
\arguments{
\item{data}{name of the initial metabolomic dataset file.  Must contain:
\itemize{
\item a column "Compound";
\item a column "Metabolite";
\item and all the columns sample from <sampleStart> to the end of the file;
}
the rest doesn't matter and the names "Compound" and "Metabolite" are optional, as long as the column position is specified in params. ;}

\item{output}{default NULL, name of the output cleaned dataset to produce if not NULL;}

\item{cohort}{default "sample", name of the sample's cohort;}

\item{compound}{default NULL, position of the compound column if named otherwise;}

\item{metabolite}{default NULL, position of the metabolite column if named otherwise;}

\item{sampleStart}{default 6, the index of the 1st sample column.}
}
\value{
a list, containing the cleaned dataset and the QC dataframe
\itemize{
\item dataset cleaned dataset
\item QC the QC dataframe containing  the columns compound, metabolite, nbna,	sd,	mu,	CV,	 and remove.

}
}
\description{
Data Quality Control
}
\examples{
for a dataset with the following header ; Compound, m/z, Metabolite, RT, Sample #1, ...
qualityControl("dummySet.txt", sampleStart = 5)

for a dataset with the following header ; compound, m/z, metabolite, RT, Sample #1, ...
qualityControl("dummySet.txt", sampleStart = 5, compound = 1, metabolite = 3)


}
