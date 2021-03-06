% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detTSR.R
\name{detTSR}
\alias{detTSR}
\alias{detTSR,tssObject,character,numeric,numeric,numeric-method}
\title{detTSR}
\usage{
detTSR(experimentName, tssSetType, tssSet = 1, tagCountThreshold, clustDist)

\S4method{detTSR}{tssObject,character,numeric,numeric,numeric}(
  experimentName,
  tssSetType,
  tssSet = 1,
  tagCountThreshold = 1,
  clustDist
)
}
\arguments{
\item{experimentName}{- a S4 object of class tssObject containing
information in slot tssTagData}

\item{tssSetType}{- specifies the set to be clustered. Options are
"replicates" or "merged"}

\item{tssSet}{- number of the dataset to be analyzed}

\item{tagCountThreshold}{- number of TSSs required at a given position}

\item{clustDist}{- maximum distance of TSSs between two TSRs (in base pairs)}
}
\value{
via the user-level function determineTSR, creates a list of
GenomicRanges objects containing TSR positions in slot 'tsrData' on
the tssObject object
}
\description{
An internal function, which is invoked using the user-level
function determineTSR that identifies TSRs from the selected tssSet 
(Internal function)
}
