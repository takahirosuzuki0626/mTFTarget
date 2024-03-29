% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMRPromoter.R
\name{DMRPromoter}
\alias{DMRPromoter}
\title{Extraction of DMR associated promoter}
\usage{
DMRPromoter(sig_targets_gr)
}
\arguments{
\item{sig_targets_gr}{a Granges object}
}
\value{
a GRanges object of DMR horboring genes at their promoter
}
\description{
This function extract a promoter horboring DMR
Gene model used in this function is Gencode_hg19(gene only)
Promoter is defined as 1.5kbp upstream and 500 bp downstream from TSS
}
\keyword{DMR}
\keyword{GRanges}
\keyword{Promoter}
