% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sum.stat.raw.R
\docType{data}
\name{sum.stat.raw}
\alias{sum.stat.raw}
\title{GIANT GWAS summary statistics on height around rs11090631 on chromosome 22}
\format{1 row per SNP. The variables are as follows:
\itemize{
   \item SNP, SNP ID
   \item A1, effect allele
   \item A2, reference allele
   \item Freq1, the allele frequency of Allele1
   \item b, estimate of marginal effect in GWAS
   \item se, standard error of the estimates of marginal effects in GWAS
   \item N, sample size
}}
\source{
GWAS Anthropometric 2014 Height by GIANT consortium. 
GIANT homepage: http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
}
\usage{
data(sum.stat.raw)
}
\description{
This data set is the GWAS summary statistics on height reported by GIANT consortium. 
963 SNPs within a 1-Mb centred at rs11090631 are included in the data.
}
\author{
Zheng Ning, 2017-05-18
}
