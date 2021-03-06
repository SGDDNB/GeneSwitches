% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_logistic_regression.R
\name{downsample_zeros}
\alias{downsample_zeros}
\title{Random downsampling of zero expression}
\usage{
downsample_zeros(glmdata, ratio_ds = 0.7)
}
\arguments{
\item{glmdata}{binary data of one gene}

\item{ratio_ds}{downsampling zeros to this proportion}
}
\value{

}
\description{
This function does randomly downsampling of cells with zero expression for one gene
}
