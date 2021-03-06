% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_switch_pathway.R
\name{merge_pathways}
\alias{merge_pathways}
\title{Merge redundant pathways}
\usage{
merge_pathways(switch_pw, pathways, ratio)
}
\arguments{
\item{switch_pw}{a data frame with significantly changed pathways}

\item{pathways}{a list of pathways with genes}

\item{ratio}{cutoff ratio for merging redundant pathways}
}
\value{

}
\description{
Merge significant pathways that are with same genes over certain ratio
}
