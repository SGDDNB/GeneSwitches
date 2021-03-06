% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_switch_pathway.R
\name{plot_pathway_density}
\alias{plot_pathway_density}
\title{Pathways ridge plot}
\usage{
plot_pathway_density(switch_pw_re, toplotgl_sig, pw_direction = c("up",
  "down"), orderbytime = TRUE)
}
\arguments{
\item{switch_pw_re}{significant pathways to plot}

\item{toplotgl_sig}{switching genes}

\item{orderbytime}{order the pathways by switching time (mean time of switching genes) if TRUE,
order the pathways by FDR if FALSE}

\item{direction}{switching direction of the pathway, up or down}
}
\value{

}
\description{
This function generates pathways ridge plots
}
