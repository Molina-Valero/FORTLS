\encoding{UTF-8}
\name{FORTLS-package}
\alias{FORTLS}
\alias{FORTLS-package}
\docType{package}
\title{FORTLS: Automatic processing of TLS point cloud data for forestry purposes}
\description{
\packageDescription{FORTLS}
}
\details{

Usage of \pkg{FORTLS} includes the following functionalities:

\itemize{

\item{Tree detection: this is the first and necessary step for the other functionalities of \pkg{FORTLS}. This can be achieved using the following functions:}

\enumerate{

\item{\code{\link{normalize}}: mandatory first step for obtaining the relative coordinates of a TLS point cloud.}
\item{\code{\link{tree.detection}}: detects as many trees as possible from a normalized TLS point cloud.}
\item{\code{\link{tree.detection.multiple}}: includes the two previous functions for a better workflow when there are several plots to be sequentially analyzed.}}

\item{Estimation of variables when no field data are available: this is the main functionality of \pkg{FORTLS} and can be achieved using the following functions:}

\enumerate{

\item{\code{\link{distance.sampling}}: optional function which can be used for considering methodologies for correcting occlusion effects in estimating variables.}

\item{\code{\link{estimation.plot.size}}: enables the best plot design to be determined on the basis of TLS data only.}

\item{\code{\link{metrics.variables}}: is used for estimating metrics and variables potentially related to forest attributes at stand level.}}

\item{Estimation of variables when field data are available: this is the main and most desirable functionality of \pkg{FORTLS} and can be achieved using the following functions:}

\enumerate{

\item{\code{\link{distance.sampling}}: as before.}

\item{\code{\link{simulations}}: computes simulations of TLS and field data for different plot designs. This is a prior step to the next functions.}

\item{\code{\link{relative.bias}}: uses \code{\link{simulations}} output to assess the accuracy of direct estimations of variables according to homologous TLS and field data.}

\item{\code{\link{correlations}}: uses \code{\link{simulations}} output to assess correlations among metrics and variables obtained from TLS data, and variables of interest estimated from field data.}

\item{\code{\link{optimize.plot.design}}: using \code{\link{correlations}} output, represents the best correlations for variables of interest according to the plot design. It is thus possible to select the best plot design for estimating forest attributes from TLS data. }

\item{\code{\link{metrics.variables}}: as before, but in this case plot parameters will be choosen on the basis of field data and better estimates will therefore be obtained.}}}}


\author{
\strong{Maintainer}: Juan Alberto Molina-Valero \email{juanalberto.molina.valero@usc.es} [copyright holder]

Authors:
\itemize{
  \item María José Ginzo Villamayor [contributor]
  \item Manuel Antonio Novo Pérez [contributor]
  \item Adela Martínez-Calvo [contributor]
  \item Juan Gabriel Álvarez-González [contributor]
  \item Fernando Montes [contributor]
  \item César Pérez-Cruzado [contributor]
}}

\references{

Molina-Valero J. A., Ginzo-Villamayor M. J., Novo Pérez M. A., Martínez-Calvo A., Álvarez-González J. G., Montes F., & Pérez-Cruzado C. (2019). FORTLS: an R package for processing TLS data and estimating stand variables in forest inventories. \emph{The 1st International Electronic Conference on Forests — Forests for a Better Future: Sustainability, Innovation, Interdisciplinarity}. \doi{https://doi.org/10.3390/IECF2020-08066}}
