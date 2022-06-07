\encoding{UTF-8}
\name{Rioja.data}
\docType{data}
\alias{Rioja.data}
\title{Inventoried Plots Data for a Stand Case Study in La Rioja}
\description{This list includes trees detected with TLS for 16 single scans corresponding to plots located in La Rioja, a region of Spain, in the north of the Iberian Peninsula (first element), as well as those inventoried in the field for these 16 plots (second element). Plot attributes related to stand stratum are also included (third element).

The elements of the list are as follows:
\enumerate{
    \item \code{tree.tls}: data frame that includes the list of trees detected with \code{\link{tree.detection.single.scan}} for 16 TLS single-scan sampling points. The following
variables are provided for each tree (see \code{\link{tree.detection.single.scan}} \sQuote{Value} for more details):

\tabular{rlll}{
[,1] \tab id   \tab character/numeric\cr
[,2] \tab file   \tab character\cr
[,3] \tab tree    \tab numeric\cr
[,4] \tab x    \tab numeric\cr
[,5] \tab y   \tab numeric\cr
[,6:8] \tab phi, phi.left, phi.right     \tab numeric\cr
[,9] \tab h.dist   \tab numeric\cr
[,10] \tab dbh   \tab numeric\cr
[,11] \tab h   \tab numeric\cr
[,12] \tab v   \tab numeric\cr
[,13:16] \tab n.pts, n.pts.red, n.pts.est, n.pts.red.est   \tab numeric\cr
[,15] \tab partial.occlusion   \tab numeric\cr}

\item \code{tree.list.field}: data frame that includes the list of trees measured in 16 circular fixed area plots of radius 20 m, whose centres coincide with TLS single-scans points. The following variables are provided for each tree:

\tabular{rlll}{
[,1] \tab id   \tab numeric \tab plot identification (coincident to TLS scans)\cr
[,2] \tab tree    \tab numeric \tab trees numbering\cr
[,3] \tab Sp    \tab numeric \tab specie code according to NFI ()\cr
[,4] \tab x    \tab numeric \tab x cartesian coordinate\cr
[,5] \tab y   \tab numeric \tab 7 cartesian coordinate\cr
[,6] \tab h.dist   \tab numeric \tab horizontal distance (m) from plot center to tree center\cr
[,7] \tab dbh   \tab numeric \tab tree diameter (cm) at breast height (1.3 m)\cr
[,8] \tab h   \tab numeric \tab tree total height (m)\cr
[,9] \tab dead   \tab numeric \tab dead (1) or not (NA)\cr
[,10] \tab v.user    \tab numeric \tab stem volume (m^3) estimated with allometric equations\cr
[,11] \tab w.user   \tab numeric \tab stem biomass (Mg) estimated with allometric equations\cr}}}
\usage{data(Rioja.data)}
\format{List with 2 data frames containing 604 observations and 17 variables (tree.tls) and 659 observations and 11 variables (tree.field).}
\keyword{datasets}
