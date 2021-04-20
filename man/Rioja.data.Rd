\encoding{UTF-8}
\name{Rioja.data}
\docType{data}
\alias{Rioja.data}
\title{Inventoried Plots Data for a Stand Case Study in La Rioja}
\description{This list includes trees detected with TLS for 16 single scans corresponding to plots located in La Rioja, a region of Spain, in the north of the Iberian Peninsula (first element), as well as those inventoried in the field for these 16 plots (second element). Plot attributes related to stand stratum are also included (third element).

The elements of the list are as follows:
\enumerate{
    \item \code{tree.list.tls}: data frame that includes the list of trees detected with \code{\link{tree.detection}} for 16 TLS single-scans. The following
variables are provided for each tree (see \code{\link{tree.detection}} \sQuote{Value} for more details):

\tabular{rlll}{
[,1] \tab id   \tab character/numeric\cr
[,2] \tab file   \tab character\cr
[,3] \tab tree    \tab numeric\cr
[,4] \tab x    \tab numeric\cr
[,5] \tab y   \tab numeric\cr
[,6:8] \tab phi, phi.left, phi.right     \tab numeric\cr
[,9] \tab horizontal.distance   \tab numeric\cr
[,10] \tab dbh   \tab numeric\cr
[,11:14] \tab num.points, num.points.hom, num.points.est, num.points.hom.est   \tab numeric\cr
[,15] \tab partial.occlusion   \tab numeric\cr}

\item \code{tree.list.field}: data frame that includes the list of trees measured in 16 circular fixed area plots of radius 20 m, whose centres coincide with TLS single-scans points. The following variables are provided for each tree:

\tabular{rlll}{
[,1] \tab id   \tab numeric \tab plot identification (coincident to TLS scans)\cr
[,2] \tab tree    \tab numeric \tab trees numbering\cr
[,3] \tab Sp    \tab numeric \tab specie code according to NFI ()\cr
[,4] \tab horizontal.distance   \tab numeric \tab horizontal distance (m) from plot center to tree center\cr
[,5] \tab dbh   \tab numeric \tab tree diameter (cm) at breast height (1.3 m)\cr
[,6] \tab total.height   \tab numeric \tab tree total height (m)\cr
[,7] \tab dead   \tab numeric \tab dead (1) or not (NA)\cr
[,8] \tab x    \tab numeric\cr
[,9] \tab y   \tab numeric\cr}

\item \code{plot.attributes}: dataframe that includes the stratum of each plot:

\tabular{rlll}{
      [,1] \tab id   \tab character/numeric \tab plot identification (coincident to TLS scans)\cr
      [,2] \tab stratum    \tab numeric \tab stratum (1-2)\cr}}}
\usage{data(Rioja.data)}
\format{List with 3 data frames containing 598 observations and 15 variables (tree.list.tls), 659 observations and 9 variables (tree.list.field), and 16 observations and 2 variables (plot.attributes).}
\keyword{datasets}
