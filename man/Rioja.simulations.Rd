\encoding{UTF-8}
\name{Rioja.simulations}
\docType{data}
\alias{Rioja.simulations}
\title{Simulated Metrics and Variables for a Stand Case Study in La Rioja}
\description{This list contains metrics and variables estimated from field data and TLS data from \code{\link{Rioja.data}}.

The elements on this list correspond to \code{\link{simulations}} \sQuote{Value}, as follows:

\enumerate{
    \item \code{simulations.fixed.area.plot}: data frame with TLS metrics and variables estimated on the basis of simulated plots in a fixed area plot design with radius increment of 0.1 m (from smallest possible radius to 20 m). The following variables are provided for each pair (plot, radius) (see \code{\link{simulations}} \sQuote{Value} for more details):
    \tabular{rlll}{
      [,1] \tab id   \tab character/numeric\cr
      [,2] \tab radius    \tab numeric\cr
      [,3:5] \tab N, G, V    \tab numeric\cr
      [,6:9] \tab d, dg, dgeom, dharm   \tab numeric\cr
      [,10:13] \tab h, hg, hgeom, hharm   \tab numeric\cr
      [,14:17] \tab d.0, dg.0, dgeom.0,  dharm.0  \tab numeric\cr
      [,18:21] \tab h.0, hg.0, hgeom.0, hharm.0   \tab numeric\cr
      [,22:27] \tab N.tls, N.hn.tls, N.hr.tls, N.hn.cov.tls, N.hr.cov.tls, N.sh.tls  \tab numeric\cr
      [,28:31] \tab num.points, num.points.est, num.points.hom, num.points.hom.est   \tab numeric\cr numeric\cr
      [,32:37] \tab G.tls, G.hn.tls, G.hr.tls, G.hn.cov.tls, G.hr.cov.tls, G.sh.tls  \tab numeric\cr
      [,38:43] \tab V.tls, V.hn.tls, V.hr.tls, V.hn.cov.tls, V.hr.cov.tls, V.sh.tls   \tab numeric\cr
      [,44:47] \tab d.tls, dg.tls, dgeom.tls, dharm.tls   \tab numeric \cr
      [,48:51] \tab h.tls, hg.tls, hgeom.tls, hharm.tls   \tab numeric\cr
      [,52:55] \tab d.0, dg.0, dgeom.0, dharm.0  \tab numeric\cr
      [,56:59] \tab h.0, hg.0, hgeom.0, hharm.0   \tab numeric\cr
      [,60:74] \tab P01, P05, P10, P20, P25, P30, P40, P50, P60, P70, P75, P80, P90, P95, P99   \tab numeric\cr}

 \item \code{simulations.k.tree.plot}: data frame with TLS metrics and variables estimated on the basis of simulated plots under k-tree plot design for incremental values of 1 tree (from 1 to largest number of trees in one plot). The following variables are provided for each pair (plot, k) (see \code{\link{simulations}} \sQuote{Value} for more details):
    \tabular{rlll}{
      [,1] \tab id   \tab character/numeric\cr
      [,2] \tab k    \tab numeric\cr
      [,3:5] \tab N, G, V    \tab numeric\cr
      [,6:9] \tab d, dg, dgeom, dharm   \tab numeric\cr
      [,10:13] \tab h, hg, hgeom, hharm   \tab numeric\cr
      [,14:17] \tab d.0, dg.0, dgeom.0,  dharm.0  \tab numeric\cr
      [,18:21] \tab h.0, hg.0, hgeom.0, hharm.0   \tab numeric\cr
      [,22:27] \tab N.tls, N.hn.tls, N.hr.tls, N.hn.cov.tls, N.hr.cov.tls, N.sh.tls  \tab numeric\cr
      [,28:31] \tab num.points, num.points.est, num.points.hom, num.points.hom.est   \tab numeric\cr numeric\cr
      [,32:37] \tab G.tls, G.hn.tls, G.hr.tls, G.hn.cov.tls, G.hr.cov.tls, G.sh.tls  \tab numeric\cr
      [,38:43] \tab V.tls, V.hn.tls, V.hr.tls, V.hn.cov.tls, V.hr.cov.tls, V.sh.tls   \tab numeric\cr
      [,44:47] \tab d.tls, dg.tls, dgeom.tls, dharm.tls   \tab numeric \cr
      [,48:51] \tab h.tls, hg.tls, hgeom.tls, hharm.tls   \tab numeric\cr
      [,52:55] \tab d.0, dg.0, dgeom.0, dharm.0  \tab numeric\cr
      [,56:59] \tab h.0, hg.0, hgeom.0, hharm.0   \tab numeric\cr
      [,60:74] \tab P01, P05, P10, P20, P25, P30, P40, P50, P60, P70, P75, P80, P90, P95, P99   \tab numeric\cr}

      \item \code{simulations.angle.count.plot}: data frame with TLS metrics and variables estimated on the basis of simulated
plots in an angle-count plot design. They plots are simulated for correlative angle-count plots and incremental values of 0.1 \eqn{{m}^{2}/ha}{m^2/ha} for BAF. The following variables are provided for each pair (plot, BAF) (see \code{\link{simulations}} \sQuote{Value} for more details):

    \tabular{rlll}{
      [,1] \tab id   \tab character/numeric\cr
      [,2] \tab BAF    \tab numeric\cr
      [,3:5] \tab N, G, V    \tab numeric\cr
      [,6:9] \tab d, dg, dgeom, dharm   \tab numeric\cr
      [,10:13] \tab h, hg, hgeom, hharm   \tab numeric\cr
      [,14:17] \tab d.0, dg.0, dgeom.0,  dharm.0  \tab numeric\cr
      [,18:21] \tab h.0, hg.0, hgeom.0, hharm.0   \tab numeric\cr
      [,22:23] \tab N.tls, N.pam.tls   \tab numeric\cr
      [,24:27] \tab num.points, num.points.est, num.points.hom, num.points.hom.est   \tab numeric\cr numeric\cr
      [,28:29] \tab G.tls, G.pam.tls   \tab numeric\cr
      [,30:31] \tab V.tls, V.pam.tls   \tab numeric\cr
      [,32:35] \tab d.tls, dg.tls, dgeom.tls, dharm.tls   \tab numeric \cr
      [,48:51] \tab h.tls, hg.tls, hgeom.tls, hharm.tls   \tab numeric\cr
      [,36:39] \tab d.0, dg.0, dgeom.0, dharm.0  \tab numeric\cr
      [,40:43] \tab h.0, hg.0, hgeom.0, hharm.0   \tab numeric\cr
      [,44:62] \tab P01, P05, P10, P20, P25, P30, P40, P50, P60, P70, P75, P80, P90, P95, P99   \tab numeric\cr}}}

\usage{data(Rioja.simulations)}
\format{List with 3 data frames containing 2224 observations and 74 variables (simulations.fixed.area.plot), 272 observations and 74 variables (simulations.k.tree.plot), and 576 observations and 62 variables (simulations.angle.count.plot).}
\keyword{datasets}
