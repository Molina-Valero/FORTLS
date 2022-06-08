\encoding{UTF-8}

\name{estimation.plot.size}

\alias{estimation.plot.size}

\title{Assess Consistency of Metrics for Simulated TLS Plots}

\description{Plots empirical linear charts of density (\emph{N}, trees/ha) and basal area (\emph{G}, \eqn{{m}^{2}/ha}{m^2/ha}) estimates (derived from simulated TLS plots) as a function of plot size (estimation-size charts) for different plot designs (circular fixed area, k-tree and angle-count), through continuous size increments (radius, k and BAF respectively). Size increments are set at 0.1 m, 1 tree and 0.1 \eqn{{m}^{2}/ha}{m^2/ha} for fixed area, k-tree and angle-count plot designs, respectively. These size-estimation line charts represent the consistency in predicting the stand variables across different values of radius, k and BAF. Size-estimation charts can be drawn for individual sample plots (including all plots together in the same charts) or for mean values (global mean computed for all the sample plots, or for group means if different strata are considered). Finally,
different plot designs can be compared if specified in the arguments, producing one size-estimation chart per variable (\emph{N} and \emph{G}).}

\usage{
estimation.plot.size(tree.tls,
                     plot.parameters = data.frame(radius.max = 25,
                                                  k.max = 50,
                                                  BAF.max = 4),
                     dbh.min = 4,
                     average = FALSE, all.plot.designs = FALSE)}

\arguments{

  \item{tree.tls}{Data frame with information of trees detected from TLS point cloud data in the same format as \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}} \sQuote{Value}.}

  \item{plot.parameters}{Optional data frame containing parameters for circular fixed area, k-tree and angle-count plot designs. The parameters are as follows:}

  \itemize{

  \item{\code{radius.max}}: maximum plot radius (m) considered for circular fixed area plots. If the \code{radius.max} specified is larger than the farthest tree from the plot centre, the horizontal distance from the farthest tree will be considered the maximum \code{radius}. By default, the \code{radius.max} will be 25 m.

\item{\code{k.tree.max}}: maximum number of trees considered for k-tree plots. If \code{k.tree.max} specified is larger than the maximum number of trees of the densest plot, this number of trees will be considered the maximum \code{k.tree.max}. By default, \code{k.tree.max} is 50.

\item{\code{BAF.max}}: maximum basal area factor (\eqn{{m}^{2}/ha}{m^2/ha}) considered for angle-count plots. By default, \code{BAF.max} is 4.}

  \item{dbh.min}{Optional minimum \emph{dbh} (cm) considered for detecting trees. By default it will be set at 4 cm.}

  \item{average}{Logical; if \code{TRUE}, plot means values and standard deviation of estimations will be represented. By default, it will be set as \code{FALSE}.}

  \item{all.plot.designs}{Logical; if \code{TRUE}, charts for each plot design are drawn together. By default, it will be set as \code{FALSE}.}}

\details{If there are strata in the \code{tree.tls} argument, they will be differentiated in charts with different colours. Strata must be specified in a numeric column named \code{stratum}.

The \code{all.plot.designs} argument only works for single strata, and therefore if there are additional strata in the \code{tree.tls} argument, they will be considered equal.

The outputs of this function are inspired by Fig. 3 of Brunner and Gizachew (2014).}

\value{\code{Invisible NULL}}

\references{
Brunner, A., & Gizachew, B. (2014). Rapid detection of stand density, tree positions, and tree diameter with a 2D terrestrial laser scanner. \emph{European Journal of Forest Research}, \bold{133(5)}, 819-831.}

\author{
Juan Alberto Molina-Valero,
María José Ginzo Villamayor,
Manuel Antonio Novo Pérez,
Adela Martínez-Calvo,
Juan Gabriel Álvarez-González,
Fernando Montes and
César Pérez-Cruzado.}

\note{Mean values are relevant when plots are representing homogenous strata.

Note that this is an option for choosing the best plot design when  field data are not available. Otherwise, using \code{\link{correlations}}, \code{\link{relative.bias}} and \code{\link{optimize.plot.design}} will be more desirable for obtaining the best possible plot design.}

\seealso{
\code{\link{tree.detection.single.scan}}, \code{\link{tree.detection.multi.scan}}, \code{\link{tree.detection.several.plots}}
}

\examples{
\donttest{

# Loading dataset with trees detected from TLS single-scans

data("Rioja.data")

tree.tls <- Rioja.data$tree.tls


# Without strata and plot parameters by default

estimation.plot.size(tree.tls)

estimation.plot.size(tree.tls, average = TRUE)

estimation.plot.size(tree.tls, all.plot.designs = TRUE)


}
}
