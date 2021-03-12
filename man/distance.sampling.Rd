\encoding{UTF-8}

\name{distance.sampling}

\alias{distance.sampling}

\title{Distance Sampling Methods for Correcting Occlusions Effects}

\description{Calculation of the probability of detection of every tree by using distance sampling methodologies (more specifically point transects methods), by means of fitting detection functions to the histogram of tree distribution according to their distance to TLS. Use both half normal and hazard rate functions, without and with \emph{dbh} as covariate. These probabilities are used for correcting estimation bias caused by lack of detection of trees due to occlusion.}

\usage{
distance.sampling(tree.list.tls,
                  id.plots = NULL,
                  strata.attributes = NULL)}

\arguments{

  \item{tree.list.tls}{Data frame with a list of trees detected and their \emph{dbh} and horizontal distances from TLS with the same structure and format as \code{\link{tree.detection}} and \code{\link{tree.detection.multiple}} \sQuote{Value}.}

  \item{id.plots}{Optional vector with plot identification encoded as character string or numeric for the plots considered. In this case, \code{tree.list.tls} argument must include a common column named \sQuote{id}. If this argument is not specified by the user, it will be set to NULL by default, and as a consequence, all plots will be considered.}

  \item{strata.attributes}{Optional data frame inluding plot radius considered at strata level. It must contain a column named \sQuote{stratum} (numeric) with encoding coinciding with that used in previous functions (\code{\link{normalize}}, \code{\link{tree.detection}} and \code{\link{tree.detection.multiple}}) for identifying strata. Therefore, strata must heve been included previously in \sQuote{tree.list.tls}. Another column named \sQuote{plot.radius} (numeric) will be required to set maximum horizontal distance (m) considered for fitting detection probability functions. If this argument is not specified by the user, it will be set to NULL by default, and as a consequence, all trees will be included.}}

\details{All internal functions related to distance sampling methodologies are fitted with the \code{\link[Distance]{ds}} function included in the \CRANpkg{Distance} package.

Detection functions are left-truncated at 1 m, according to Astrup et al., (2014).

Same warning messages as \code{\link[Distance]{ds}} function are provided when fits do not converge or another warnings occur.

For further details on these point transects methods  and similar sampling methodologies, as well as their application with R, see Buckland et al., (2001); Marques & Buckland, (2003); Miller & Thomas, (2015) and Clark (2016). Examples of distance sampling analyses, as well as lectures, are available at \url{http://examples.distancesampling.org/} and \url{http://workshops.distancesampling.org/}.}

\value{List containing the following elements:

\item{tree}{Data frame with detection probabilities for every tree and method.}

\itemize{

  \item {\code{stratum}}: stratum identification (coincident with strata of \code{tree.list.tls}). If there are not strata, it will be set as a single stratum encoded as 1 (numeric).

  \item {\code{id}}: plot identification (coincident with id of \code{tree.list.tls}).

  \item {\code{tree}}: tree numbering (coincident with \code{tree} of \code{tree.list.tls}).

  \item {\code{P.hn}}: tree detection probability according to half normal function.

  \item {\code{P.hn.cov}}: tree detection probability according to half normal function with \emph{dbh} as covariate.

  \item {\code{P.hr}}: tree detection probability according to half rate function.

  \item {\code{P.hr.cov}}: tree detection probability according to half rate function with \emph{dbh} as covariate.}

\item{parameters}{Data frame with parameters estimated for detection functions (see references for understanding their meaning).}

\itemize{

  \item {\code{P.hn.scale}}: scale parameter for half normal function (sigma).

  \item {\code{P.hn.cov.scale.intercept}}: alpha.0 parameter of scale parameter for half normal function with \emph{dbh} as covariate.

  \item {\code{P.hn.cov.dbh}}: alpha.1 parameter of scale parameter for half normal function with \emph{dbh} as covariate.

  \item {\code{P.hr.scale}}: scale parameter for half rate function (sigma).

  \item {\code{P.hr.shape}}: shape parameter for half rate function (b).

  \item {\code{P.hr.cov.scale.intercept}}: alpha.0 parameter of scale parameter for half normal function with \emph{dbh} as covariate.

  \item {\code{P.hr.cov.dbh}}: alpha.1 parameter of scale parameter for half normal function with \emph{dbh} as covariate.

  \item {\code{P.hr.cov.shape}}: shape parameter for half rate function with \emph{dbh} as covariate (b).}


\item{AIC}{Data frame with Akaike information criterions (AIC) of every detection function fit.}

\itemize{

  \item {\code{P.hn}}: AIC of half normal function fit.

  \item {\code{P.hn.cov}}:  AIC of half normal function with \emph{dbh} as covariate fit.

  \item {\code{P.hr}}:  AIC of half rate function fit.

  \item {\code{P.hr.cov}}:  AIC of half rate function with \emph{dbh} as covariate fit.}}

  \note{Although this step is optional for other functionalities of \pkg{FORTLS}, such as obtaining metrics and assessing the best plot designs (implemented in \code{\link{metrics.variables}}, \code{\link{correlations}}, \code{\link{relative.bias}} and \code{\link{optimize.plot.design}}), its inclusion is highly recommended, especially with high rates of occlusions.

  Note that this function could be more useful after assessing the best possible plot design with \code{\link{estimation.plot.size}}, \code{\link{correlations}}, \code{\link{relative.bias}} or \code{\link{optimize.plot.design}} functions.}

\references{
Astrup, R., Ducey, M. J., Granhus, A., Ritter, T., & von Lüpke, N. (2014). Approaches for estimating stand-level volume using terrestrial laser scanning in a single-scan mode. \emph{Canadian Journal of Forest Research}, \bold{44(6)}, 666-676. \doi{https://doi.org/10.1139/cjfr-2013-0535}.

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers, D. L., & Thomas, L. (2001). \emph{Introduction to distance sampling: estimating abundance of biological populations}, Oxford, United Kindown, Oxford University Press.

Clark, R. G. (2016). Statistical efficiency in distance sampling. \emph{PloS one}, \bold{11(3)}, e0149298. \doi{https://doi.org/10.1371/journal.pone.0149298}.

Marques, F. F., & Buckland, S. T. (2003). Incorporating covariates into standard line transect analyses. \emph{Biometrics}, \bold{59(4)}, 924-935. \doi{https://doi.org/10.1111/j.0006-341X.2003.00107.x}.

Miller, D. L., & Thomas, L. (2015). Mixture models for distance sampling detection functions. \emph{PloS one}, \bold{10(3)}, e0118726. \doi{https://doi.org/10.1371/journal.pone.0118726}.}

\author{
Juan Alberto Molina-Valero,
María José Ginzo Villamayor,
Manuel Antonio Novo Pérez,
Adela Martínez-Calvo,
Juan Gabriel Álvarez-González,
Fernando Montes and
César Pérez-Cruzado.}

\seealso{
\code{\link{tree.detection}}, \code{\link{tree.detection.multiple}}, \code{\link{metrics.variables}}, \code{\link{simulations}}.

\code{\link{ds}} in \CRANpkg{Distance} package.
}


\examples{

\donttest{

# Loading example data

data(Rioja.data)

tree.list.tls <- Rioja.data$tree.list.tls

# Whithout considering maximum distance

ds <- distance.sampling(tree.list.tls)


# Considering only some plots (first 8 plots)

ds <- distance.sampling(tree.list.tls, id.plots = 1:8)


# Considering strata
# Loading dataset with strata

plot.attributes <- Rioja.data$plot.attributes

# Merging the plot.attributes data set with strata information

tree.list.tls <- merge(tree.list.tls, plot.attributes, by = "id")


# Considering maximum distances of 10 and 15 m for stratum 1 and 2 respectively

strata.attributes = data.frame(stratum = c(1, 2),
                               plot.radius = c(10, 15))

ds <- distance.sampling(tree.list.tls, strata.attributes = strata.attributes)


}}

