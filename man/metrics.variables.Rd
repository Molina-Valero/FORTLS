\encoding{UTF-8}

\name{metrics.variables}

\alias{metrics.variables}

\title{Compute Metrics and Variables for Terrestrial-Based Technologies Point Clouds}

\description{This function computes a set of metrics and variables from the Terrestrial-Based Technologies point cloud data, which have a high potential to be related or used as direct estimates (in the case of variables) of forest attributes at plot level. These can be obtained for different plot designs (circular fixed area, k-tree and angle-count plots). This function also includes methodologies for correcting occlusions generated in TLS single-scan point clouds.}

\usage{
metrics.variables(tree.tls, tree.ds = NULL, tree.field = NULL,
                  plot.design = c("fixed.area", "k.tree", "angle.count"),
                  plot.parameters = data.frame(radius = 25, k = 10, BAF = 2),
                  scan.approach = "single",
                  var.metr = list(tls = NULL, field = NULL),
                  v.calc = "parab", dbh.min = 4, h.min = 1.3,
                  max.dist = Inf, dir.data = NULL,
                  save.result = TRUE, dir.result = NULL)}

\arguments{

  \item{tree.tls}{Data frame with information about trees detected from terrestrial-based technologies point clouds data in the same format as \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multi.scan}} \sQuote{Values}.}

  \item{tree.ds}{Optional list containing detection probabilities of trees with distance sampling methods. The format must be the same as the \sQuote{Value} obtained with \code{\link{distance.sampling}}. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, TLS metrics using distance sampling based correction will not be calculated for circular fixed area and k-tree plot designs.}

    \item{tree.field}{
  Data frame with information about trees measured in the field plots. Each
  row must correspond to a (plot, tree) pair, and it must include at least the
  following columns:
    \itemize{
      \item \code{id}: plot identification encoded as character string or
      numeric. Plot identifications must coincide with those
      included in \code{id} column of the \code{tree.tls} argument.
      \item \code{tree}: trees numbering.
      \item \code{h.dist}: horizontal distance (m) from plot
      centre to tree centre. Centres of the field plots must coincide with
      centres of their corresponding Terrestrial-based technologies plots.
      \item \code{dbh}: tree diameter (cm) at breast height (1.3 m).
      \item \code{h}: tree total height (m).
      \item \code{h.blc}: height based live crown (m) (optional).
      \item \code{h.bdc}: height based dead crown (m) (optional).
      \item \code{v.user}: tree volume (m^3) (optional).
      \item \code{w.user}: tree biomass (Mg) (optional).
      \item \code{dead}: integer value indicating for each tree if it is dead
      (1) or not (NA).
    }
  }

  \item{plot.design}{Vector containing the plot designs considered. By default, all plot designs will be considered (circular fixed area, k-tree and angle-count plots).}

  \item{plot.parameters}{Data frame containing parameters for circular fixed area, k-tree and angle-count plot designs. If there is a \code{stratum} column in the \code{tree.list.tls} argument, it must have the same number of rows as strata values and they must be named using strata encoding. If plot parameters are not specified, the corresponding plot designs will not be considered in the function. If no parameter is specified, the function will stop giving an error message! The parameters are as follows:}

  \itemize{

  \item{\code{radius}}: plot radius (m) considered for circular fixed area plots. Absence of this argument rules out this plot design.

\item{\code{k.tree}}: number of trees (trees) considered for k-tree plots. Absence of this argument rules out this plot design.

\item{\code{BAF}}: basal area factor (\eqn{{m}^{2}/ha}{m^2/ha}) considered for angle-count plots. Absence of this argument rules out this plot design.

\item{\code{num.trees}}: number of dominant trees per ha (tree/ha), i.e. those with largest \emph{dbh}, considered for calculating dominant diameters and heights. In the absence of this argument, the number will be assumed to be 100 trees/ha.}

  \item{scan.approach}{Character parameter indicating TLS single-scan (\sQuote{single}) or TLS multi-scan approach or SLAM point clouds (\sQuote{multi}) approaches. If this argument is not specified by the user, it will be set to \sQuote{multi} approach.}

  \item{var.metr}{Optional vector containing all the metrics and variables of interest. By default it will be set as NULL and thus, all the metrics and variables available will be generated.}

  \item{v.calc}{Optional parameter to calculate volume when is not included in tree.tls input data.}

  \item{dbh.min}{Optional minimum \emph{dbh} (cm) considered for detecting trees. By default it will be set at 4 cm.}

  \item{h.min}{Optional minimum \emph{h} (m) considered for detecting trees. By default it will be set at 1.3 m.}

  \item{max.dist}{Optional argument to specify the maximum horizontal distance considered in which trees will be included.}

  \item{dir.data}{Optional character string naming the absolute path of the directory where TXT
  files containing TLS point clouds are located. \code{.Platform$file.sep} must
  be used as the path separator in \code{dir.data}, and TXT files in the
  directory must have the same description and format as indicated for TXT
  files in \code{\link{normalize}} \sQuote{Output Files}. If this
  argument is not specified by the user, it will be set to \code{NULL} by default
  and, as consequence, the current working directory of the \R process will be
  assigned to \code{dir.dat} during the execution.}

  \item{save.result}{Optional logical which indicates whether or not the output files described in \sQuote{Output Files} section must be saved in \code{dir.result}. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as consequence, the output files will be saved.}

  \item{dir.result}{Optional character string naming the absolute path of an existing directory where files described in \sQuote{Output Files} section will be saved. \code{.Platform$file.sep} must be used as the path separator in \code{dir.result}. If this argument is not specified by the user, and \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.result} during the execution.}}

\details{This function also works for several plots. In this case, a column named "id" to identify plots (string character or numeric) must be included in the \code{tree.list.tls} database argument. This must coincide with the corresponding "id" assigned in \code{\link{normalize}} to TXT files saved in \code{dir.data} (for more details see \code{\link{normalize}}). In addition, if there are several strata, they can be processed separately according to \code{plot.parameters} values (where each row will represent one stratum). If \code{tree.list.tls} does not include a specific "stratum" column, it will be assumed to have only one stratum, which will be encoded according to \code{rownames(plot.parameters)[1]}.

Using TLS data, this function computes metrics and estimations of variables (see \sQuote{Value}) for plot design specified in the \code{plot.parameters} argument. The notation used for variables is based on IUFRO (1959).

At this stage, three plot designs are available:

\itemize{
    \item Circular fixed area plots, simulated only if a \code{radius} value is
    specified in the \code{plot.parameters} argument.
    \item k-tree plots, simulated only if a \code{k.tree} value is
    specified in the \code{plot.parameters} argument.
    \item Angle-count plots, simulated only if a \code{BAF} value is
    specified in the \code{plot.parameters} argument.
  }

Volume is estimated modelling stem profile as a paraboloid and calculating the volumes of revolution; where trees \emph{dbh} are estimated in \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}}, and total heights are estimated as percentile 99 of z coordinate of points delimited by Voronoi polygons.

Regarding occlusion corrections for estimating variables, apart from distance sampling methods considered in \code{\link{distance.sampling}}, another occlusion corection based on correcting the shadowing effect (Seidel & Ammer, 2014) has been included to estimate some variables in circular fixed area and k-tree plots. In the case of angle-count plots, occlusion corrections are based on gap probability attenuation with distance to TLS depending on Poisson model (Strahler et al., 2008; Lovell et al., 2011).}


\value{List including TLS-based metrics and variables computed for plot designs considered. The list will contain one element per plot design considered (fixed.area.plot, k.tree.plot and angle.count.plot):

\item{fixed.area.plot}{If no value for \code{plot.radius} is specified in the \code{plot.parameters} argument, \code{NULL}; otherwise, the data frame will include TLS metrics and variables estimated in a circular fixed area plot design (rows represent simulations). The data frame will have the following columns:

Stratum, plot identification and radius:
    \itemize{
    \item \code{stratum}: stratum identification encoded as a character string or
      numeric. It must coincide with those included in the \code{stratum} column of \code{tree.list.tls}.
      \item \code{id}: plot identification encoded as character string or
      numeric. It will coincide with those included in the \code{id}
      column of \code{tree.list.tls} and, if applicable, the
      \code{distance.sampling} argument.
      \item \code{radius}: radius (m) of the simulated plot.}

      Terrestrial-Based Technologies variables:

      \itemize{

      \item \code{N.tls}: stand density (trees/ha) without occlusion
      corrections.
      \item \code{N.hn}, \code{N.hr}, \code{N.hn.cov}, \code{N.hr.cov}: stand
      density (trees/ha) with occlusion corrections based on distance sampling
      methodologies (see \code{\link{distance.sampling}}). These columns will not appear if the \code{distance.sampling}
      argument is \code{NULL}.
      \item \code{N.sh}: stand density (trees/ha) with correction of the
      shadowing effect.
      respectively).
      \item \code{G.tls}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}) without occlusion
      corrections.
      \item \code{G.hn}, \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov}: stand
      basal area (\eqn{{m}^{2}/ha}{m^2/ha}) with occlusion corrections based on distance sampling
      methodologies (see \code{\link{distance.sampling}}). These columns are missing if \code{distance.sampling}
      argument is \code{NULL}.
      \item \code{G.sh}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}) with correction of the
      shadowing effect.
      \item \code{V.tls}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha}) without occlusion corrections.
      \item \code{V.hn}, \code{V.hr}, \code{V.hn.cov}, \code{V.hr.cov}: stand
      volume (\eqn{{m}^{3}/ha}{m^3/ha}) with occlusion corrections based on distance sampling
      methodologies (see \code{\link{distance.sampling}}). These  columns will be missing if the \code{distance.sampling}
      argument is \code{NULL}.
      \item \code{V.sh}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha}) with correction of the
      shadowing effect.
      \item \code{d.tls}, \code{dg.tls}, \code{dgeom.tls}, \code{dharm.tls}:
      mean tree diameters (cm) at breast height (1.3 m), calculated from the arithmetic mean,
      quadratic mean, geometric mean and harmonic mean, respectively.
      \item \code{h.tls}, \code{hg.tls}, \code{hgeom.tls}, \code{hharm.tls}:
      mean tree heights (m), calculated from the arithmetic mean, quadratic mean, geometric
      mean and harmonic mean, respectively.
      \item \code{d.0.tls}, \code{dg.0.tls}, \code{dgeom.0.tls},
      \code{dharm.0.tls}: dominant mean tree diameters (cm) at breast height
      (1.3 m), calculated from the arithmetic mean, quadratic mean, geometric mean and
      harmonic mean, respectively.
      \item \code{h.0.tls}, \code{hg.0.tls}, \code{hgeom.0.tls},
      \code{hharm.0.tls}: dominant mean tree heights (m), calculated from the arithmetic
      mean, quadratic mean, geometric mean and harmonic mean, respectively.}

      Terrestrial-Based Technologies metrics:

      \itemize{

      \item \code{n.pts}, \code{n.pts.est}, \code{n.pts.red}, \code{n.pts.red.est}: number of points and estimated number of
      points (points) belonging to trees normal sections (+/- 5 cm) in the
      original point cloud (\code{n.pts} and \code{n.pts.est},
      respectively); and number of points and estimated number of points
      (points) belonging to trees normal sections (+/- 5 cm) in the reduced
      point cloud (\code{n.pts.red} and \code{n.pts.red.est},
      respectively).
      \item \code{P01}, \code{P05}, \code{P10}, \code{P20}, \code{P25},
      \code{P30}, \code{P40}, \code{P50}, \code{P60}, \code{P70}, \code{P75},
      \code{P80}, \code{P90}, \code{P95}, \code{P99}: height percentiles
      derived from z coordinates of TLS point clouds relative to ground level (m).
      \item \code{mean.arit.z/rho/phi}, \code{mean.qua.z/rho/phi}, \code{mean.geom.z/rho/phi}, \code{mean.harm.z/rho/phi}, \code{median.z/rho/phi}, \code{mode.z/rho/phi}: central tendency statistics for z, rho (horizontal distance) and r (euclidean distance) coordinates (arithmetic, quadratic, geometrical and harmonic means, median and mode, respectively).
      \item \code{var.z/rho/phi}, \code{sd.z/rho/phi}, \code{CV.z/rho/phi}, \code{D.z/rho/phi}, \code{ID.z/rho/phi}, \code{max.z/rho/phi}, \code{min.z/rho/phi}: dispersion tendency statistics for z, rho and r coordinates (variance, standar deviation, coefficient of variation, range, interquartile range, maximum and minimum, respectively).
      \item \code{kurt.z/rho/phi}, \code{skw.z/rho/phi}: curtosis and skewness, respectively.
      \item \code{L1.z/rho/phi}, \code{L2.z/rho/phi}, \code{L3.z/rho/phi}, \code{L4.z/rho/phi}, \code{L-CV.z/rho/phi}, \code{MAD-median.z/rho/phi}, \code{MAD-mode.z/rho/phi}: L-moments of order 2, 3 and 4, median of the absolute deviations from the overall median and mode of the absolute deviations from the overall mode, respectively.
      \item \code{L3.mu.z/rho/phi}, \code{L4.mu.z/rho/phi}: third and fourth central moments, respectively.
      \item \code{PA.2m}, \code{PA.mean.z}, \code{PA.mode.z}: percentage of laser returns above 2 m, mean and mode, respectively.
      \item \code{PB.2m}, \code{PB.mean.z}, \code{PB.median.z}, \code{PB.mode.z}: percentage of laser returns below 2 m, mean and mode, respectively.
      \item \code{weibull.b.z/rho/phi}, \code{weibull.c.z/rho/phi}: scale and shape parameters, respectively, for Weibull distribution fitted for z coordinates of TLS point clouds relative to ground level.}}


\item{k.tree.plot}{If no value for \code{k.tree} is specified in the \code{plot.parameters} argument, \code{NULL}; otherwise, the data frame will include TLS metrics and variables estimated in the k-tree plot design (rows represent simulations). The data frame will include the following columns (same description and format as indicated in \code{fixed.area.plot} element):

Stratum, plot identification and k:
    \itemize{
    \item \code{stratum}: stratum identification encoded as character string or
      numeric. It must coincide with those included in the \code{stratum} column of \code{tree.list.tls}.
      \item \code{id}: plot identification encoded as character string or numeric value. It will coincide with those included in \code{id}
      column of \code{tree.list.tls} or, if applicable, the
      \code{distance.sampling} argument.
      \item \code{k}: number of trees in the simulated plot.}

      Terrestrial-Based Technologies variables:
    \itemize{
      \item \code{N.tls}, \code{N.hn}, \code{N.hr}, \code{N.hn.cov},
      \code{N.hr.cov}, \code{N.sh}, \cr
      \code{G.tls}, \code{G.hn}, \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov},
      \code{G.sh}, \cr \code{V.tls}, \code{V.hn}, \code{V.hr}, \code{V.hn.cov},
      \code{V.hr.cov}, \code{V.sh}, \cr \code{d.tls}, \code{dg.tls},
      \code{dgeom.tls}, \code{dharm.tls}, \cr \code{h.tls}, \code{hg.tls},
      \code{hgeom.tls}, \code{hharm.tls}, \cr \code{d.0.tls}, \code{dg.0.tls},
      \code{dgeom.0.tls}, \code{dharm.0.tls}, \cr \code{h.0.tls}, \code{hg.0.tls},
      \code{hgeom.0.tls}, \code{hharm.0.tls}}

    Terrestrial-Based Technologies metrics:
    \itemize{
      \item \code{num.points},
      \code{num.points.est}, \code{num.points.hom}, \code{num.points.hom.est}, \cr          \code{P01},\code{P05},
      \code{P10}, \code{P20}, \code{P25}, \code{P30}, \code{P40}, \code{P50},
      \code{P60}, \code{P70}, \code{P75}, \code{P80}, \code{P90}, \code{P95},
      \code{P99}: same description and format
as indicated in fixed.area.plot element.}}


\item{angle.count.plot}{If no value for \code{BAF} is specified in the \code{plot.parameters} argument, \code{NULL}; otherwise, the data frame will include TLS metrics and variables estimated in the angle-count plot design (rows represent simulations). The data frame will include the following columns:

Stratum, plot identification and BAF:
    \itemize{
    \item \code{stratum}: stratum identification encoded as character string or
      numeric. It must coincide with those included in \code{stratum} column of \code{tree.list.tls}.
      \item \code{id}: plot identification encoded as character string or numeric. It will coincide with those included in the \code{id} column of \code{tree.list.tls}.
      \item \code{BAF}: BAF (\eqn{{m}^{2}/ha}{m^2/ha}) of the simulated plot.}


Terrestrial-Based Technologies variables:
    \itemize{
      \item \code{N.tls}: same description and format as indicated in the
      \code{fixed.area.plot} element.
      \item \code{N.pam}: stand density (trees/ha) with occlusion correction
      based on gap probability attenuation with distance to TLS.
      \code{fixed.area.plot} element.
      \item \code{G.tls}: same description and format as indicated for the
      \code{fixed.area.plot} element.
      \item \code{G.pam}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}) with occlusion correction
      based on gap probability attenuation with distance from TLS.
      \item \code{V.tls}: same description and format as indicated for the
      \code{fixed.area.plot} element.
      \item \code{V.pam}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha}) with occlusion correction based
      on gap probability attenuation with distance from TLS.
      \item \code{d.tls}, \code{dg.tls}, \code{dgeom.tls}, \code{dharm.tls}, \cr
      \code{h.tls}, \code{hg.tls}, \code{hgeom.tls}, \code{hharm.tls}, \cr
      \code{d.0.tls}, \code{dg.0.tls}, \code{dgeom.0.tls}, \code{dharm.0.tls}, \cr
      \code{h.0.tls}, \code{hg.0.tls}, \code{hgeom.0.tls}, \code{hharm.0.tls}}

      Terrestrial-Based Technologies metrics:
    \itemize{
      \item \code{num.points}, \code{num.points.est}, \code{num.points.hom},
      \code{num.points.hom.est}, \cr
      \code{P01}, \code{P05}, \code{P10}, \code{P20}, \code{P25}, \code{P30},
      \code{P40}, \code{P50}, \code{P60}, \code{P70}, \code{P75}, \code{P80},
      \code{P90}, \code{P95}, \code{P99}: same description and format as
      indicated for the \code{fixed.area.plot} element.}}}



\section{Output Files}{

After computing metrics and variables, if the \code{save.result} argument is \code{TRUE}, the function will save the elements in the list described in \sQuote{Value} (\code{fixed.area.plot}, \code{k.tree.plot} and \code{angle.count.plot}), which are different from \code{NULL} as CSV files.  Data frames are written without row names in the \code{dir.result} directory by using \code{\link{write.csv}} function from the \pkg{utils} package. The pattern used for naming these files is ‘metrics.variables.<plot design>.csv’,
being ‘<plot design>’ equal to “fixed.area.plot”, “k.tree.plot” or “angle.count.plot”
according to the plot design.}



\note{

In order to optimize plot designs and, therefore, for better use of \code{\link{metrics.variables}}, other functions such as \code{\link{correlations}}, \code{\link{relative.bias}} and \code{\link{estimation.plot.size}} should be used.

This function will be updated as new metrics are developed.}


\references{IUFRO (1959). Standarization of symbols in forest mensuration. Vienna, Austria, IUFRO.

Lovell, J. L., Jupp, D. L. B., Newnham, G. J., & Culvenor, D. S. (2011). Measuring tree stem diameters using intensity profiles from ground-based scanning lidar from a fixed viewpoint. \emph{ISPRS Journal of Photogrammetry and Remote Sensing}, \bold{66(1)}, 46-55. \doi{10.1016/j.isprsjprs.2010.08.006}

Seidel, D., & Ammer, C. (2014). Efficient measurements of basal area in short rotation forests based on terrestrial laser scanning under special consideration of shadowing. \emph{iForest-Biogeosciences and Forestry}, \bold{7(4)}, 227. \doi{10.3832/ifor1084-007}

Strahler, A. H., Jupp, D. L. B., Woodcock, C. E., Schaaf, C. B., Yao, T., Zhao,
F., Yang, X., Lovell, J., Culvenor, D., Newnham, G., Ni-Miester, W., &
Boykin-Morris, W. (2008). Retrieval of forest structural parameters using a
ground-based lidar instrument (Echidna®). \emph{Canadian Journal of Remote
Sensing}, \bold{34(sup2)}, S426-S440. \doi{10.5589/m08-046}}

\author{
Juan Alberto Molina-Valero,
Adela Martínez-Calvo.}

\seealso{
\code{\link{tree.detection.single.scan}}, \code{\link{tree.detection.multi.scan}}, \code{\link{tree.detection.several.plots}}, \code{\link{distance.sampling}},
\code{\link{normalize}}.
}

\examples{

\donttest{

# Establishment of working directories (optional)
# By default here we propose the current working directory of the R process

dir.data <- tempdir()
dir.result <- tempdir()


# Loading example data included in FORTLS

data("Rioja.data")
tree.tls <- Rioja.data$tree.tls
tree.tls <- tree.tls[tree.tls$id == "1", ]

# Download example of TXT file corresponding to plot 1 from Rioja data set

download.file(url = "https://www.dropbox.com/s/w4fgcyezr2olj9m/Rioja_1.txt?dl=1",
              destfile = file.path(dir.data, "1.txt"), mode = "wb")


# Considering distance sampling methods (only for single-scan approaches)

# ds <- distance.sampling(tree.tls)

met.var.TLS <- metrics.variables(tree.tls = tree.tls,
                                 # tree.ds = ds,
                                 plot.parameters = data.frame(radius = 10, k = 10, BAF = 2),
                                 dir.data = dir.data, dir.result = dir.result)


  }

}


