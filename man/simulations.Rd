
\encoding{UTF-8}


\name{simulations}


\alias{simulations}


\title{
Compute Metrics and Variables for Simulated TLS and Field Plots
}


\description{
Computes TLS metrics derived from simulated TLS plots and variables estimated on the basis of simulated field plots. Real TLS and field data from the same
set of plots are required in order to build simulated plots. Three different
plot designs are currently available: circular fixed area, k-tree and angle-count.
During the simulation process, plots with incremental values for radius, k and
BAF are simulated for circular fixed area, k-tree and angle-count designs, respectively,
according to the parameters specified in the \code{plot.parameters} argument. For
TLS metrics, different method are included for correcting occlusions generated in
TLS point clouds.
}


\usage{
simulations(tree.list.tls, distance.sampling = NULL, tree.list.field,
            plot.parameters = list(radius.max = 25, k.tree.max = 50,
                                   BAF.max = 4),
            dir.data = NULL, save.result = TRUE, dir.result = NULL)
}


\arguments{

  \item{tree.list.tls}{
  Data frame with information about trees detected from TLS point cloud data. The
  structure and format must be analogous to output  returned by
  \code{\link{tree.detection}} and \code{\link{tree.detection.multiple}} functions. In particular, each row must
  correspond to a (plot, tree) pair, and it must include at least the following columns:
    \itemize{
      \item \code{id}, \code{file}, \code{tree}, \code{x}, \code{y},
      \code{phi.left}, \code{phi.right}, \code{horizontal.distance},
      \code{dbh}, \code{num.points}, \code{num.points.hom},
      \code{num.points.est}, \code{num.points.hom.est},
      \code{partial.occlusion}: same description and format as indicated for the
      same named columns in \code{\link{tree.detection}} \sQuote{Value}.
    }
  }

  \item{distance.sampling}{
  An optional list containing results arises from the application of distance
  sampling methodologies. The structure and format must be analogous to output
  returned by \code{\link{distance.sampling}} function. In particular, it must
  include at least the following elements:
    \itemize{
      \item \code{tree}: data frame with detection probabilities for each tree
      using distance sampling methodologies. Each row must correspond
      to a (plot, tree) pair, and it must include at least the following columns:
        \itemize{
          \item \code{id}, \code{tree}, \code{P.hn}, \code{P.hn.cov},
          \code{P.hr}, \code{P.hr.cov}: same description and format as
          indicated for same named columns of \code{tree} in \code{\link{distance.sampling}} \sQuote{Value}. In addition, plot
          identification and tree numbering included in \code{id} and
          \code{tree} columns must coincide with those included in the same
          named columns of \code{tree.list.tls} argument.
        }

  If this argument is not specified by the user, it will be set to \code{NULL} by
  default and, as a consequence, the TLS metrics using distance sampling based
  correction will not be calculated for a circular fixed area or k-tree plot designs.
}}

  \item{tree.list.field}{
  Data frame with information about trees measured in the field plots. Each
  row must correspond to a (plot, tree) pair, and it must include at least the
  following columns:
    \itemize{
      \item \code{id}: plot identification encoded as character string or
      numeric. Plot identifications must coincide with those
      included in \code{id} column of the \code{tree.list.tls} argument.
      \item \code{tree}: trees numbering.
      \item \code{horizontal.distance}: horizontal distance (m) from plot
      centre to tree centre. Centres of the field plots must coincide with
      centres of their corresponding TLS plots.
      \item \code{dbh}: tree diameter (cm) at breast height (1.3 m).
      \item \code{total.height}: tree total height (m).
      \item \code{dead}: integer value indicating for each tree if it is dead
      (1) or not (NA).
    }
  }

  \item{plot.parameters}{
  Optional list containing parameters for circular fixed area, k-tree and angle-count
  plot designs. User can set all or any of the following parameters specifying
  them as named elements of the list:
    \itemize{
      \item \code{radius.max}: maximum radius (m) allowed for the increasing
      sequence of radius values to be used in the simulation process of TLS and
      field plots under circular fixed area plot design. If this element is not included
      in the argument, circular fixed area plots will not be simulated.
      \item \code{radius.increment}: positive increment (m) for the increasing
      sequence of radius values to be used in the simulation process of TLS and
      field plots in a circular fixed area plot design. If this element is not included
      in the argument, it is set to 0.1 m by default.
      \item \code{k.tree.max}: maximum number of trees (trees) allowed for the
      increasing sequence of k values to be used in the simulation process of
      TLS and field plots under k-tree plot design. If this element is not
      included in the argument, k-tree plots will not be simulated.
      \item \code{BAF.max}: maximum BAF (\eqn{{m}^{2}/ha}{m^2/ha}) allowed for the increasing
      sequence of BAF values to be used in the simulation process of TLS and
      field plots in an angle-count plot design. If this element is not
      included in the argument, angle-count plots will not be simulated.
      \item \code{BAF.increment}: positive increment (\eqn{{m}^{2}/ha}{m^2/ha}) for the increasing
      sequence of BAF values to be used in the simulation process of TLS and
      field plots under angle-count plot design. If this element is not
      included in the argument, it is set to 0.1 \eqn{{m}^{2}/ha}{m^2/ha} by default.
      \item \code{num.trees}: number of dominant trees per ha (tree/ha) to be
      used for calculating dominant diameters and heights during the simulation
      process of TLS and field plots under all the available plot designs.
      Dominant trees are those with the largest diameter at breast height. If
      this element is not included in \code{plot.parameters} argument, it is
      set to 100 trees/ha by default.
    }
  If this argument is specified by the user, it must include at least one of
  the following elements: \code{radius.max}, \code{k.tree.max} or
  \code{BAF.max}. If this argument is not specified by the user, it is set to
  \code{list(radius.max = 25, k.tree.max = 50, BAF.max = 4)} by default and, as a
  consequence, the three available plot designs will be simulated.
  }

  \item{dir.data}{
  Optional character string naming the absolute path of the directory where TXT
  files containing TLS point clouds are located. \code{.Platform$file.sep} must
  be used as the path separator in \code{dir.dat}, and TXT files in the
  directory must have the same description and format as indicated for TXT
  files in \code{\link{normalize}} \sQuote{Output Files}. If this
  argument is not specified by the user, it will be set to \code{NULL} by default
  and, as a consequence, the current working directory of the \R process will be
  assigned to \code{dir.dat} during the execution.
  }

  \item{save.result}{
  Optional logical which indicates wheter or not the output files described in
  \sQuote{Output Files} section must be saved in \code{dir.result}. If
  this argument is not specified by the user, it will be set to \code{TRUE} by
  default and, as a consequence, the output files are saved.
  }

  \item{dir.result}{
  Optional character string naming the absolute path of an existing directory
  where files described in \sQuote{Output Files} section will be saved.
  \code{.Platform$file.sep} must be used as the path separator in
  \code{dir.result}. If this argument is not specified by the user, and
  \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and,
  as a consequence, the current working directory of the \R process will be
  assigned to \code{dir.result} during the execution.
  }
}


\details{
Using real TLS and field data from the same set of plots, this function enables construction
of simulated plots under different plot designs and computation of the
corresponding TLS metrics and estimated variables. The notation used for
variables is based on IUFRO (1959).

At this stage, three plot designs are available:
  \itemize{
    \item Circular fixed area plots, simulated only if a \code{radius.max} value is
    specified in the \code{plot.parameters} argument.
    \item k-tree plots, simulated only if a \code{k.tree.max} value is
    specified in the \code{plot.parameters} argument.
    \item Angle-count plots, simulated only if a \code{BAF.max} value is
    specified in the \code{plot.parameters} argument.
  }
For each real plot, a simulation process is run under each of the plot designs
specified by means of elements of the \code{plot.parameters} argument. Although there
are some minor differences depending on the plot design, the rough outline of
the simulation process is similar for all, and it consists of the
following main steps:
  \enumerate{
    \item Define an increasing sequence of the plot design parameter (radius,
    k or BAF) according to the maximum value and, if applicable, the positive
    increment set in \code{plot.parameters} argument.
    \item Build simulated plots for each parameter value in the previous
    sequence based on either TLS or field data.
    \item Compute either TLS metrics or variables estimated on the basis of
    simulated plots for each parameter value (see \sQuote{Value} section for
    details). For the simulated TLS plots, note that in
    addition to the counterparts of variables computed for the simulated field
    plots, the function also computes the following:
      \itemize{
        \item Metrics related to the number of points belonging to normal tree
        sections.
        \item Metrics with occlusion corrections based on the following:
        \itemize{
          \item Distance sampling methodologies (Astrup et al., 2014) for
          circular fixed area and k-tree plot designs, if the \code{distance.sampling}
          argument is not \code{NULL}.
          \item Correction of the shadowing effect (Seidel & Ammer, 2014) for
          circular fixed area and k-tree plot designs.
          \item Gap probability attenuation with distance to TLS (Strahler et
          al., 2008; Lovell et al., 2011) for angle-count plot design.
        }
        \item Height percentiles derived from z coordinates of TLS point clouds
         relative to ground level.
      }
  }
}


\value{List with field estimates and TLS metrics for plot designs considered. It will contain one element per plot design considered (fixed.area.plot, k.tree.plot and angle.count.plot)

  \item{fixed.area.plot}{
  If no value for \code{radius.max} is specified in the \code{plot.parameters}
  argument, \code{NULL}; otherwise, data frame with TLS metrics and variables will be estiamted on the basis of simulated plots in a circular fixed area plot design. Each row will correspond to a (plot, radius) pair, and the
  following columns will be included:

  Plot identification and radius:
    \itemize{
      \item \code{id}: plot identification encoded as a character string or
      numeric. It will coincide with those included in the \code{id}
      column of \code{tree.list.tls}, \code{tree.list.field} or, if applicable,
      \code{distance.sampling} arguments.
      \item \code{radius}: radius (m) of the simulated plot.
    }
  Variables estimated on the basis of simulated field plots:
    \itemize{
      \item \code{N}: stand density (trees/ha).
      \item \code{G}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}).
      \item \code{V}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha}).
      \item \code{d}, \code{dg}, \code{dgeom}, \code{dharm}: mean tree
      diameters (cm) at breast height (1.3 m), calculated from the arithmetic mean, quadratic
      mean, geometric mean and harmonic mean, respectively.
      \item \code{h}, \code{hg}, \code{hgeom}, \code{hharm}: mean tree
      heights (m), calculated from the arithmetic mean, quadratic mean, geometric mean and
      harmonic mean, respectively.
      \item \code{d.0}, \code{dg.0}, \code{dgeom.0}, \code{dharm.0}: dominant
      mean tree diameters (cm) at breast height (1.3 m), calcualted from the arithmetic mean,
      quadratic mean, geometric mean and harmonic mean, respectively.
      \item \code{h.0}, \code{hg.0}, \code{hgeom.0}, \code{hharm.0}: dominant
      mean tree heights (m), calculated from the arithmetic mean, quadratic mean, geometric
      mean and harmonic mean, respectively.
      }
  TLS metrics derived from simulated TLS plots:
    \itemize{
      \item \code{N.tls}: stand density (trees/ha) without occlusion
      corrections.
      \item \code{N.hn}, \code{N.hr}, \code{N.hn.cov}, \code{N.hr.cov}: stand
      density (trees/ha) with occlusion corrections based on distance sampling
      methodologies. These columns will be missing if the \code{distance.sampling}
      argument is \code{NULL}.
      \item \code{N.sh}: stand density (trees/ha) with correction of the
      shadowing effect.
      \item \code{num.points}, \code{num.points.est}, \code{num.points.hom},
      \code{num.points.hom.est}: number of points and estimated number of
      points (points) belonging to trees with normal sections (+/- 5 cm) in the
      original point cloud (\code{num.points} and \code{num.points.est},
      respectively); and number of points and estimated number of points
      (points) belonging to trees normal sections (+/- 5 cm) in the reduced
      point cloud (\code{num.points.hom} and \code{num.points.hom.est},
      respectively).
      \item \code{G.tls}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}) without occlusion
      corrections.
      \item \code{G.hn}, \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov}: stand
      basal area (\eqn{{m}^{2}/ha}{m^2/ha}) with occlusion corrections based on distance sampling
      methodologies. These columns will be missing if the \code{distance.sampling}
      argument is \code{NULL}.
      \item \code{G.sh}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}) with correction of the
      shadowing effect.
      \item \code{V.tls}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha}) without occlusion corrections.
      \item \code{V.hn}, \code{V.hr}, \code{V.hn.cov}, \code{V.hr.cov}: stand
      volume (\eqn{{m}^{3}/ha}{m^3/ha}) with occlusion corrections based on distance sampling
      methodologies. These  columns will be missing if the \code{distance.sampling}
      argument is \code{NULL}.
      \item \code{V.sh}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha}) with correction of the
      shadowing effect.
      \item \code{d.tls}, \code{dg.tls}, \code{dgeom.tls}, \code{dharm.tls}:
      mean tree diameters (cm) at breast height (1.3 m), calculated from the arithmetic mean,
      quadratic mean geometric mean, and harmonic mean, respectively.
      \item \code{h.tls}, \code{hg.tls}, \code{hgeom.tls}, \code{hharm.tls}:
      mean tree heights (m), calculated from the arithmetic mean, quadratic mean, geometric
      mean and harmonic mean, respectively.
      \item \code{d.0.tls}, \code{dg.0.tls}, \code{dgeom.0.tls},
      \code{dharm.0.tls}: dominant mean tree diameters (cm) at breast height
      (1.3 m), calculated from the arithmetic mean, quadratic mean, geometric mean and
      harmonic mean, respectively.
      \item \code{h.0.tls}, \code{hg.0.tls}, \code{hgeom.0.tls},
      \code{hharm.0.tls}: dominant mean tree heights (m), calculated from the arithmetic
      mean, quadratic mean, geometric mean and harmonic mean, respectively.
      \item \code{P01}, \code{P05}, \code{P10}, \code{P20}, \code{P25},
      \code{P30}, \code{P40}, \code{P50}, \code{P60}, \code{P70}, \code{P75},
      \code{P80}, \code{P90}, \code{P95}, \code{P99}: height percentiles
      derived from z coordinates of TLS point clouds relative to ground level.
    }
  }

  \item{k.tree.plot}{
  If no value for \code{k.tree.max} is specified in the \code{plot.parameters}
  argument, \code{NULL}; otherwise the data frame with TLS metrics and estimations
  of variables will be based on simulated plots in the k-tree plot design. Each of row will correspond to a (plot, k) pair, and the following
  columns will be included:

  Plot identification and k:
    \itemize{
      \item \code{id}: plot identification encoded as character string or
      numeric. The \code{id} will coincide with those included in the \code{id}
      column of \code{tree.list.tls}, \code{tree.list.field} or, if applicable,
      \code{distance.sampling} arguments.
      \item \code{k}: number of trees (trees) in the simulated plot.
    }
  Estimated variables based on simulated field plots:
    \itemize{
      \item \code{N}, \code{G}, \code{V}, \code{d}, \code{dg}, \code{dgeom},
      \code{dharm}, \code{h}, \code{hg}, \code{hgeom}, \code{hharm},
      \code{d.0}, \code{dg.0}, \code{dgeom.0}, \code{dharm.0}, \code{h.0},
      \code{hg.0}, \code{hgeom.0}, \code{hharm.0}: same description and format
      as indicated in the \code{fixed.area.plot} element.
    }
  TLS metrics derived from simulated TLS plots:
    \itemize{
      \item \code{N.tls}, \code{N.hn}, \code{N.hr}, \code{N.hn.cov},
      \code{N.hr.cov}, \code{N.sh}, \code{num.points}, \cr
      \code{num.points.est}, \code{num.points.hom}, \code{num.points.hom.est},
      \code{G.tls}, \code{G.hn}, \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov},
      \code{G.sh}, \code{V.tls}, \code{V.hn}, \code{V.hr}, \code{V.hn.cov}, \cr
      \code{V.hr.cov}, \code{V.sh}, \code{d.tls}, \code{dg.tls},
      \code{dgeom.tls}, \code{dharm.tls}, \code{h.tls}, \code{hg.tls},
      \code{hgeom.tls}, \code{hharm.tls}, \code{d.0.tls}, \code{dg.0.tls},
      \code{dgeom.0.tls}, \code{dharm.0.tls}, \code{h.0.tls}, \code{hg.0.tls},
      \code{hgeom.0.tls}, \code{hharm.0.tls}, \code{P01}, \code{P05},
      \code{P10}, \code{P20}, \code{P25}, \code{P30}, \code{P40}, \code{P50},
      \code{P60}, \code{P70}, \code{P75}, \code{P80}, \code{P90}, \code{P95},
      \code{P99}: same description and format as indicated in
      \code{fixed.area.plot} element.
    }
  }


  \item{angle.count.plot}{
  If no value for \code{BAF.max} is specified in the \code{plot.parameters}
  argument, \code{NULL}; otherwise the data frame will include TLS metrics and estimated
  variables based on simulated plots in the angle-count plot design. Each
  row will correspond to a (plot, BAF) pair, and the
  following columns will be included:

  Plot identification and BAF:
    \itemize{
      \item \code{id}: plot identification encoded as character string or
      numeric. The \code{id} will coincide with those included in the \code{id}
      column of \code{tree.list.tls} and \code{tree.list.field}.
      \item \code{BAF}: BAF (\eqn{{m}^{2}/ha}{m^2/ha}) of the simulated plot.
    }
  Estimated variables based on simulated field plots:
    \itemize{
      \item \code{N}, \code{G}, \code{V}, \code{d}, \code{dg}, \code{dgeom},
      \code{dharm}, \code{h}, \code{hg}, \code{hgeom}, \code{hharm},
      \code{d.0}, \code{dg.0}, \code{dgeom.0}, \code{dharm.0}, \code{h.0},
      \code{hg.0}, \code{hgeom.0}, \code{hharm.0}: same description and format
      as indicated in the \code{fixed.area.plot} element.
    }
  TLS metrics derived from simulated TLS plots:
    \itemize{
      \item \code{N.tls}: same description and format as indicated in the
      \code{fixed.area.plot} element.
      \item \code{N.pam}: stand density (trees/ha) with occlusion correction
      based on gap probability attenuation with distance to TLS.
      \item \code{num.points}, \code{num.points.est}, \code{num.points.hom},
      \code{num.points.hom.est}: same description and format as
      indicated in \code{fixed.area.plot} element.
      \item \code{G.tls}: same description and format as indicated in
      \code{fixed.area.plot} element.
      \item \code{G.pam}: stand basal area (\eqn{{m}^{2}/ha}{m^2/ha}) with occlusion correction
      based on gap probability attenuation with distance to TLS.
      \item \code{V.tls}: same description and format as indicated in
      \code{fixed.area.plot} element.
      \item \code{V.pam}: stand volume (\eqn{{m}^{3}/ha}{m^3/ha})with occlusion correction based
      on gap probability attenuation with distance to TLS.
      \item \code{d.tls}, \code{dg.tls}, \code{dgeom.tls}, \code{dharm.tls},
      \code{h.tls}, \code{hg.tls}, \code{hgeom.tls}, \cr \code{hharm.tls},
      \code{d.0.tls}, \code{dg.0.tls}, \code{dgeom.0.tls}, \code{dharm.0.tls},
      \code{h.0.tls}, \code{hg.0.tls}, \code{hgeom.0.tls}, \code{hharm.0.tls},
      \code{P01}, \code{P05}, \code{P10}, \code{P20}, \code{P25}, \code{P30},
      \code{P40}, \code{P50}, \code{P60}, \code{P70}, \code{P75}, \code{P80},
      \code{P90}, \code{P95}, \code{P99}: same description and format as
      indicated in \code{fixed.area.plot} element.
    }
  }
}


\section{Output Files}{
At the end of the simulation process, if the \code{save.result} argument is
\code{TRUE}, the function will print all the elements described in
\sQuote{Value} section and which are different from \code{NULL} to files. Data frames are
written without row names in \code{dir.result} directory using
the \code{\link{write.csv}} function from the \pkg{utils} package. The pattern used
for naming these files is \file{simulations.<plot design>.csv}, where
\file{<plot design>} is equal to \dQuote{\code{fixed.area.plot}},
\dQuote{\code{k.tree.plot}} or \dQuote{\code{angle.count.plot}} according to
plot design.
}


\references{
Astrup, R., Ducey, M. J., Granhus, A., Ritter, T., & von Lüpke, N. (2014).
Approaches for estimating stand level volume using terrestrial laser scanning
in a single-scan mode. \emph{Canadian Journal of Forest Research},
\bold{44(6)}, 666-676. \doi{10.1139/cjfr-2013-0535}

IUFRO (1959). \emph{Standarization of symbols in forest mensuration}. IUFRO,
Wien, 32 pp.

Lovell, J. L., Jupp, D. L. B., Newnham, G. J., & Culvenor, D. S. (2011).
Measuring tree stem diameters using intensity profiles from ground-based
scanning lidar from a fixed viewpoint. \emph{ISPRS Journal of Photogrammetry
and Remote Sensing}, \bold{66(1)}, 46-55. \doi{10.1016/j.isprsjprs.2010.08.006}

Seidel, D., & Ammer, C. (2014). Efficient measurements of basal area in short
rotation forests based on terrestrial laser scanning under special
consideration of shadowing. \emph{iForest-Biogeosciences and Forestry},
\bold{7(4)}, 227. \doi{10.3832/ifor1084-007}

Strahler, A. H., Jupp, D. L. B., Woodcock, C. E., Schaaf, C. B., Yao, T., Zhao,
F., Yang, X., Lovell, J., Culvenor, D., Newnham, G., Ni-Miester, W., &
Boykin-Morris, W. (2008). Retrieval of forest structural parameters using a
ground-based lidar instrument (Echidna®). \emph{Canadian Journal of Remote
Sensing}, \bold{34(sup2)}, S426-S440. \doi{10.5589/m08-046}
}


\note{
The simulation process implemented in this function is computationally intensive.
Although the function currently uses the \code{\link{vroom}} function from the \CRANpkg{vroom} package for reading large files and contains fast implementations of
several critical calculations (C++ via \CRANpkg{Rcpp} package), long
computation times may be required when a large number of plots are considered,
number of points in TLS point clouds are very high, or the radius, k
or BAF sequences used in the simulation process are very long.

Using reduced point clouds (according to point cropping process implemented in the
\code{\link{normalize}} function), rather than original ones, may be
recommended in order to cut down on computing time. Another possibility would
be to specify large increments for radius
and BAF, and/or low maximum values for radius, number of trees and BAF in the \code{plot.parameters} argument. This
would make the function more efficient, though there may be a notable
loss of detail in the results generated.
}


\author{
Juan Alberto Molina-Valero, María José Ginzo Villamayor,
Manuel Antonio Novo Pérez, Adela Martínez-Calvo, Juan Gabriel Álvarez-González,
Fernando Montes and César Pérez-Cruzado.
}


\seealso{
\code{\link{tree.detection}}, \code{\link{tree.detection.multiple}}, \code{\link{distance.sampling}},
\code{\link{normalize}}.
}


\examples{
\donttest{

# Load information of trees detected from TLS point clouds data corresponding to
# plot 1 from Rioja data set

data("Rioja.data")
example.tls <- subset(Rioja.data$tree.list.tls, id == 1)

# Compute detection probabilities using distance sampling methods

example.ds <- distance.sampling(example.tls)

# Load information of trees measured in field plots corresponding to plot 1 from
# Rioja data set

example.field <- subset(Rioja.data$tree.list.field, id == 1)

# Establish directory where TXT file containing TLS point cloud corresponding to
# plot 1 from Rioja data set is located. For instance, current working directory

dir.data <- getwd()

# Download example of TXT file corresponding to plot 1 from Rioja data set

download.file(url = "https://www.dropbox.com/s/w4fgcyezr2olj9m/1.txt?raw=1",
              destfile = file.path(dir.data, "1.txt"), method = "wininet",
              mode = "wb")

# Establish directory where simulation results corresponding to plot 1 from
# Rioja data set will be saved. For instance, current working directory

dir.result <- getwd()

# Compute metrics and variables for simulated TLS and field plots corresponding
# to plot 1 from Rioja data set
# Without occlusion correction based on distance sampling methods

sim <- simulations(tree.list.tls = example.tls, tree.list.field = example.field,
                   plot.parameters = list(radius.max = 20, k.tree.max = 30,
                                          BAF.max = 4),
                   dir.data = dir.data, dir.result = dir.result)

# With occlusion correction based on distance sampling methods

sim <- simulations(tree.list.tls = example.tls, distance.sampling = example.ds,
                   tree.list.field = example.field,
                   plot.parameters = list(radius.max = 20, k.tree.max = 30,
                                          BAF.max = 4),
                   dir.data = dir.data,
                   dir.result = dir.result)
# Only circular fixed area plot design with non-default parameters, and with occlusion
# correction based on distance sampling methods

sim <- simulations(tree.list.tls = example.tls, distance.sampling = example.ds,
                   tree.list.field = example.field,
                   plot.parameters <- list(radius.max = 20,
                                           radius.increment = .25,
                                           num.trees = 50),
                   dir.data = dir.data, dir.result = dir.result)
# Only k-tree plot design with non-default parameters, and with occlusion
# correction based on distance sampling methods

sim <- simulations(tree.list.tls = example.tls, distance.sampling = example.ds,
                   tree.list.field = example.field,
                   plot.parameters <- list(k.tree.max = 30, num.trees = 50),
                   dir.data = dir.data, dir.result = dir.result)
# Only angle-count plot design with non-default parameters, and with occlusion
# correction based on distance sampling methods

sim <- simulations(tree.list.tls = example.tls, distance.sampling = example.ds,
                   tree.list.field = example.field,
                   plot.parameters <- list(BAF.max = 5, BAF.increment = 0.5,
                                           num.trees = 50),
                   dir.data = dir.data, dir.result = dir.result)

}
}
