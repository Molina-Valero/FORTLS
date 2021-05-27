
\encoding{UTF-8}


\name{relative.bias}


\alias{relative.bias}


\title{
Relative Bias Between Field Estimations and TLS metrics
}


\description{
Computes relative bias between variables estimated from field data and
their TLS counterparts derived from TLS data. Field estimates and TLS metrics
for a common set of plots are required in order to compute relative bias. These
data must come from any of the three different plot designs currently available
(circular fixed area, k-tree and angle-count) and correspond to plots with incremental
values for the plot design parameter (radius, k and BAF, respectively). In
addition to computing relative bias, interactive line charts graphically
representing the values obtained between each field estimate and its related
TLS metrics are also generated.
}


\usage{
relative.bias(simulations,
              variables = c("N", "G", "V", "d", "dg", "d.0", "h", "h.0"),
              save.result = TRUE, dir.result = NULL)
}


\arguments{

  \item{simulations}{
  List containing variables estimated from field data and metrics
  derived from TLS data. The structure and format must be analogous to output
  returned by the \code{\link{simulations}} function. In particular, teh list must include
  at least one of the following named elements:
    \itemize{
      \item \code{fixed.area.plot}: same description and format as indicated
      for same named element in \code{simulations} argument of
      \code{\link{correlations}} function. The only difference is the
      columns required when it is included in the argument: in adittion to \code{id},
      \code{radius} and at least one of the field estimate columns, it must
      include at least one TLS counterpart for each field estimate condidered.

      \item \code{k.tree.plot}: same description and format as indicated for
      same named element in \code{simulations} argument of
      \code{\link{correlations}} function. The only difference is the columns required when it is included in the argument: in adittion to \code{id}, \code{k}
      and at least one of the field estimate columns, it must include at
      least one TLS counterpart for each field estimate considered.

      \item \code{angle.count.plot}: same description and format as indicated
      for the same named element in the \code{simulations} argument of
      \code{\link{correlations}} function. The only difference is the columns that are required when this element is included in the argument: in adittion to \code{id},
      \code{BAF} and at least one of the field estimate columns, it must
      contain at least one TLS counterpart for each field estimate considered.
    }
  }

  \item{variables}{
  Optional character vector naming field estimates for which the relative bias
  between them and all their available TLS counterparts will be computed. If
  this argument is specified by the user, it must contain at least one of the
  following character strings: \dQuote{\code{N}}, \dQuote{\code{G}},
  \dQuote{\code{V}}, \dQuote{\code{d}}, \dQuote{\code{dg}},
  \dQuote{\code{dgeom}}, \dQuote{\code{dharm}}, \dQuote{\code{d.0}},
  \dQuote{\code{dg.0}}, \dQuote{\code{dgeom.0}}, \dQuote{\code{dharm.0}},
  \dQuote{\code{h}}, \dQuote{\code{hg}}, \dQuote{\code{hgeom}},
  \dQuote{\code{hharm}}, \dQuote{\code{h.0}}, \dQuote{\code{hg.0}}, \cr
  \dQuote{\code{hgeom.0}}, or \dQuote{\code{hharm.0}}. If this argument is not
  specified by the user, it will be set to
  \code{c("N", "G", "V", "d", "dg", "d.0", "h", "h.0")} by default. In both
  cases, all the elements in \code{simulations} argument
  must include at least the columns corresponding to the field estimations specified in the
  \code{variables} argument.
  }

  \item{save.result}{
  Optional logical which indicates whether or not the output files described in
  \sQuote{Output Files} section must be saved in \code{dir.result}. If
  this argument is not specified by the user, it will be set to \code{TRUE} by
  default and, as a consequence, the output files will be saved.
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
For each radius, k or BAF value (according to the currently available plot
designs: circular fixed area, k-tree and angle-count), this function computes the relative
bias between each variable estimated from field data, and specified in the
\code{variables} argument, and their counterparts derived from TLS data, and
existing in the data frames included in the \code{simulations} argument. TLS
metrics considered \emph{counterparts} for each field estimate are
detailed below (see \code{\link{simulations}} \sQuote{Value} function
for details about used notation):
  \itemize{
    \item TLS counterparts for \code{N} are \code{N.tls}, \code{N.hn},
    \code{N.hr}, \code{N.hn.cov}, \code{N.hr.cov} and \code{N.sh} in the fixed
    area and k-tree plot designs; and \code{N.tls} and \code{N.pam} in the
    angle-count plot design.
    \item TLS counterparts for \code{G} are \code{G.tls}, \code{G.hn},
    \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov} and \code{G.sh} in the fixed
    area and k-tree plot designs; and \code{G.tls} and \code{G.pam} in the
    angle-count plot design.
    \item TLS counterparts for \code{V} are \code{V.tls}, \code{V.hn},
    \code{V.hr}, \code{V.hn.cov}, \code{V.hr.cov} and \code{V.sh} in the fixed
    area and k-tree plot designs; and \code{V.tls} and \code{V.pam} in the
    angle-count plot design.
    \item TLS counterparts for \code{d}, \code{dg}, \code{dgeom}, \code{dharm},
    \code{d.0}, \code{dg.0}, \code{dgeom.0}, and \code{dharm.0} are,
    respectively: \code{d.tls}, \code{dg.tls}, \code{dgeom.tls},
    \code{dharm.tls}, \code{d.0.tls}, \code{dg.0.tls}, \code{dgeom.0.tls}, \cr
    and \code{dharm.0.tls} in any of the three available plot designs.
    \item TLS counterparts for \code{h}, \code{hg}, \code{hgeom}, \code{hharm},
    \code{h.0}, \code{hg.0}, \code{hgeom.0}, and \code{hharm.0} are,
    respectively \code{h.tls}, \code{hg.tls}, \code{hgeom.tls},
    \code{hharm.tls}, \code{h.0.tls}, \code{hg.0.tls}, \code{hgeom.0.tls}, and
    \cr \code{hharm.0.tls} in any of the three available plot designs. In adittion, \code{P99} is also taken into account as a counterpart for all
    these field estimates.
  }

The relative bias between a field estimation and any of its TLS counterparts is
estimated as follows
\deqn{\frac{\frac{1}{n}\sum_{i = 1}^{n}{y_{i}} -
\frac{1}{n}\sum_{i = 1}^{n}{x_{i}}}{\frac{1}{n}\sum_{i = 1}^{n}{y_{i}}} * 100}{
((\sum(y_i)/n - \sum(x_i)/n) / (\sum(y_i)/n) * 100,}
where \eqn{y_{i}}{y_i} and \eqn{x_{i}}{x_i}  are the values of the field
estimate and its TLS counterpart, respectively, corresponding to plot
\eqn{i}{i} for \eqn{i = 1, \ldots, n}.
}


\value{

  \item{fixed.area.plot}{
  If no \dQuote{\code{fixed.area.plot}} element exists in \code{simulations}
  argument, missing; otherwise, the matrix will include estimates of the the relative bias
  under circular fixed area plot design between each field estimation specified in
  \code{variables} argument and its TLS counterpart(s) existing in the
  \dQuote{\code{fixed.area.plot}} element in the \code{simulations} argument. Each
  row will correspond to a radius value, and the following columns will be included:
    \itemize{
      \item \code{radius}: radius (m) of the simulated plots used for computing
      the estimated relative bias.
      \item Column(s) \sQuote{\code{<x>.<y>}}: numeric column(s) containing
      estimated relative bias between \sQuote{\code{<x>}}, a field estimation,
      and \sQuote{\code{<y>}}, a TLS counterpart.
    }
  }

  \item{k.tree.plot}{
  If no \dQuote{\code{k.tree.plot}} element exists in the \code{simulations}
  argument, missing; otherwise, the matrix will include the relative bias estimated
  in the k-tree plot design between each field estimation specified in
  \code{variables} argument and its TLS counterpart(s) existing in
  \dQuote{\code{k.tree.plot}} element in \code{simulations} argument. Each row will correspond to a k value, and the following columns will be included:
    \itemize{
      \item \code{k}: number of trees (trees) of the simulated plots used for
      computing the estimated relative bias.
      \item Column(s) \sQuote{\code{<x>.<y>}}: numeric column(s) containing
      estimated relative bias between \sQuote{\code{<x>}}, a field estimation,
      and \sQuote{\code{<y>}}, a TLS counterpart.
    }
  }

  \item{angle.count.plot}{
  If no \dQuote{\code{angle.count}} element exists in \code{simulations}
  argument, missing; otherwise, the matrix will contain estimated relative bias under
  angle-count plot design between each field estimation specified in the
  \code{variables} argument and its TLS counterpart(s) existing in the
  \dQuote{\code{angle.count.plot}} element in the \code{simulations} argument. Each
  row will correspond to a BAF value, and the following
  columns will be included:
    \itemize{
      \item \code{BAF}: BAF (\eqn{{m}^{2}/ha}{m^2/ha}) of the simulated plots used for computing
      the estimated relative bias.
      \item Column(s) \sQuote{\code{<x>.<y>}}: numeric column(s) containing
      estimated relative bias between \sQuote{\code{<x>}}, a field estimation,
      and \sQuote{\code{<y>}}, a TLS counterpart.
    }
  }
}


\section{Output Files}{
During the execution, if the \code{save.result} argument is \code{TRUE}, the
function will print the matrices described in the \sQuote{Value} section to files. These
are written without row names in \code{dir.result} directory using
\code{\link{write.csv}} function from the \pkg{utils} package. The pattern used
for naming these files is \file{RB.<plot design>.csv}, where
\file{<plot design>} is equal to \dQuote{\code{fixed.area.plot}},
\dQuote{\code{k.tree.plot}} or \dQuote{\code{angle.count.plot}} is according to the
plot design.

Furthermore, if the \code{save.result} argument is \code{TRUE}, interactive line
charts graphically representing relative bias values will also be created and saved
in the \code{dir.result} directory by means of the \code{\link{saveWidget}} function
in the \CRANpkg{htmlwidgets} package. Generated widgets allow users to
consult relative bias data directly on the plots, select/deselect different
sets of traces, to zoom and scroll, and so on. The pattern used for naming
these files is \file{RB.<x>.<plot design>.html}, where \file{<plot design>} is
indicated for the previously described files, and \file{<x>} equals \code{N},
\code{G}, \code{V}, \code{d} and/or \code{h} according to the \code{variables}
argument. All relative biases related to diameters are
plotted in the same chart (files named as \file{RB.d.<plot design>.html}), and
the same applies to those related to heights (files named as
\file{RB.h.<plot design>.html}).
}



\note{
The results obtained using this function are
merely descriptive, and they do not guarantee any type of statistical accuracy
in using TLS metrics instead of field estimations in order to estimate forest
attributes of interest.
}


\author{
Juan Alberto Molina-Valero, María José Ginzo Villamayor,
Manuel Antonio Novo Pérez, Adela Martínez-Calvo, Juan Gabriel Álvarez-González,
Fernando Montes and César Pérez-Cruzado.
}


\seealso{
\code{\link{simulations}}, \code{\link{correlations}}.
}


\examples{
\donttest{

# Load variables estimated from field data, and TLS metrics
# corresponding to Rioja data set

data("Rioja.simulations")


# Establish directory where relative bias results corresponding to Rioja example
# will be saved. For instance, current working directory

dir.result <- getwd()


# Compute relative bias between field-based estimates of TLS metrics
# corresponding to Rioja example
# Relative bias for variables by default

rb <- relative.bias(simulations = Rioja.simulations, dir.result = dir.result)


# Relative bias for variable 'N'

rb <- relative.bias(simulations = Rioja.simulations, variables = "N",
                    dir.result = dir.result)


# Relative bias corresponding to angle-count design for all available variables

rb <- relative.bias(simulations = Rioja.simulations["angle.count.plot"],
                    variables <- c("N", "G", "V", "d", "dg", "dgeom", "dharm",
                                   "d.0", "dg.0", "dgeom.0", "dharm.0", "h",
                                   "hg", "hgeom", "hharm", "h.0", "hg.0",
                                   "hgeom.0", "hharm.0"),
                    dir.result = dir.result)
}
}
