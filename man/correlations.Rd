
\encoding{UTF-8}


\name{correlations}


\alias{correlations}


\title{
Correlation Between Field Estimations and TLS Metrics
}


\description{
Computes correlations between variables estimates from field data and
metrics derived from TLS data. Field estimates and TLS metrics for a common
set of plots are required in order to compute correlations. These data must
be obtained from any of the three different plot designs currently available (fixed
area, k-tree and angle-count), and correspond to plots with incremental values
for the plot design parameter (radius, k and BAF, respectively). Two
correlation measures are implemented: Pearson's correlation coefficient and
Spearman's \emph{rho}. In addition to estimating these measures, tests for
association are also executed, and interactive line charts graphically
representing correlations are generated.
}


\usage{
correlations(simulations,
             variables = c("N", "G", "V", "d", "dg", "d.0", "h", "h.0"),
             method = c("pearson", "spearman"), save.result = TRUE,
             dir.result = NULL)
}


\arguments{

  \item{simulations}{
  List including estimated variables based on field data and metrics
  derived from TLS data. The structure and format must be analogous to output
  returned by the \code{\link{simulations}} function. Specifically, it must have
  at least one of the following named elements:
    \itemize{
      \item \code{fixed.area.plot}: data frame with field estimates and TLS
      metrics under a circular fixed area plot design. Each row corresponds to
      a (plot, radius) pair, and all or any of the following
      columns are included:

      Plot identification and radius:
        \itemize{
          \item \code{id}, \code{radius}: same description and format as
          indicated for same named columns of \code{fixed.area.plot} in
          \code{\link{simulations}} \sQuote{Value}.
        }
      Variables estimated on the basis of simulated field plots:
        \itemize{
          \item \code{N}, \code{G}, \code{V}, \code{d}, \code{dg},
          \code{dgeom}, \code{dharm}, \code{h}, \code{hg}, \code{hgeom},
          \code{hharm}, \code{d.0}, \code{dg.0}, \code{dgeom.0},
          \code{dharm.0}, \code{h.0}, \code{hg.0}, \code{hgeom.0},
          \code{hharm.0}: same description and format as indicated for same
          named columns of \code{fixed.area.plot} in \code{\link{simulations}} \sQuote{Value}.
        }
      TLS metrics derived from simulated TLS plots:
        \itemize{
          \item \code{N.tls}, \code{N.hn}, \code{N.hr}, \code{N.hn.cov},
          \code{N.hr.cov}, \code{N.sh}, \code{num.points}, \cr
          \code{num.points.est}, \code{num.points.hom},
          \code{num.points.hom.est}, \code{G.tls}, \cr \code{G.hn},
          \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov}, \code{G.sh},
          \code{V.tls}, \code{V.hn}, \code{V.hr}, \code{V.hn.cov}, \cr
          \code{V.hr.cov}, \code{V.sh}, \code{d.tls}, \code{dg.tls},
          \code{dgeom.tls}, \code{dharm.tls}, \code{h.tls}, \cr \code{hg.tls},
          \code{hgeom.tls}, \code{hharm.tls}, \code{d.0.tls}, \code{dg.0.tls},
          \code{dgeom.0.tls}, \code{dharm.0.tls}, \code{h.0.tls},
          \code{hg.0.tls}, \code{hgeom.0.tls}, \code{hharm.0.tls}, \code{P01},
          \code{P05}, \code{P10}, \code{P20}, \code{P25}, \code{P30},
          \code{P40}, \code{P50}, \code{P60}, \code{P70}, \code{P75},
          \code{P80}, \code{P90}, \code{P95}, \code{P99}: same description and
          format as indicated for same named columns of \code{fixed.area.plot}
          in \code{\link{simulations}} \sQuote{Value}.
        }
      If the \code{fixed.area.plot} element is included in \code{simulations}
      argument, it must contain at least \code{id} and \code{radius} columns,
      one of the field estimates columns and one of the TLS metrics
      columns.

      \item \code{k.tree.plot}: data frame with field estimates and TLS
      metrics under the k-tree plot design. Each row corresponds to a
      (plot, k) pair, and all or any of the following columns are included:

      Plot identification and k:
        \itemize{
          \item \code{id}, \code{k}: same description and format as indicated
          for same named columns of \code{k.tree.plot} in code{\link{simulations}} \sQuote{Value}.
        }
      Variables estimates on the basis of simulated field plots:
        \itemize{
          \item \code{N}, \code{G}, \code{V}, \code{d}, \code{dg},
          \code{dgeom}, \code{dharm}, \code{h}, \code{hg}, \code{hgeom},
          \code{hharm}, \code{d.0}, \code{dg.0}, \code{dgeom.0},
          \code{dharm.0}, \code{h.0}, \code{hg.0}, \code{hgeom.0},
          \code{hharm.0}: same description and format as indicated for same
          named columns of \code{k.tree.plot} in \code{\link{simulations}} \sQuote{Value}.
        }
      TLS metrics derived from simulated TLS plots:
        \itemize{
          \item \code{N.tls}, \code{N.hn}, \code{N.hr}, \code{N.hn.cov},
          \code{N.hr.cov}, \code{N.sh}, \code{num.points}, \cr
          \code{num.points.est}, \code{num.points.hom},
          \code{num.points.hom.est}, \code{G.tls}, \cr \code{G.hn},
          \code{G.hr}, \code{G.hn.cov}, \code{G.hr.cov}, \code{G.sh},
          \code{V.tls}, \code{V.hn}, \code{V.hr}, \cr \code{V.hn.cov},
          \code{V.hr.cov}, \code{V.sh}, \code{d.tls}, \code{dg.tls},
          \code{dgeom.tls}, \code{dharm.tls}, \code{h.tls}, \code{hg.tls},
          \code{hgeom.tls}, \code{hharm.tls}, \code{d.0.tls}, \code{dg.0.tls},
          \cr \code{dgeom.0.tls}, \code{dharm.0.tls}, \code{h.0.tls},
          \code{hg.0.tls}, \code{hgeom.0.tls}, \cr \code{hharm.0.tls},
          \code{P01}, \code{P05}, \code{P10}, \code{P20}, \code{P25},
          \code{P30}, \code{P40}, \code{P50}, \code{P60}, \code{P70},
          \code{P75}, \code{P80}, \code{P90}, \code{P95}, \code{P99}: same
          description and format as indicated for same named columns of
          \code{k.tree.plot} in \code{\link{simulations}} \sQuote{Value}.
        }
      If a \code{k.tree.plot} element is included in the \code{simulations} argument,
      it must include at least \code{id} and \code{k} columns, one of the field
      estimate columns, and one of the TLS metrics columns.

      \item \code{angle.count.plot}: data frame with field estimates and TLS
      metrics under the angle-count plot design. Each row corresponds
      to a (plot, BAF) pair, and all or one of the following
      columns are included:

      Plot identification and BAF:
        \itemize{
          \item \code{id}, \code{BAF}: same description and format as indicated
          for same named columns of \code{angle.count.plot} in \code{\link{simulations}} \sQuote{Value}.
        }
      Variables estimated on the basis of simulated field plots:
        \itemize{
          \item \code{N}, \code{G}, \code{V}, \code{d}, \code{dg},
          \code{dgeom}, \code{dharm}, \code{h}, \code{hg}, \code{hgeom},
          \code{hharm}, \code{d.0}, \code{dg.0}, \code{dgeom.0},
          \code{dharm.0}, \code{h.0}, \code{hg.0}, \code{hgeom.0},
          \code{hharm.0}: same description and format as indicated for same
          named columns of \code{angle.count.plot} in \code{\link{simulations}} \sQuote{Value}.
        }
      TLS metrics derived from simulated TLS plots:
        \itemize{
          \item \code{N.tls}, \code{N.pam}, \code{num.points},
          \code{num.points.est}, \code{num.points.hom}, \cr
          \code{num.points.hom.est}, \code{G.tls}, \code{G.pam}, \code{V.tls},
          \code{V.pam}, \code{d.tls}, \code{dg.tls}, \code{dgeom.tls},
          \code{dharm.tls}, \code{h.tls}, \code{hg.tls}, \code{hgeom.tls},
          \code{hharm.tls}, \cr \code{d.0.tls}, \code{dg.0.tls},
          \code{dgeom.0.tls}, \code{dharm.0.tls}, \code{h.0.tls}, \cr
          \code{hg.0.tls}, \code{hgeom.0.tls}, \code{hharm.0.tls}, \code{P01},
          \code{P05}, \code{P10}, \code{P20}, \code{P25}, \code{P30},
          \code{P40}, \code{P50}, \code{P60}, \code{P70}, \code{P75},
          \code{P80}, \code{P90}, \code{P95}, \code{P99}: same description and
          format as indicated for same named columns of \code{angle.count.plot}
          in \code{\link{simulations}} \sQuote{Value}.
        }
      If the \code{angle.count.plot} element is included in the \code{simulations}
      argument, it must contain at least \code{id} and \code{BAF} columns, one
      of the field estimates columns and one of the TLS metrics columns.
    }
  }

  \item{variables}{
  Optional character vector naming field estimates for which correlations
  between these and all the available TLS metrics will be computed. If
  this argument is specified by the user, it must include at least one of the
  following character strings: \dQuote{\code{N}}, \dQuote{\code{G}},
  \dQuote{\code{V}}, \dQuote{\code{d}}, \dQuote{\code{dg}},
  \dQuote{\code{dgeom}}, \dQuote{\code{dharm}}, \dQuote{\code{d.0}},
  \dQuote{\code{dg.0}}, \dQuote{\code{dgeom.0}}, \dQuote{\code{dharm.0}},
  \dQuote{\code{h}}, \dQuote{\code{hg}}, \dQuote{\code{hgeom}},
  \dQuote{\code{hharm}}, \dQuote{\code{h.0}}, \dQuote{\code{hg.0}},
  \dQuote{\code{hgeom.0}}, or \dQuote{\code{hharm.0}}. If this argument is not
  specified by the user, it will be set to

  \code{c("N", "G", "V", "d", "dg", "d.0", "h", "h.0")} by default. In both
  cases all the elements in the \code{simulations} argument
  must include at least the columns corresponding to the field estimates specified in the
  \code{variables} argument.
  }

  \item{method}{
  Optional character vector naming which correlation measurements will be used.
  If this argument is specified by the user, it must include at least one of
  the following character strings: \dQuote{\code{pearson}} or
  \dQuote{\code{spearman}}. If this argument is not specified by the user, it
  will be set to \code{c("pearson", "spearman")} by default.
  }

  \item{save.result}{
  Optional logical indicating wheter or not the output files described in
  \sQuote{Output Files} section must be saved in \code{dir.result} or not. If
  this argument is not specified by the user, it will be set to \code{TRUE} by
  default and, as a consequence, the output files will be saved.
  }

  \item{dir.result}{
  Optional character string naming the absolute path of an existing directory
  where files described in the \sQuote{Output Files} section will be saved.
  \code{.Platform$file.sep} must be used as the path separator in
  \code{dir.result}. If this argument is not specified by the user, and
  \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and,
  as consequence, the current working directory of the \R process will be
  assigned to \code{dir.result} during execution.
  }
}


\details{
For each radius, k or BAF value (according to the currently available plot
designs: circular fixed area, k-tree and angle-count), this function computes
correlations between each variable estimated from field data
specified in the \code{variables} argument and all the metrics derived from TLS
data existing in the data frames included in the \code{simulations} argument.

Two correlation measures are implemented at present: Pearson's correlation
coefficient and Spearman's \emph{rho}. For each method, in addition to the
estimated measure, the p-value of a test for association is also returned.
The \code{\link{cor.test}} function from the \pkg{utils} package is used to compute
both the estimated correlations and the p-values of the associated tests; more
details about these measures and their tests for association can be found in
the corresponding documentation. There cannot be
missing data for three or more plots, and there cannot be zero standard deviation,
in order to prevent missing correlation values for each (field estimation, TLS
metric) pair and plot design parameter (radius, k or BAF).

Apart from estimated correlations and their corresponding p-values, for ecah method, the
function also returns the plot design parameter and field
estimates, the value of the optimal correlation (i.e. the maximum of the
absolute value of available correlations) and the TLS metric to which it
corresponds.
}


\value{

  \item{correlations}{
  A list including the estimated correlations for each measure specified in the
  \code{method} argument will be generated. This will include all or any of the following named
  elements:
    \itemize{
      \item \code{pearson}: if \dQuote{\code{pearson}} is not included in
      \code{method} parameter, missing; otherwise, the list will include the
      estimated Pearson's correlations for each plot design specified in the
      \code{simulations} argument. In the latter case, the list will include all or any
      of the following named elements:
        \itemize{
          \item \code{fixed.area.plot}: if the \code{simulations} argument does not
          have an element named \dQuote{\code{fixed.area.plot}}, missing;
          otherwise, the matrix will include the estimated Pearson's correlations
          for a circular fixed area plot design. Each row will correspond to a
          radius value, and the following columns will be included:
            \itemize{
              \item \code{radius}: radius (m) of the simulated plots used for
              computing the estimated correlations.
              \item Column(s) \sQuote{\code{<x>.<y>}}: numeric column(s)
              containing estimated Pearson's correlations between
              \sQuote{\code{<x>}}, a field estimate, and \sQuote{\code{<y>}},
              a TLS metric.
            }
          \item \code{k.tree.plot}: if the \code{simulations} argument does not
          include an element named \dQuote{\code{k.tree.plot}}, missing;
          otherwise, the matrix will include the estimated Pearson's correlations
          for a k-tree plot design. Each row will correspond to a k
          value, and the following columns will be included:
            \itemize{
              \item \code{k}: number of trees (trees) of the simulated plots
              used for computing the estimated correlations.
              \item Column(s) \sQuote{\code{<x>.<y>}}: same description and
              format as indicated in
              \code{correlations$pearson$fixed.area.plot} element.
            }
          \item \code{angle.count.plot}: if the \code{simulations} argument does
          not have any element named \dQuote{\code{angle.count}},  missing;
          otherwise, the matrix will include the estimated Pearson's correlations
          for the angle-count plot design. Each row will correspond to a
          BAF value, and the following columns will be included:
            \itemize{
              \item \code{BAF}: BAF (\eqn{{m}^{2}/ha}{m^2/ha}) of the simulated plots used for
              computing the estimated correlations.
              \item Column(s) \sQuote{\code{<x>.<y>}}: same description
              and format as indicated in
              \code{correlations$pearson$fixed.area.plot} element.
            }
        }
      \item \code{spearman}: if \dQuote{\code{spearman}} is not included in
      \code{method} parameter, missing; otherwise, the list will include the
      estimated Spearman's correlations for each plot design specified in
      \code{simulations} argument. In the latter case, the structure and format
      will be analogous to that indicated for the previous element but
      estimated Pearson's correlations will be replaced by Spearman's correlations.
    }
  }

  \item{correlations.pval}{
  List containing the p-value of the test for association corresponding to each
  measure specified in \code{method} argument. The structure and format will be
  the same as indicated for the previous element but estimated
  correlations will be replaced by p-values for their corresponding tests for association.
  }

  \item{opt.correlations}{
  List containing the optimal correlations, and the names of the TLS metrics to
  which they correspond, for each measure specified in \code{method} argument.
  The list will will include all or any of the following named elements:
    \itemize{
      \item \code{pearson}: if \dQuote{\code{pearson}} is not included in
      \code{method} parameter, missing; otherwise, the list will include the optimal
      Pearson's correlations, and the names of the TLS metrics to which they
      correspond, for each plot design specified in \code{simulations}
      argument. In the latter case, it will include all or any of the following
      named elements:
        \itemize{
          \item \code{fixed.area.plot}: if \code{simulations} argument does not
          have any element named \dQuote{\code{fixed.area.plot}}, missing;
          otherwise, the data frame will include the optimal Pearson's correlations,
          and the names of the TLS metrics to which they correspond, for a
          circular fixed area plot design. Each row will correspond to a radius
          value, and the following columns will be included:
            \itemize{
              \item \code{radius}: radius (m) of the simulated plots used for
              computing the estimated correlations.
              \item Columns \sQuote{\code{<x>.cor}} and
              \sQuote{\code{<x>.metric}}: the former, numeric column(s)
              including optimal Pearson's correlations between
              \sQuote{\code{<x>}}, a field estimate, and all the available
              TLS metrics; and the latter, character column(s) will include names
              of the TLS metrics to which they correspond.
            }
          \item \code{k.tree.plot}: if the \code{simulations} argument does not
          have any element named \dQuote{\code{k.tree.plot}}, missing;
          otherwise, the data frame will include the optimal Pearson's correlations
          and the names of the TLS metrics to which they correspond for the
          k-tree plot design. Each row will correspond to a k value,
          and the following columns will be included:
            \itemize{
              \item \code{k}: number of trees (trees) of the simulated plots
              used for computing the estimated correlations.
              \item Columns \sQuote{\code{<x>.cor}} and
              \sQuote{\code{<x>.metric}}: same description and format as
              indicated in \code{opt.correlations$pearson$fixed.area.plot}
              element.
            }
          \item \code{angle.count.plot}: if the \code{simulations} argument does
          not have any element named \dQuote{\code{angle.count}}, missing;
          otherwise, the data frame will include the optimal Pearson's correlations,
          and the names of the TLS metrics to which they correspond for the
          angle-count plot design. Each row will correspond to a BAF
          value, and the following columns will be included:
            \itemize{
              \item \code{BAF}: BAF (\eqn{{m}^{2}/ha}{m^2/ha}) of the simulated plots used for
              computing the estimated correlations.
              \item Columns \sQuote{\code{<x>.cor}} and
              \sQuote{\code{<x>.metric}}: same description and format as
              indicated in \code{opt.correlations$pearson$fixed.area.plot}
              element.
            }
        }
      \item \code{spearman}: if \dQuote{\code{spearman}} is not included in
      \code{method} parameter, missing; otherwise, the list will include the
      optimal Spearman's correlations, and the names of the TLS metrics to
      which they correspond, for each plot design specified in
      \code{simulations} argument. In the latter case, the structure and
      format will be analogous to that indicated for the previous element, but
      optimal Pearson's correlations will be replaced by Spearman's correlations.
    }
  }
}


\section{Output Files}{
During the execution, if the \code{save.result} argument is \code{TRUE}, the
function prints to files the matrices and data frames included in
\code{correlations} and \code{opt.correlations} elements described in
\sQuote{Value}. Both are written without row names in
\code{dir.result} directory by using the \code{\link{write.csv}} function in the
\pkg{utils} package. The patterns used for naming these files are
\file{correlations.<plot design>.<method>.csv} and
\file{opt.correlations.<plot design>.plot.<method>.csv} for correlation
matrices and optimal correlation data frames, respectively, where
\file{<plot design>} is equal to \dQuote{\code{fixed.area.plot}},
\dQuote{\code{k.tree.plot}} or \dQuote{\code{angle.count.plot}} according to
plot design, and \file{<method>} equals \dQuote{\code{pearson}} or
\dQuote{\code{spearman}} according to the correlation measure.

Furthermore, if the \code{save.result} argument is \code{TRUE}, interactive line
charts graphically representing correlations will also be created and saved
in the \code{dir.result} directory by means of \code{\link{saveWidget}} function
in the \CRANpkg{htmlwidgets} package. Generated widgets enable users to
consult correlation data directly on the plots, select/deselect different
sets of traces, to zoom and scroll, etc. The pattern used for naming
these files is \file{correlations.<x>.<plot design>.<method>.html}, where both
\file{<plot design>} and \file{<method>} are as indicated for the previous
described files, and \file{<x>} is equal to any of elements specified in the
\code{variables} argument.
}


\note{
This function is particularly useful for further steps related to model-based
and model-assisted approaches, as correlations measure the strength of a relationship between two variables (linear for
Pearson's correlation, monotonic for Spearman's correlation).
}


\author{
Juan Alberto Molina-Valero and
Adela Martínez-Calvo.}


\seealso{
\code{\link{simulations}}, \code{\link{optimize.plot.design}}.

\code{\link{cor.test}} in \pkg{utils} package.
}


\examples{
\donttest{

# Load field estimates and TLS metrics corresponding to Rioja data set

data("Rioja.simulations")


# Establish directory where correlation results corresponding to the Rioja example
# will be saved. For instance, current working directory

# dir.result <- getwd()


# Compute correlations between field estimates and TLS metrics corresponding
# to Rioja example
# Pearson's and Spearman's correlations for variables by default

# corr <- correlations(simulations = Rioja.simulations, dir.result = dir.result)


# Pearson's and Spearman's correlations for variable 'N'

# corr <- correlations(simulations = Rioja.simulations, variables = "N",
#                      dir.result = dir.result)


# Only Pearson's correlations for variables by default

# corr <- correlations(simulations = Rioja.simulations, method = "pearson",
#                      dir.result = dir.result)


# Pearson's and Spearman's correlations corresponding to angle-count design for
# all available variables

# corr <- correlations(simulations = Rioja.simulations["angle.count"],
#                      variables <- c("N", "G", "V", "d", "dg", "dgeom", "dharm",
#                                     "d.0", "dg.0", "dgeom.0", "dharm.0", "h",
#                                     "hg", "hgeom", "hharm", "h.0", "hg.0",
#                                     "hgeom.0", "hharm.0"),
#                      dir.result = dir.result)

  }

}

