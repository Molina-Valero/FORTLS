
\encoding{UTF-8}


\name{optimize.plot.design}


\alias{optimize.plot.design}


\title{
Optimize Plot Design Based on Optimal Correlations
}


\description{
Generation of interactive heatmaps graphically represent the optimal
correlations between variables estimated from field data, and metrics
derived from TLS data. These data must be derived from any of the three different
plot designs currently available (circular fixed area, k-tree and angle-count) and
correspond to plots with incremental values for the plot design parameter
(radius, k and BAF, respectively). In addition, correlation measures that are currently
admissible are Pearson's correlation coefficient and/or Spearman's \emph{rho}.
}


\usage{
optimize.plot.design(correlations,
                     variables = c("N", "G", "V", "d", "dg", "d.0", "h", "h.0"),
                     dir.result = NULL)
}


\arguments{

  \item{correlations}{
  List including the optimal correlations between field estimations and TLS
  metrics. The structure and format must be analogous to the
  \code{opt.correlations} element in the output returned by the
  \code{\link{correlations}} function. In particular, the list must include at least one
  of the following named elements:
    \itemize{
      \item \code{pearson}: list containing the optimal Pearson's correlations,
      and the names of the TLS metrics to which they correspond. It must include
      at least one of the following named elements:
        \itemize{
          \item \code{fixed.area.plot}: data frame containing the optimal
          Pearson's correlations, and the names of the TLS metrics to which
          they correspond, under circular fixed area plot design. Each row will
          correspond to a radius value, and the data frame will include the following columns:
            \itemize{
              \item \code{radius}: radius (m) of the simulated plots used for
              computing the estimated correlations.
              \item Columns \sQuote{\code{<x>.cor}} and
              \sQuote{\code{<x>.metric}}: the former, numeric column(s)
              containing optimal Pearson's correlations between
              \sQuote{\code{<x>}}, a field estimate, and the available TLS
              metrics when \code{\link{correlations}} function was executed;
              and the latter, character column(s) containing names of the TLS
              metrics to which they correspond.
            }
          If \code{fixed.area.plot} is included in \code{pearson} element, it
          must contain at least the \code{radius} column and a
          (\sQuote{\code{<x>.cor}}, \sQuote{\code{<x>.metric}}) pair of
          columns corresponding to the same field estimation.
          \item \code{k.tree.plot}: data frame including the optimal Pearson's
          correlations and the names of the TLS metrics to which they
          correspond, under the k-tree plot design. Each row will
          correspond to a k value, and the following columns will be included:
            \itemize{
              \item \code{k}: number of trees (trees) of the simulated plots
              used for computing the estimated correlations.
              \item Columns \sQuote{\code{<x>.cor}} and
              \sQuote{\code{<x>.metric}}: same description and format as
              indicated in \code{correlations$pearson$fixed.area.plot} element.
            }
          If \code{k.tree.plot} is included in \code{pearson} element, it must
          contain at least \code{k} column, and a (\sQuote{\code{<x>.cor}},
          \sQuote{\code{<x>.metric}}) pair of columns corresponding to the same
          field estimation.
          \item \code{angle.count.plot}: data frame including the optimal
          Pearson's correlations and the names of the TLS metrics to which
          they correspond, for the angle-count plot design. Each row will
          correspond to a BAF value, and the data frame will include the following columns:
            \itemize{
              \item \code{BAF}: BAF (\eqn{{m}^{2}/ha}{m^2/ha}) of the simulated plots used to
              compute the estimated correlations.
              \item Columns \sQuote{\code{<x>.cor}} and
              \sQuote{\code{<x>.metric}}: same description and format as
              indicated in \code{correlations$pearson$fixed.area.plot} element.
            }
          If the \code{angle.count.plot} is included in the \code{pearson} element, it
          must include at least the \code{BAF} column and a
          (\sQuote{\code{<x>.cor}}, \sQuote{\code{<x>.metric}}) pair of columns
          corresponding to the same field estimation.
        }

      \item \code{spearman}: list containing the optimal Spearman's
      correlations, and the names of the TLS metrics to which they correspond.
      The structure and format will be analogous to that indicated for the
      previous element but with optimal Pearson's correlations replaced by
      Spearman's correlations.
    }
  }

  \item{variables}{
  Optional character vector naming field estimations whose optimal correlations
  will be represented graphically in the heatmaps generated during the
  execution. If this argument is specified by the user, it must include at
  least one of the following character strings: \dQuote{\code{N}},
  \dQuote{\code{G}}, \dQuote{\code{V}}, \dQuote{\code{d}}, \dQuote{\code{dg}},
  \dQuote{\code{dgeom}}, \dQuote{\code{dharm}}, \dQuote{\code{d.0}},
  \dQuote{\code{dg.0}}, \dQuote{\code{dgeom.0}}, \dQuote{\code{dharm.0}},
  \dQuote{\code{h}}, \dQuote{\code{hg}}, \dQuote{\code{hgeom}},
  \dQuote{\code{hharm}}, \dQuote{\code{h.0}}, \dQuote{\code{hg.0}},
  \dQuote{\code{hgeom.0}}, or \dQuote{\code{hharm.0}}. If this argument is not
  specified by the user, it will be set to
  \code{c("N", "G", "V", "d", "dg", "d.0", "h", "h.0")} by default. In both
  cases, all data frames in the \code{correlations}
  argument must have at least the (\sQuote{\code{<x>.cor}},
  \sQuote{\code{<x>.metric}}) pairs corresponding to the field estimations
  specified in the \code{variables} argument.
  }

  \item{dir.result}{
  Optional character string naming the absolute path of an existing directory
  where files described in \sQuote{Output Files} section will be saved.
  \code{.Platform$file.sep} must be used as the path separator in
  \code{dir.result}. If this argument is not specified by the user, it will be set
  to \code{NULL} by default and, as consequence, the current working directory
  of the \R process will be assigned to \code{dir.result} during the execution.
  }
}


\details{
This function represents graphically, by means of interactive heatmaps, the strongest correlations (positive or negative) for each plot design and size simulated, between the
estimated variables based on field data specified in the \code{variables}
argument, and metrics derived from TLS data, under circular fixed area, k-tree and/or
angle-count plot designs.

Two correlation measures are implemented at present: Pearson’s correlation
coefficient and Spearman’s rho. Hence, only optimal correlations based on \code{correlations} arguments will be taken into account during the
execution.

For each correlation measure and plot design, at least one no missing value for optimal correlations must be represented;
otherwise, execution will be stopped, and an error message will appear. In
addition, at least two different no missing values for
optimal correlations are required to ensure that the colour palette is
correctly applied when the heatmap is generated.
}


\value{
Invisible \code{NULL}.
}


\section{Output Files}{
During the execution, interactive heatmaps graphically representing optimal
correlations values between field estimations and TLS metrics are created and
saved in \code{dir.result} directory by means of the \code{\link{saveWidget}}
function in the \CRANpkg{htmlwidgets} package. The widgets generated allow users
to consult optimal correlations values and TLS metrics to which they correspond
directly on the plots, to zoom and scroll, and so on. The pattern used for
naming these files is \file{opt.correlations.<plot design>.<method>.html},
where \file{<plot design>} equals \dQuote{\code{fixed.area.plot}},
\dQuote{\code{k.tree.plot}} or \dQuote{\code{angle.count.plot}} according to
plot design, and \file{<method>} equals \dQuote{\code{pearson}} or
\dQuote{\code{spearman}} according to correlation measure.
}


\note{
This function is key to choosing the best possible plot design (in terms of
correlation measures) considering all variables of interest before establishing
definitive sampling design.
}


\author{
Juan Alberto Molina-Valero, María José Ginzo Villamayor,
Manuel Antonio Novo Pérez, Adela Martínez-Calvo, Juan Gabriel Álvarez-González,
Fernando Montes and César Pérez-Cruzado.
}


\seealso{
\code{\link{correlations}}.
}


\examples{
\donttest{

# Load field estimations and TLS metrics corresponding to Rioja data set

data("Rioja.simulations")


# Compute correlations between field estimations and TLS metrics corresponding
# to Rioja example, and select optimal correlations results

corr <- correlations(simulations = simulations,
                     variables = c("N", "G", "V", "d", "dg", "dgeom","dharm",
                                   "d.0", "dg.0", "dgeom.0", "dharm.0", "h",
                                   "hg", "hgeom", "hharm", "h.0", "hg.0",
                                   "hgeom.0", "hharm.0"),
                     save.result = FALSE)

opt.corr <- corr$opt.correlations


# Establish directory where optimal correlations heatmaps corresponding to Rioja
# example will be saved. For instance, current working directory

dir.result <- getwd()


# Generate heatmaps for optimal correlations between field estimations and TLS
# metrics corresponding to Rioja example
# Optimal Pearson's and Spearman's correlations for variables by default

optimize.plot.design(correlations = opt.corr, dir.result = dir.result)


# Optimal Pearson's and Spearman's correlations for variables 'N', 'G' and 'V'

optimize.plot.design(correlations = opt.corr, variables = c("N", "G", "V"),
                     dir.result = dir.result)


# Only optimal Pearson's correlations for variables by default

optimize.plot.design(correlations = opt.corr["pearson"],
                     dir.result = dir.result)


# Optimal Pearson's and Spearman's correlations corresponding to angle-count
# design for all available variables

optimize.plot.design(
  correlations = list(pearson = opt.corr$pearson["angle.count.plot"],
                      spearman = opt.corr$spearman["angle.count.plot"]),
  variables <- c("N", "G", "V", "d", "dg", "dgeom", "dharm", "d.0", "dg.0",
                 "dgeom.0", "dharm.0", "h", "hg", "hgeom", "hharm", "h.0",
                 "hg.0", "hgeom.0", "hharm.0"),
  dir.result = dir.result)
}
}
