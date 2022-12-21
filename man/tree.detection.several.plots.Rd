\encoding{UTF-8}

\name{tree.detection.several.plots}

\alias{tree.detection.several.plots}

\title{Tree-Level Variables Estimation for Several Plots}

\description{This function integrates both, the \code{\link{normalize}} and \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multi.scan}} functions, generating the same \sQuote{Output Files} as indicated for these, and it returs the same \sQuote{Value} as described for \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multi.scan}} respectively. However, this function is designed for working with several plots, producing a list of all scans considered automatically from LAS files.}

\usage{
tree.detection.several.plots(las.list, id = NULL, file = NULL,

                             scan.approach = "single",
                             pcd.red = NULL, normalized = NULL,
                             x.center = NULL, y.center = NULL,
                             x.side = NULL, y.side = NULL,
                             max.dist = NULL, min.height = NULL, max.height = 50,
                             algorithm.dtm = "knnidw", res.dtm = 0.2,
                             csf = list(cloth_resolution = 0.5),
                             RGB = NULL,

                             single.tree = NULL,
                             dbh.min = 4, dbh.max = 200, h.min = 1.3,
                             ncr.threshold = 0.1,
                             tls.resolution = NULL, tls.precision = NULL,
                             stem.section = NULL, breaks = NULL,
                             slice = 0.1, understory = NULL, bark.roughness = 2,
                             den.type = 1, d.top = NULL,

                             plot.attributes = NULL,

                             dir.data = NULL, save.result = TRUE, dir.result = NULL)}

\arguments{

  \item{las.list}{Character vector containing the names of all LAS files for analysis and belonging to TLS point cloud, including .las extension (see \sQuote{Examples})}

  \item{id}{Optional vector with plots identification encoded as character string or numeric. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the plots will be encoded with correlative numbers from 1 to n plots.}

  \item{file}{Optional vector containing files name identification encoded as character string or numeric value. If it is null, file will be encoded as \code{id} by default.}

  \item{scan.approach}{Character parameter indicating TLS single-scan (\sQuote{single}) or TLS multi-scan approach or SLAM point clouds (\sQuote{multi}) approaches. If this argument is not specified by the user, it will be set to \sQuote{multi} approach.}

  \item{pcd.red}{Optional argument to indicate if point cloud density must be reduced to detect trees.}

  \item{normalized}{Optional argument to establish as \code{TRUE} when point cloud is already normalized.}

  \item{x.center}{Planimetric x center coordinate of point cloud data.}

  \item{y.center}{Planimetric y center coordinate of point cloud data.}

  \item{x.side}{x-side (m) of the plot when the plot is square or rectangular.}

  \item{y.side}{y-side (m) of the plot when the plot is square or rectangular.}

  \item{max.dist}{Optional maximum horizontal distance (m) considered from the plot centre. All points farther than \code{max.dist} will be discarded after the normalization process. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, all points will be used in processing, with \code{max.dist} representing the farthest point.}

  \item{min.height}{Optional minimum height (m) considered from ground level. All points below \code{min.height} will be discarded after the normalization process. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, all points will be used in processing, with \code{min.height} representing the lowest point.}

  \item{max.height}{Optional maximum height (m) considered from ground level. All points above \code{max.height} will be discarded after the normalization process. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, all points will be used in processing, with \code{max.height} representing the highest point.}

    \item{algorithm.dtm}{Algorithm used to generate the digital terrain model (DTM) from the TLS point cloud. There are two posible options based on spatial interpolation: \sQuote{tin} and \sQuote{knnidw} (see \sQuote{Details}). If this argument is not specified by the user, it will be set to \sQuote{knnidw} algorithm.}

    \item{res.dtm}{Numeric parameter. Resolution of the DTM generated to normalize point cloud (see \sQuote{Details}). If this argument is not specified by the user, it will be set to 0.2 m.}

    \item{csf}{List containing parameters of CSF algorithm:}

\itemize{

  \item {\code{cloth_resolution}}: by default 0.5.}

  \item{RGB}{Logical parameter useful when point clouds are colorized, thus including values of RGB colors. It is based on the Green Leaf Algorithm (GLA) (see \sQuote{Details}).}

  \item{single.tree}{Optional argument to indicate if there is only one tree.}

  \item{dbh.min}{Optional minimum \emph{dbh} (cm) considered for detecting trees. By default it will be set at 4 cm.}

  \item{dbh.max}{Optional maximum \emph{dbh} (cm) considered for detecting trees. By default it will be set at 200 cm.}

  \item{h.min}{Optional minimum \emph{h} (m) considered for detecting trees. By default it will be set at 1.3 m.}

  \item{ncr.threshold}{Local surface variation (also known as normal change rate, NCR). By default it will be set as 0.1. For better understanding of this argument see \sQuote{Details}.}

  \item{tls.resolution}{List containing parameters of TLS resolution. This can be defined by the angle aperture:}

\itemize{

  \item {\code{horizontal.angle}}: horizontal angle aperture (degrees).

  \item {\code{vertical.angle}}: vertical angle aperture (degrees).

  \item {\code{point.dist}}: distance (mm) between two consecutive points.

  \item {\code{tls.dist}}: distance (m) from TLS at which two consecutive points are separated by \code{point.dist}.

  If this argument is not specified by the user, it will be set to NULL by default and, as a consequence the function will stop giving an error message.}

  \item{tls.precision}{Optional argument indicating the average point cloud precision in cm.}

  \item{stem.section}{Section free of noise (shurb, branches, etc.) considered to detect trees. If not specified, an automatic internal algorithm will be applied (see \sQuote{Details}).}

  \item{breaks}{Height above ground level (m) of slices considered for detecting trees. By default it will be considered all possible sections from 0.1 m to maximum height by 0.3 m intervals (+/- 5 cm).}

  \item{slice}{Slice width considered for detecting trees. By default it will be considered as 0.1 m.}

  \item{understory}{Optional argument to indicate if there is dense understory vegetation.}

  \item{bark.roughness}{Bark roughness established in 3 degrees (1 < 2 < 3). By default it will be considered as 2.}

  \item{den.type}{Numeric argument indicating the dendrometic type used to estimate volumen when there are not sections enough to fit a taper equation. Dendrometrics types available are the following: cylinder = 0, paraboloid = 1 (by default), cone = 2 and neiloid = 3.}

  \item{d.top}{Top stem diameter (cm) considered to estimate commercial timber volume.}

  \item{plot.attributes}{Data frame with attributes at plot level. It must contain a column named \code{id} (character string or numeric value) with encoding coinciding with that used in \code{id} argument of \code{\link{normalize}} for identifying plots. If there are strata, another column named \sQuote{stratum} (numeric) will be required for other functionalities of \pkg{FORTLS} (see, for instance, \code{\link{estimation.plot.size}} or \code{\link{metrics.variables}}). If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the function will not add these possible plot attributes.}

  \item{dir.data}{Optional character string naming the absolute path of the directory where LAS files containing TLS point clouds are located. \code{.Platform$file.sep} must be used as the path separator in \code{dir.dat}. If this argument is not specified by the user, it will be set to \code{NULL} by default and, as consequence, the current working directory of the \R process will be assigned to \code{dir.dat} during the execution.}

  \item{save.result}{Optional logical which indicates whether or not the output files described in \sQuote{Output Files} section should be saved in the \code{dir.result}. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as a consequence, the output files will be saved.}

  \item{dir.result}{Optional character string naming the absolute path of an existing directory where files described in \sQuote{Output Files} section will be saved. \code{.Platform$file.sep} must be used as the path separator in \code{dir.result}. If this argument is not specified by the user, and \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.result} during the execution.}}


\details{

See \code{\link{normalize}}, \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}} for further details.}

\value{Data frame with the same description and format as \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}} \sQuote{Values}. In this case, the \code{id} of plots will be encoded with correlative numbers from 1 to n, where n is the number of LAS files included in \code{files} argument, and \code{file} column will be encoded as \code{id}, but including .las extension.}


\section{Output Files}{

At the end of the tree detection process, if the \code{save.result} argument is \code{TRUE}, the function will save both, the reduced point clouds as TXT files encoded according to \code{file} column of \sQuote{Value}; and the data frame with the tree list described in \sQuote{Value} as CSV file (see \code{\link{normalize}} and \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multi.scan}} \sQuote{Output files}). All outputs are written without row names in the \code{dir.result} directory using \code{\link{vroom_write}} function from \CRANpkg{vroom} package.}


\note{

This function has been developed for working with several plots, which will be the most common situation in forest inventory approaches. Nevertheless, several LAS files are not provided as examples due to problems with memory capacity.}


\author{
Juan Alberto Molina-Valero
and Adela Martínez-Calvo.}

\seealso{
\code{\link{normalize}}, \code{\link{tree.detection.single.scan}},\code{\link{tree.detection.multi.scan}}, \code{\link{distance.sampling}}, \code{\link{estimation.plot.size}}, \code{\link{simulations}}, \code{\link{metrics.variables}}.}

\examples{

\donttest{

# Establishment of working directories (optional)
# By default here we propose the current working directory of the R process

dir.data <- getwd()
dir.result <- getwd()

# Loading example data (LAZ files) to dir.data

download.file("https://www.dropbox.com/s/2nem6tfbsw7eyd7/Galicia_multi_scan.laz?dl=1",
              destfile = file.path(dir.data, "Galicia_multi_scan.laz"),
              mode = "wb")

download.file("https://www.dropbox.com/s/7pnzxsjf8d2x84f/Go_Fagus_SLAM.laz?dl=1",
              destfile = file.path(dir.data, "Go_Fagus_SLAM.laz"),
              mode = "wb")


# Obtaining a vector with names of LAZ files located in dir.data

files <- list.files(pattern = "laz$", path = dir.data)


# Tree detection

tree.tls <- tree.detection.several.plots(las.list =
                                         c("Galicia_multi_scan.laz",
                                         "Go_Fagus_SLAM.laz"),

                                         scan.approach = "multi",

                                         breaks = 1.3, stem.section = c(1,5))

  }
}

