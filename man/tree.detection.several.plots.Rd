\encoding{UTF-8}

\name{tree.detection.several.plots}

\alias{tree.detection.several.plots}

\title{Tree-Level Variables Estimation for Several Plots}

\description{This function integrates both, the \code{\link{normalize}} and \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multi.scans}} functions, generating the same \sQuote{Output Files} as indicated for these, and it returs the same \sQuote{Value} as described for \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multi.scans}} respectively. However, this function is designed for working with several plots, producing a list of all scans considered automatically from LAS files.}

\usage{
tree.detection.several.plots(las.list, id = NULL, file = NULL,


                             scan.approach = "single",


                             x.center = NULL, y.center = NULL,

                             max.dist = NULL, min.height = NULL, max.height = NULL,

                             algorithm.dtm = "tin", res.dtm = 0.2, csf = list(cloth_resolution = 0.5),



                             dbh.min = 4, dbh.max = 200, h.min = 1.3,

                             ncr.threshold = 0.1,

                             tls.resolution = NULL, tls.precision = NULL,

                             breaks = NULL, stem.section = NULL,

                             d.top = NULL,


                             plot.attributes = NULL,


                             dir.data = NULL, save.result = TRUE, dir.result = NULL)}

\arguments{

  \item{las.list}{Character vector containing the names of all LAS files for analysis and belonging to TLS point cloud, including .las extension (see \sQuote{Examples})}

  \item{id}{Optional vector with plots identification encoded as character string or numeric. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the plots will be encoded with correlative numbers from 1 to n plots.}

  \item{file}{Optional vector containing files name identification encoded as character string or numeric value. If it is null, file will be encoded as \code{id} by default.}

  \item{scan.approach}{Character parameter indicating TLS single-scan (\sQuote{single}) or TLS multi-scan approach or SLAM point clouds (\sQuote{multi}) approaches. If this argument is not specified by the user, it will be set to \sQuote{multi} approach.}

  \item{Normalize:}

  \itemize{

    \item \code{max.dist}, \code{min.height}, \code{max.height}, \code{algorithm.dtm}, \code{res.dtm}, \code{csf}, \code{RGB}: same description and format as indicated for same named \code{\link{normalize}} arguments.}


  \item{Tree detection:}

  \itemize{

    \item \code{dbh.min}, \code{dbh.max}, \code{h.min}, \code{ncr.threshold}, \code{tls.resolution}, \code{tls.precision}, \code{breaks}, \code{stem.section}, \code{d.top}, \code{plot.attributes}: same description and format as indicated for same named \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}} arguments}

  \item{dir.data}{Optional character string naming the absolute path of the directory where LAS files containing TLS point clouds are located. \code{.Platform$file.sep} must be used as the path separator in \code{dir.dat}. If this argument is not specified by the user, it will be set to \code{NULL} by default and, as consequence, the current working directory of the \R process will be assigned to \code{dir.dat} during the execution.}

  \item{save.result}{Optional logical which indicates whether or not the output files described in \sQuote{Output Files} section should be saved in the \code{dir.result}. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as a consequence, the output files will be saved.}

  \item{dir.result}{Optional character string naming the absolute path of an existing directory where files described in \sQuote{Output Files} section will be saved. \code{.Platform$file.sep} must be used as the path separator in \code{dir.result}. If this argument is not specified by the user, and \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.result} during the execution.}}


\details{

See \code{\link{normalize}}, \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}} for further details.}

\value{Data frame with the same description and format as \code{\link{tree.detection.single.scan}} and \code{\link{tree.detection.multi.scan}} \sQuote{Values}. In this case, the \code{id} of plots will be encoded with correlative numbers from 1 to n, where n is the number of LAS files included in \code{files} argument, and \code{file} column will be encoded as \code{id}, but including .las extension.}


\section{Output Files}{

At the end of the tree detection process, if the \code{save.result} argument is \code{TRUE}, the function will save both, the reduced point clouds as TXT files encoded according to \code{file} column of \sQuote{Value}; and the data frame with the tree list described in \sQuote{Value} as CSV file (see \code{\link{normalize}} and \code{\link{tree.detection}} \sQuote{Output files}). All outputs are written without row names in the \code{dir.result} directory using \code{\link{vroom_write}} function from \CRANpkg{vroom} package.}


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

download.file("https://www.dropbox.com/s/iaf4sf9z6w4ga14/Pinus_radiata_multi_scan.laz?dl=1",
              destfile = file.path(dir.data, "Pinus_radiata_multi_scan.laz"),
              mode = "wb")

download.file("https://www.dropbox.com/s/ypar7fwqcht8aw6/Fagus_sylvatica_SLAM.laz?dl=1",
              destfile = file.path(dir.data, "Fagus_sylvatica_SLAM.laz"),
              mode = "wb")


# Obtaining a vector with names of LAZ files located in dir.data

files <- list.files(pattern = "laz$", path = dir.data)


# Tree detection

tree.tls <- tree.detection.several.plots(las.list = files[1:2],

                                         id = c("Pinus_radiata_multi_scan",
                                                "Fagus_sylvatica_SLAM"),

                                         scan.approach = "multi",

                                         breaks = 1.3, stem.section = c(1,5))

  }
}

