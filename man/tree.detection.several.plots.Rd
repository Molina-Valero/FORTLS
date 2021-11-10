\encoding{UTF-8}

\name{tree.detection.several.plots}

\alias{tree.detection.several.plots}

\title{Tree Detection and Cross Section Estimation for Several Plots}

\description{This function integrates both, the \code{\link{normalize}} and \code{\link{tree.detection.single.scan}} or \code{\link{tree.detection.multiple.scans}} functions, generating the same \sQuote{Output Files} as indicated for these, and it returs the same \sQuote{Value} as described for \code{\link{tree.detection}}. However, this function is designed for working with several plots, producing a list of all scans considered automatically from LAS files.}

\usage{
tree.detection.multiple(las.list, id = NULL, file = NULL,
                        normalize.arguments = list(max.dist = NULL,
                                                   min.height = NULL,
                                                   max.height = NULL,
                                                   algorithm.dtm = "tin",
                                                   res.dtm = 0.2,
                                                   multiple.scans = NULL),
                        tree.detection.arguments = list(dbh.min = 7.5,
                                                        dbh.max = 200,
                                                        ncr.threshold = 0.1,
                                                        tls.resolution = list(),
                                                        breaks = c(1.0, 1.3, 1.6),
                                                        plot.attributes = NULL),
                        dir.data = NULL, save.result = TRUE, dir.result = NULL)}

\arguments{

  \item{las.list}{Character vector containing the names of all LAS files for analysis and belonging to TLS point cloud, including .las extension (see \sQuote{Examples})}

  \item{id}{Optional vector with plots identification encoded as character string or numeric. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the plots will be encoded with correlative numbers from 1 to n plots.}

  \item{file}{Optional vector containing files name identification encoded as character string or numeric value. If it is null, file will be encoded as \code{id} by default.}

  \item{normalize.arguments}{Optional list including the following arguments related to the normalization process (see \code{\link{normalize}}):}

  \itemize{

    \item \code{max.dist}, \code{min.height}, \code{max.height}, \code{algorithm.dtm}, \code{res.dtm}: same description and format as indicated for same named \code{\link{normalize}} arguments, \code{multiple.scans}: same description and format as indicated for same named \code{\link{normalize}} arguments}


  \item{tree.detection.arguments}{Optional list including the following arguments related to tree detection process (see \code{\link{tree.detection}}):}

  \itemize{

    \item \code{dbh.min}, \code{dbh.max}, \code{ncr.threshold}, \code{tls.resolution}, \code{breaks}, \code{plot.attributes}: same description and format as indicated for same named \code{\link{tree.detection}} arguments}

  \item{dir.data}{Optional character string naming the absolute path of the directory where LAS files containing TLS point clouds are located. \code{.Platform$file.sep} must be used as the path separator in \code{dir.dat}. If this argument is not specified by the user, it will be set to \code{NULL} by default and, as consequence, the current working directory of the \R process will be assigned to \code{dir.dat} during the execution.}

  \item{save.result}{Optional logical which indicates whether or not the output files described in \sQuote{Output Files} section should be saved in the \code{dir.result}. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as a consequence, the output files will be saved.}

  \item{dir.result}{Optional character string naming the absolute path of an existing directory where files described in \sQuote{Output Files} section will be saved. \code{.Platform$file.sep} must be used as the path separator in \code{dir.result}. If this argument is not specified by the user, and \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.result} during the execution.}}


\details{

See \code{\link{normalize}} and \code{\link{tree.detection}} for further details.}

\value{Data frame with the same description and format as \code{\link{tree.detection}} \sQuote{Value}. In this case, the \code{id} of plots will be encoded with correlative numbers from 1 to n, where n is the number of LAS files included in \code{files} argument, and \code{file} column will be encoded as \code{id}, but including .las extension.}


\section{Output Files}{

At the end of the tree detection process, if the \code{save.result} argument is \code{TRUE}, the function will save both, the reduced point clouds as TXT files encoded according to \code{file} column of \sQuote{Value}; and the data frame with the tree list described in \sQuote{Value} as CSV file (see \code{\link{normalize}} and \code{\link{tree.detection}} \sQuote{Output files}). All outputs are written without row names in the \code{dir.result} directory using \code{\link{vroom_write}} function from \CRANpkg{vroom} package.}


\note{

This function has been developed for working with several plots, which will be the most common situation in forest inventory approaches. Nevertheless, several LAS files are not provided as examples due to problems with memory capacity.}


\author{
Juan Alberto Molina-Valero,
María José Ginzo Villamayor,
Manuel Antonio Novo Pérez,
Adela Martínez-Calvo,
Juan Gabriel Álvarez-González,
Fernando Montes and
César Pérez-Cruzado.}

\seealso{
\code{\link{normalize}}, \code{\link{tree.detection}}, \code{\link{distance.sampling}}, \code{\link{estimation.plot.size}}, \code{\link{simulations}},
\code{\link{metrics.variables}}.}

\examples{

\donttest{

# Establishment of working directories (optional)
# By default here we propose the current working directory of the R process

dir.data <- getwd()
dir.result <- getwd()

# Loading example data (LAS files) to dir.data

download.file("https://www.dropbox.com/s/2c3d320o3srcawb/1.las?raw=1",
              destfile = file.path(dir.data, "1.las"), mode = "wb")

download.file("https://www.dropbox.com/s/9k8zn5dt0xcxfof/2.las?raw=1",
              destfile = file.path(dir.data, "2.las"), mode = "wb")


# Obtaining a vector with names of LAS files located in dir.data

files <- list.files(pattern = "las$", path = dir.data)


# Tree detection

tree.list.tls <- tree.detection.multiple(las.list = files,
                                         normalize.arguments = list(max.dist = 15,
                                                                    min.height = 0.25,
                                                                    max.height = 25),
                                         tree.detection.arguments = list(dbh.min = 7.5,
                                                                         dbh.max = 50,
                                                                         breaks = 1.3,
                                         tls.resolution = list(point.dist = 7.67,
                                                               tls.dist = 10)),
                                         dir.data = dir.data, dir.result = dir.result)



}}

