\encoding{UTF-8}

\name{tree.detection}

\alias{tree.detection}

\title{Tree Detection and Cross Section Estimation}

\description{Detects trees from TLS point clouds corresponding to a single scan. For each tree detected, the function calculates the central coordinates and estimates the diameter at 1.3 m above ground level (which is known as \emph{dbh}, diameter at breast height) and classifies it as fully visible or partially occluded. Finally, the function obtains the number of points belonging to normal sections of trees (those corresponding to \emph{dbh} +/- 5 cm) and estimates them for both original and reduced (with point cropping process) point clouds.}

\usage{
tree.detection(data,
               dbh.min = 7.5, dbh.max = 200,
               ncr.threshold = 0.1,
               tls.resolution = list(),
               breaks = c(1.0, 1.3, 1.6),
               plot.attributes = NULL,
               dir.result = NULL, save.result = TRUE)}

\arguments{

  \item{data}{Data frame with same description and format as indicated for \code{\link{normalize}} \sQuote{Value}.}

  \item{dbh.min}{Optional minimum \emph{dbh} (cm) considered for detecting trees. By default it will be set at 7.5 cm.}

  \item{dbh.max}{Optional maximum \emph{dbh} (cm) considered for detecting trees. By default it will be set at 200 cm.}

  \item{ncr.threshold}{Local surface variation (also known as normal change rate, NCR). By default it will be set as 0.1. For better understanding of this argument see \sQuote{Details}.}

  \item{tls.resolution}{List containing parameters of TLS resolution. This can be defined by the angle aperture:}

\itemize{

  \item {\code{horizontal.angle}}: horizontal angle aperture (degrees).

  \item {\code{vertical.angle}}: vertical angle aperture (degrees).}

  \item{}{or separation between two consecutive points at a certain distance from TLS:}

  \itemize{

  \item {\code{point.dist}}: distance (mm) between two consecutive points.

  \item {\code{tls.dist}}: distance (m) from TLS at which two consecutive points are separated by \code{point.dist}.

  If this argument is not specified by the user, it will be set to NULL by default and, as a consequence the function will stop giving an error message.}

  \item{breaks}{Height above ground level (m) of slices considered for detecting trees. By default it will be 1.0, 1.3 and 1.6 m (+/- 5 cm).}

  \item{plot.attributes}{Data frame with attributes at plot level. It must contain a column named \code{id} (character string or numeric value) with encoding coinciding with that used in \code{id} argument of \code{\link{normalize}} for identifying plots. If there are strata, another column named \sQuote{stratum} (numeric) will be required for other functionalities of \pkg{FORTLS} (see, for instance, \code{\link{estimation.plot.size}} or \code{\link{metrics.variables}}). If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the function will not add these possible plot attributes.}

  \item{save.result}{Optional logical which indicates whether or not the output files described in \sQuote{Output Files} section should be saved in \code{dir.result}. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as a consequence, the output files will be saved.}

  \item{dir.result}{Optional character string naming the absolute path of an existing directory where the files described in \sQuote{Output Files} section will be saved. \code{.Platform$file.sep} must be used as the path separator in \code{dir.result}. If this argument is not specified by the user, and the \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.result} during the execution.}}

\details{

Slices determined by \code{breaks} argument are clustered using the DBSCAN algorithm (Ester et al., 1996) on the horizontal plane according to Cartesian coordinates (x, y). Before and after this process, several algorithms are used to remove noisy points and apply classification criteria to select the clusters of trees.

\emph{dbh} is directly estimated for the section of 1.3 m above ground level, and estimated from other sections using \emph{dbh}~\emph{breaks} linear regression. Finally, the mean value of all estimates is provided in \sQuote{Value} as the \emph{dbh} of the tree section.

The number of points corresponding to a normal section (+/- 5 cm) is estimated in proportion to \emph{dbh}, using the average number of points per radius unit as reference. In this respect, only tree sections fully visible at 1.3 m above ground level will be considered for estimating the average number of points.

Local surface variation (also known as normal change rate ,NCR), is a quantitative measure of curvature feature (Pauly et al., 2002). This is useful for distinguishing points belonging to fine branches and foliage (e.g. leaves, shrubs) and stem points (e.g. Jin et al., 2016; Zhang et al., 2019). Just as we considered 5 cm as suitable for calculating local surface variation for the stem separation in forests, according to other authors (Ma et al., 2015; Xia et al., 2015), we also established the NCR threshold as 0.1, according to Zhang et al. (2019). However, this argument (\code{ncr.threshold}) may be modified in order to use more appropriate values.}

\value{Data frame with the following columns for every tree detected (each row corresponds to one tree detected):

  \item{id}{Optional plot identification encoded as a character string or numeric. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the plot will be encoded as 1.}

  \item{file}{Optional file name identification encoded as character string or numeric. If it is null, the file will be encoded as \code{id} by default.}

  \item{tree}{tree numbering}


\item{Coordinates}{

  Cartesian (according to \url{https://en.wikipedia.org/wiki/Cartesian_coordinate_system} notation):}
    \itemize{
      \item \code{x}: distance on x axis (m) of tree centre.
      \item \code{y}: distance on y axis (m) of tree centre.}

  Azimuthal angles:
    \itemize{
      \item \code{phi}: angular coordinate (rad) of tree centre.
      \item \code{phi.left}: angular coordinate (rad) of left border of tree section.
      \item \code{phi.right}: angular coordinate (rad) of right border of tree section}

\item{horizontal.distance}{horizontal distance (m) from plot centre to tree centre.}

\item{\emph{dbh}}{estimated tree diameter (cm) at breast height (1.3 m).}

\item{num.points}{number of points corresponding to a normal section (+/- 5 cm) in the original point cloud.}

\item{num.points.hom}{number of points corresponding to a normal section (+/- 5 cm) in the point cloud reduced by the point cropping process.}

\item{num.points.est}{number of points estimated for a normal section (+/- 5 cm) in the original point cloud.}

\item{num.points.hom.est}{number of points estimated for a normal section (+/- 5 cm) in the point cloud reduced by the point cropping process.}

\item{partial.occlusion}{yes (1) or no (0)}}

\section{Output Files}{

At the end of the tree detection process, if the \code{save.result} argument is \code{TRUE}, the function will save the data frame described in \sQuote{Value} as a CSV file named \sQuote{tree.list.tls.csv}. The data frame will be written without row names in the \code{dir.result} directory by using \code{\link[utils]{write.csv}} function from the \pkg{utils} package.}


\note{
Although \code{tree.detection} also works with reduced point clouds, thus reducing the computing time, use of the original point cloud is recommended in order to detect more trees. This will also depend on forest conditions, especially those related to visibility. The more distant the trees are, the lower the density of points will be, and using reduced point clouds will therefore complicate detection of the most distant trees.

Note that \code{dbh.min} and \code{dbh.max} are important for avoiding outlier values when inventory data are used for reference purposes. Otherwise, knowledge about the autoecology of species could be used for filtering anomalous values of \emph{dbh}.

The argument \code{breaks = 1.3} could be sufficient for detecting trees visible at \emph{dbh}, involving lower computational cost. However, those trees not detected at \emph{dbh}, may be estimated from lower and/or higher sections. Considering the three default sections in the argument \code{breaks = c(1.0, 1.3, 1.6)} maintains a good balance in the case study of this package.}


\references{
Ester, M., Kriegel, H. P., Sander, J., & Xu, X. (1996). A density-based algorithm for discovering clusters in large spatial databases with noise. In Kdd (Vol. 96, No. 34, pp. 226-231).

Jin, S., Tamura, M., & Susaki, J. (2016). A new approach to retrieve leaf normal distribution using terrestrial laser scanners. J. \emph{Journal of Forestry Research}, \bold{27(3)}, 631-638. \doi{https://doi.org/10.1007/s11676-015-0204-z}

Ma, L., Zheng, G., Eitel, J. U., Moskal, L. M., He, W., & Huang, H. (2015). Improved salient feature-based approach for automatically separating photosynthetic and nonphotosynthetic components within terrestrial lidar point cloud data of forest canopies. \emph{IEEE Transactions Geoscience Remote Sensing}, \bold{54(2)}, 679-696. \doi{https://doi.org/10.1109/TGRS.2015.2459716}

Pauly, M., Gross, M., & Kobbelt, L. P., (2002). Efficient simplification of point-sampled surfaces. In IEEE Conference on Visualization. (pp. 163-170). Boston, USA. \doi{https://doi.org/10.1109/VISUAL.2002.1183771}

Xia, S., Wang, C., Pan, F., Xi, X., Zeng, H., & Liu, H. (2015). Detecting stems in dense and homogeneous forest using single-scan TLS. \emph{Forests}. \bold{6(11)}, 3923-3945. \doi{https://doi.org/10.3390/f6113923}

Zhang, W., Wan, P., Wang, T., Cai, S., Chen, Y., Jin, X., & Yan, G. (2019). A novel approach for the detection of standing tree stems from plot-level terrestrial laser scanning data. \emph{Remote Sens}. \bold{11(2)}, 211. \doi{https://doi.org/10.3390/rs11020211}}

\author{
Juan Alberto Molina-Valero,
María José Ginzo Villamayor,
Manuel Antonio Novo Pérez,
Adela Martínez-Calvo,
Juan Gabriel Álvarez-González,
Fernando Montes and
César Pérez-Cruzado.}

\seealso{
\code{\link{normalize}}, \code{\link{tree.detection.multiple}}, \code{\link{distance.sampling}}, \code{\link{estimation.plot.size}}, \code{\link{simulations}}, \code{\link{metrics.variables}}
}

\examples{

\donttest{

# Establishment of working directories (optional)
# By default here we propose the current working directory of the R process

dir.data <- getwd()
dir.result <- getwd()


# Loading example data (LAS file) to dir.data

download.file("https://www.dropbox.com/s/2c3d320o3srcawb/1.las?raw=1",
              destfile = file.path(dir.data, "1.las"),
              method = "wininet", mode = "wb")

# Normalizing the point cloud data as a neccesary step for detecting trees

pcd <- normalize(las = "1.las",
                 max.dist = 15, min.height = 0.25, max.height = 25,
                 id = "1", file = "1.txt",
                 dir.data = dir.data, dir.result = dir.result)


# Tree detection without considering arguments
# For this case study, TLS resolution was established as:
# point.dist = 7.67 mm and tls.dist = 10 m

tree.list.tls <- tree.detection(data = pcd,
                                tls.resolution = list(point.dist = 7.67, tls.dist = 10),
                                dir.result = dir.result)



# Tree detection considering several arguments

tree.list.tls <- tree.detection(data = pcd,
                                dbh.min = 7.5, dbh.max = 60,
                                tls.resolution = list(point.dist = 7.67, tls.dist = 10),
                                breaks = 1.3,
                                dir.result = dir.result)

  }

}

