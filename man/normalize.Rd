\encoding{UTF-8}

\name{normalize}

\alias{normalize}

\title{Relative Coordinates and Density Reduction for Terrestrial-Based Technologies Point Clouds}

\description{This function obtains coordinates relative to the plot centre speciefied for Terrestrial Laser Scanner (TLS) and Mobile Laser Scanner (MLS) point clouds (supplied as LAS or LAZ files). Point clouds obtained from other devices/approaches (e.g. photogrammetry) can be also used, but the guarantee of good performance is likely to be lower. In addition, the point cropping process developed by Molina-Valero et al., (2019) is applied as a criterion for reducing point density homogeneously in space and proportionally to object size when TLS single-scans are provided.}

\usage{
normalize(las, normalized = NULL,
          x.center = NULL, y.center = NULL,
          x.side = NULL, y.side = NULL,
          xpoly = NULL, ypoly = NULL,
          max.dist = NULL, min.height = NULL, max.height = 50,
          algorithm.dtm = "knnidw", res.dtm = 0.2,
          csf = list(cloth_resolution = 0.5),
          intensity = NULL, RGB = NULL,
          scan.approach = "single",
          voxel_size = NULL,
          id = NULL, file = NULL, plot = TRUE,
          dir.data = NULL, save.result = TRUE, dir.result = NULL,
          save.las = NULL)}

\arguments{

  \item{las}{Character string containing the name of LAS/LAZ file belonging to point cloud, including .las/.laz extensions (see \sQuote{Examples}). Planimetric coordinates of point cloud data must be in local, representing TLS scan point the origin with Cartesian coordinates x and y as (0, 0).}

  \item{normalized}{Optional argument to establish as \code{TRUE} when point cloud is already normalized.}

  \item{x.center}{Planimetric x center coordinate of point cloud data.}

  \item{y.center}{Planimetric y center coordinate of point cloud data.}

  \item{x.side}{x-side (m) of the plot when the plot is square or rectangular.}

  \item{y.side}{y-side (m) of the plot when the plot is square or rectangular.}

  \item{xpoly}{Vector containing the x-coordinates of the vertices of the clipping polygon.}

  \item{ypoly}{Vector containing the y-coordinates of the vertices of the clipping polygon.}

  \item{max.dist}{Optional maximum horizontal distance (m) considered from the plot centre. All points farther than \code{max.dist} will be discarded after the normalization process. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, all points will be used in processing, with \code{max.dist} representing the farthest point.}

  \item{min.height}{Optional minimum height (m) considered from ground level. All points below \code{min.height} will be discarded after the normalization process. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, all points will be used in processing, with \code{min.height} representing the lowest point.}

  \item{max.height}{Optional maximum height (m) considered from ground level. All points above \code{max.height} will be discarded after the normalization process. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, all points will be used in processing, with \code{max.height} representing the highest point.}

    \item{algorithm.dtm}{Algorithm used to generate the digital terrain model (DTM) from the TLS point cloud. There are two posible options based on spatial interpolation: \sQuote{tin} and \sQuote{knnidw} (see \sQuote{Details}). If this argument is not specified by the user, it will be set to \sQuote{knnidw} algorithm.}

      \item{res.dtm}{Numeric parameter. Resolution of the DTM generated to normalize point cloud (see \sQuote{Details}). If this argument is not specified by the user, it will be set to 0.2 m.}

      \item{csf}{List containing parameters of CSF algorithm:}

\itemize{

  \item {\code{cloth_resolution}}: by default 0.5.}

  \item{intensity}{Logical parameter useful when point clouds have LiDAR intesinty values.}

  \item{RGB}{Logical parameter useful when point clouds are colorized, thus including values of RGB colors. It is based on the Green Leaf Algorithm (GLA) (see \sQuote{Details}).}

  \item{scan.approach}{Character parameter indicating TLS single-scan (\sQuote{single}) or TLS multi-scan approach or SLAM point clouds (\sQuote{multi}) approaches. If this argument is not specified by the user, it will be set to \sQuote{multi} approach.}

  \item{voxel_size}{Defines the size of the 3D grid cells used for downsampling.}

  \item{id}{Optional plot identification encoded as character string or numeric. If this argument is not specified by the user, it will be set to NULL by default and, as a consequence, the plot will be encoded as 1.}

  \item{file}{Optional file name identification encoded as character string or numeric value. If it is null, file will be encoded as \code{id} by default.}

  \item{plot}{Optional logical which indicates whether or not the normalized point cloud will be plot. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as consequence, the normalized point cloud will be plot.}

  \item{dir.data}{Optional character string naming the absolute path of the directory where LAS files containing TLS point clouds are located. \code{.Platform$file.sep} must be used as the path separator in \code{dir.data}. If this argument is not specified by the user, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.data} during the execution.}

  \item{save.result}{Optional logical which indicates whether or not the output files described in \sQuote{Output Files} section must be saved in \code{dir.result}. If this argument is not specified by the user, it will be set to \code{TRUE} by default and, as consequence, the output files will be saved.}

  \item{dir.result}{Optional character string naming the absolute path of an existing directory where files described in \sQuote{Output Files} section will be saved. \code{.Platform$file.sep} must be used as the path separator in \code{dir.result}. If this argument is not specified by the user, and \code{save.result} is \code{TRUE}, it will be set to \code{NULL} by default and, as a consequence, the current working directory of the \R process will be assigned to \code{dir.result} during the execution.}

  \item{save.las}{Optional logical which indicates whether or not the imput point cloud must be saved in \code{dir.result} as LAZ file.}}

\details{

Relative coordinates are obtained by means of a normalization process, generating a digital terrain model (DTM) from the TLS point cloud, with the ground height set at 0 m. The DTM is generated by spatial interpolation of ground points classified with the CSF algorithm (Zhang et al., (2016)). Two algorithms are available for that purpose: (i) spatial interpolation based on a Delaunay triangulation, which performs a linear interpolation within each triangle (\sQuote{tin}); (ii) spatial interpolation using a k-nearest neighbour (KNN) approach with inverse-distance weighting (IDW) (\sQuote{knnidw}). Note that normalization process is based on \CRANpkg{lidR} package functions: \link[lidR]{classify_ground}, \link[lidR]{grid_terrain} and \link[lidR]{normalize_height}.

The point cropping process reduces the point cloud density proportionally to the likelihood that objects will receive points according to their distance from TLS and their size, which is determined by angle aperture (the farther they are, the lower the density is). The result is an approximately homogeneous point cloud in three-dimensional space (for more details see Molina-Valero et al., (2019)).

The Green Leaf Algorithm (GLA) is calculated according to Louhaichi et al., (2001)as follows:

Those points with values below 0 are clasified as woody parts, thus retained for tree detection in further functions.}

\value{Data frame of normalized point cloud including the following columns (each row corresponds to one point):

\item{id}{Plot identification encoded as a character string or numeric in the argument \code{id}.}

\item{file}{File name identification encoded as character string or numeric, corresponding to the normalized and reduced point clouds saved. This coincides with the TXT file in the absolute path specified in \code{dir.result} (if \code{save.result} is set to \code{TRUE}).}

\item{Coordinates}{

  Cartesian (according to \url{https://en.wikipedia.org/wiki/Cartesian_coordinate_system} notation):
    \itemize{
      \item \code{x}: x axis distance (m).
      \item \code{y}: y axis distance (m).
      \item \code{z}: height (m).
    }

  Cylindrical (according to \url{https://en.wikipedia.org/wiki/Cylindrical_coordinate_system} notation):
    \itemize{
      \item \code{rho}: horizontal distance (m).
      \item \code{phi}: angular coordinate (rad).
      \item \code{z}: height (m).
    }

  Spherical (according to \url{https://en.wikipedia.org/wiki/Spherical_coordinate_system} notation):
    \itemize{
      \item \code{r}: radial distance (m).
      \item \code{theta}: polar angle (rad).
      \item \code{phi}: azimuthal angle (rad)
    }}

\item{slope}{Slope of the terrain (rad).}

\item{intensity}{Intensity (only if point cloud has intensity values and specified in arguments).}

\item{R}{Red (only if point cloud is colorized and specified in arguments).}

\item{G}{Green (only if point cloud is colorized and specified in arguments).}

\item{B}{Blue (only if point cloud is colorized and specified in arguments).}

\item{GLA}{Green Leaf Algorithm (only if point cloud is colorized and specified in arguments).}

\item{prob}{selection probability assigned in point cropping process (0-1]. Only the farthest will have probability of 1.}

\item{prob.select}{final selection probability assigned in point cropping process. Selected (1) and discarded point (0).}}

\section{Output Files}{

At the end of the normalization process, if the \code{save.result} argument is \code{TRUE}, the function will save the reduced point cloud as TXT file and encoded according to \code{file} argument. The format is the same as data frame described in \sQuote{Value}, except for a \code{prob} column, which is removed because all points selected after the point cropping process have a final selection probability of 1. The data frame is written without row names in \code{dir.result} directory using the \code{\link[vroom]{vroom_write}} function in the \CRANpkg{vroom} package.}

\note{
Note that \code{max.dist}, \code{min.height} and \code{max.height} arguments may be useful for optimizing computing time as well as for removing unnecessary and/or outlier points. These values may be selected more appropriately when inventory data are already available, or the user has some knowledge about autoecology of scanned tree species.

Note also that the linear interpolation algorithm (\sQuote{tin} in this package) showed the highest accuracy in Liang et al., (2018) in DTM generation with single-scans. In this work a DTM resolution of 0.2 m was also considered adequately for square plots of 32 x 32 m.}

\references{
Liang, X., Hyyppä, J., Kaartinen, H., Lehtomäki, M., Pyörälä, J., Pfeifer, N., ... & Wang, Y. (2018). International benchmarking of terrestrial laser scanning approaches for forest inventories. \emph{ISPRS journal of photogrammetry and remote sensing}, \bold{144}, 137-179. \doi{https://doi.org/10.1016/j.isprsjprs.2018.06.021}

Louhaichi, M., Borman, M. M., & Johnson, D. E. (2001). Spatially located platform and aerial photography for documentation of grazing impacts on wheat. \emph{Geocarto International}, \bold{16(1)}, 65-70. \doi{https://doi.org/10.1080/10106040108542184}

Molina-Valero J. A., Ginzo-Villamayor M. J., Novo Pérez M. A., Álvarez-González J. G., & Pérez-Cruzado C. (2019). Estimación del área basimétrica en masas maduras de Pinus sylvestris en base a una única medición del escáner laser terrestre (TLS). \emph{Cuadernos de la Sociedad Espanola de Ciencias Forestales}, \bold{45(3)}, 97-116. \doi{https://doi.org/10.31167/csecfv0i45.19887}.

Zhang, W., Qi, J., Wan, P., Wang, H., Xie, D., Wang, X., & Yan, G. (2016). An easy-to-use airborne LiDAR data filtering method based on cloth simulation. \emph{Remote Sensing}, \bold{8(6)}, 501. \doi{https://doi.org/10.3390/rs8060501}.}

\seealso{
\code{\link{tree.detection.single.scan}}, \code{\link{tree.detection.multi.scan}}, \code{\link{tree.detection.several.plots}}.
}

\author{
Juan Alberto Molina-Valero and
Adela Martínez-Calvo.}

\examples{

\donttest{

# Establishment of working directories (optional)
# By default here we propose the current working directory of the R process

dir.data <- getwd()
dir.result <- getwd()


# TLS SINGLE-SCAN APPROACH

# Loading example TLS data (LAZ file) to dir.data

download.file("https://www.dropbox.com/s/17yl25pbrapat52/PinusRadiata.laz?dl=1",
              destfile = file.path(dir.data, "PinusRadiata.laz"),
              mode = "wb")


# Normalizing the whole point cloud data without considering arguments

# Note that in the case of TLS single-scans, the coordinates of the plot center
# must be specified by the x.center and y.center arguments, which are
# x.center = 0 and y.center = 0 in this case

pcd <- normalize(las = "PinusRadiata.laz",

                 id = "PinusRadiata",

                 x.center = 0, y.center = 0,

                 dir.data = dir.data, dir.result = dir.result)


# MLS OR TLS MULTI-SCAN APPROACHES

# Loading example MLS data (LAZ file) to dir.data

download.file(
"www.dropbox.com/scl/fi/es5pfj87wj0g6y8414dpo/PiceaAbies.laz?rlkey=ayt21mbndc6i6fyiz2e7z6oap&dl=1",
              destfile = file.path(dir.data, "PiceaAbies.laz"),
              mode = "wb")

# Normalizing the whole point cloud data without considering arguments

pcd <- normalize(las = "PiceaAbies.laz",

                 id = "PiceaAbies",

                 scan.approach = "multi",

                 dir.data = dir.data, dir.result = dir.result)

  }

}
