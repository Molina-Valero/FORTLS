\encoding{UTF-8}

\name{geometric.features}

\alias{geometric.features}

\title{This function fit a circle based on 3 points}

\description{This function fit a circle based on 3 points}


\usage{
geometric.features(data,
                   grid_method,
                   voxel_resolution = 6.0,
                   data_table_grid_cells = 20,
                   features = c("PCA1",
                                "PCA2",
                                "anisotropy",
                                "planarity",
                                "linearity",
                                "surface_variation",
                                "sphericity",
                                "verticality",
                                "number_neighbors",
                                "omnivariance",
                                "eigenentropy",
                                "surface_density",
                                "volume_density"),
                   dist = 0.05,
                   threads = 1,
                   keep_NaN = FALSE,
                   verbose = FALSE,
                   solver_threshold = 50000)}

\arguments{

  \item{data}{data}

  \item{grid_method}{grid_method}

  \item{voxel_resolution}{voxel_resolution}

  \item{data_table_grid_cells}{data_table_grid_cells}

  \item{features}{features}

  \item{threads}{threads}

  \item{keep_NaN}{keep_NaN}

  \item{verbose}{verbose}

  \item{solver_threshold}{solver_threshold}}

\value{Data frame}


\author{
Juan Alberto Molina-Valero}

\keyword{internal}



