\encoding{UTF-8}

\name{geometric_features_dist}

\alias{geometric_features_dist}

\title{This function obtains geometric features at point level}

\description{This function obtains geometric features at point level}


\usage{
geometric_features_dist(points,
                        dist,
                        Anisotropy = FALSE,
                        Eigenentropy = FALSE,
                        Eigenvalue_sum = FALSE,
                        First_eigenvalue = FALSE,
                        Linearity = FALSE,
                        Normal_x = FALSE,
                        Normal_y = FALSE,
                        Normal_z = FALSE,
                        Number_of_points = FALSE,
                        Omnivariance = FALSE,
                        PCA_1 = FALSE,
                        PCA_2 = FALSE,
                        Planarity = FALSE,
                        Second_eigenvalue = FALSE,
                        Sphericity = FALSE,
                        Surface_density = FALSE,
                        Surface_variation = FALSE,
                        Third_eigenvalue = FALSE,
                        Verticality = FALSE,
                        Volume_density = FALSE,
                        num_threads = 0,
                        use_spatial_index = TRUE,
                        grid_cell_size = -1.0)}

\arguments{

  \item{points}{points}

  \item{dist}{dist}

  \item{Anisotropy}{Anisotropy}

  \item{Eigenentropy}{Eigenentropy}

  \item{Eigenvalue_sum}{Eigenvalue_sum}

  \item{First_eigenvalue}{First_eigenvalue}

  \item{Linearity}{Linearity}

  \item{Normal_x}{Normal_x}

  \item{Normal_y}{Normal_y}

  \item{Normal_z}{Normal_z}

  \item{Number_of_points}{Number_of_points}

  \item{Omnivariance}{Omnivariance}

  \item{PCA_1}{PCA_1}

  \item{PCA_2}{PCA_2}

  \item{Planarity}{Planarity}

  \item{Second_eigenvalue}{Second_eigenvalue}

  \item{Sphericity}{Sphericity}

  \item{Surface_density}{Surface_density}

  \item{Surface_variation}{Surface_variation}

  \item{Third_eigenvalue}{Third_eigenvalue}

  \item{Verticality}{Verticality}

  \item{Volume_density}{Volume_density}

  \item{num_threads}{num_threads}

  \item{use_spatial_index}{use_spatial_index}

  \item{grid_cell_size}{grid_cell_size}}

\value{Data frame}


\author{
Juan Alberto Molina-Valero}

\keyword{internal}



