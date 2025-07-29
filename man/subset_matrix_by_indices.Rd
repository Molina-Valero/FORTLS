\encoding{UTF-8}

\name{subset_matrix_by_indices}

\alias{subset_matrix_by_indices}

\title{This function obtains geometric features at point level}

\description{This function obtains geometric features at point level}


\usage{
subset_matrix_by_indices(mat,
                         index_lists,
                         dist,
                         First_eigenvalue = TRUE,
                         Second_eigenvalue = TRUE,
                         Third_eigenvalue = TRUE,
                         Sum_of_eigenvalues = TRUE,
                         PCA_1 = TRUE,
                         PCA_2 = TRUE,
                         Anisotropy = TRUE,
                         Planarity = TRUE,
                         Linearity = TRUE,
                         Surface_variation = TRUE,
                         Normal_change_rate = TRUE,
                         Verticality = TRUE,
                         Number_of_points = TRUE,
                         omnivariance = TRUE,
                         eigenentropy = TRUE,
                         surface_density = TRUE,
                         volume_density = TRUE,
                         threads = 1,
                         solver_thresh = 50000)}

\arguments{

  \item{mat}{mat}

  \item{index_lists}{index_lists}

  \item{dist}{dist}

  \item{First_eigenvalue}{First_eigenvalue}

  \item{Second_eigenvalue}{Second_eigenvalue}

  \item{Third_eigenvalue}{Third_eigenvalue}

  \item{Sum_of_eigenvalues}{Sum_of_eigenvalues}

  \item{PCA_1}{PCA_1}

  \item{PCA_2}{PCA_2}

  \item{Anisotropy}{Anisotropy}

  \item{Planarity}{Planarity}

  \item{Linearity}{Linearity}

  \item{Surface_variation}{Surface_variation}

  \item{Normal_change_rate}{Normal_change_rate}

  \item{Number_of_points}{Number_of_points}

  \item{omnivariance}{omnivariance}

  \item{eigenentropy}{eigenentropy}

  \item{surface_density}{surface_density}

  \item{volume_density}{volume_density}

  \item{threads}{threads}

  \item{solver_thresh}{solver_thresh}}

\value{Data frame}


\author{
Juan Alberto Molina-Valero}

\keyword{internal}



