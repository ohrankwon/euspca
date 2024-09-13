\name{euspca}
\alias{euspca}
\title{Exactly uncorrelated sparse principal component analysis (EUSPCA) with \eqn{\ell_1}{l1} regularization}
\description{Finding uncorrelated sparse principal components by solving EUSPCA with \eqn{\ell_1}{l1} regularization using a combination of augmented Lagrangian and non-monotone proximal gradient methods.}
\usage{
euspca(X, is.data.mat, k, lambda, scale=FALSE, \dots)
}
\arguments{
	\item{X}{A matrix. It can be a data matrix or an empirical covariance/correlation matrix.}
	\item{is.data.mat}{Is \code{X} a data matrix? If \code{is.data.mat=TRUE}, \code{X} is a data matrix. If \code{is.data.mat=FALSE}, \code{X} is a cov/cor matrix.}
	\item{k}{The number of uncorrelated principal components.}
  \item{lambda}{A non-negative scalar. It is a value of tuning parameter \eqn{\lambda}. }
	\item{scale}{Perform EUSPCA on the empirical correlation matrix? If \code{scale=FALSE}, the empirical covariance matrix is used. Default is \code{FALSE}.}
	\item{\dots}{Other arguments being passed to \code{\link{euspca_sig}} or \code{\link{euspca_dat}}.}
}

\details{The function \code{euspca} solves EUSPCA with \eqn{\ell_1}{l1} regularization, that is to find a minimizer \eqn{\mathbf V}{V} of
\deqn{- \mbox{tr} (\mathbf V \mathbf
        \Sigma_n^2 \mathbf V^T ) + \lambda \sum_{ij} | v_{ij} |
        ~~\mbox{ s.t.} ~  \mathbf V \mathbf \Sigma_n \mathbf V^T = \mathbf{I},}{ -tr (V \Sigma_n^2 V^T) + \lambda \sum |v_ij| s.t. V \Sigma_n V^T = I,}
where \eqn{\Sigma_n}{\Sigma_n} is an empirical covariance/correlation matrix and \eqn{v_{ij}}{v_ij} is the (i,j)th component of \eqn{\mathbf V}{V}. The user sepcifies a value of \eqn{\lambda}, i.e., \code{lambda}. 
}

\value{
An ``euspca" object is returned. The below are some quantities which the user may be interested in: \cr
  \item{V}{The solution of EUSPCA with \eqn{\ell_1}{l1} regularization.}
  \item{loadings}{The loadings of the unccorrelated sparse PCs. \eqn{\ell_2}{l2} normalize the rows of \code{V}. }
  \item{p.nz}{The percentage of non-zero components in \code{V} or \code{loadings}.}  
  \item{p.ev}{The percentage of explained adjusted total variance.}
  \item{pc.cor}{The correlation matrix of uncorrelated sparse components.}
}

\seealso{\code{\link{euspca_sig}}, \code{\link{euspca_dat}}, and \code{\link{print.euspca}}.}

\examples{
data(syn)
# fit a euspca on covariance matrix
euspca_syn = euspca(syn, is.data.mat=FALSE, k=2, lamb=1000, scale=FALSE)
euspca_syn
}

\keyword{euspca}
\keyword{principal component analysis}
\keyword{sparse principal component}
\keyword{uncorrelated component}
\keyword{package}
