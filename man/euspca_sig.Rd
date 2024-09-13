\name{euspca_sig}
\alias{euspca_sig}
\title{solves EUSPCA with \eqn{\ell_1}{l1} regularization when an empirical covariance/correlation matrix is given}
\description{Solving EUSPCA with \eqn{\ell_1}{l1} regularization when an empirical covariance/correlation matrix is explictly given.}
\usage{
euspca_sig(Sig, k, scale, lambda,
  eps1 = 1e-4, eps2 = 1e-5, eps.sub = 1e-4,
  max.iter.outer = 1000, max.iter.inner = 2000, track=5,
  parm.outer = list(sig = 0.2, tau = 1.1), 
  parm.inner = list(eta = 10, gam = 1e-4, M = 5, beta.min = 1e-15, beta.max = 1e+100) )
}  
\arguments{
	\item{Sig}{An empirical covariance/correlation matrix.}
	\item{k}{The number of uncorrelated principal components.}
  \item{lambda}{A non-negative scalar. It is a value of tuning parameter \eqn{\lambda}. }
	\item{scale}{Perform EUSPCA on the empirical correlation matrix? If \code{scale=FALSE}, the empirical covariance matrix is used. Default is \code{FALSE}.}
	\item{eps1}{An accuracy parameter for the objective function.}
	\item{eps2}{An accuracy parameter for the uncorrelated constraint.}
	\item{eps.sub}{An accuracy parameter for the objective function of the subproblem.}
	\item{max.iter.outer}{The maximum number of iterations of the outer algorithm.}
	\item{max.iter.inner}{The maximum number of iterations of the inner algorithm.}
	\item{parm.outer}{A list with parameters in the outer algorithm, \eqn{\sigma} and \eqn{\tau}. Default is \code{list(sig = 0.2, tau = 1.1)}.}
	\item{parm.inner}{A list with parameters in the inner algorithm, \eqn{\eta}, \eqn{\gamma}, \eqn{M}, \eqn{\beta_{min}}{\beta_min}, and \eqn{\beta_{max}}{\beta_
	max}. Default is \code{list(eta = 10, gam = 10^-4, M = 5, beta.min = 10^-15, beta.max = 10^100)}.}
	\item{track}{Prints out its progress every \code{track} iterations. Default is \code{track=5}. If the user does not want to prinit out the progress, use \code{track=NULL}.}
}

\value{
An ``euspca" object is returned. The below are some quantities which the user may be interested in: \cr
  \item{V}{The solution of EUSPCA with \eqn{\ell_1}{l1} regularization.}
  \item{loadings}{The loadings of the unccorrelated sparse PCs. \eqn{\ell_2}{l2} normalize the rows of \code{V}. }
  \item{p.nz}{The percentage of non-zero components in \code{V} or \code{loadings}.}  
  \item{p.ev}{The percentage of explained adjusted total variance.}
  \item{pc.cor}{The correlation matrix of uncorrelated sparse components.}
}

\seealso{\code{\link{euspca}} and \code{\link{print.euspca}}}
