\name{syn}
\alias{syn}

\title{Covariance matrix of Synthetic example}
\usage{
data(syn)
}

\description{
Covariance matrix of Synethetic example consisting of 9 variables with 3 hidden factors.
}

\details{This \eqn{9 \times 9}{9*9} matrix is the covariance matrix of 9 variables,
\deqn{\xi_i = \eta_1 + \epsilon_i, ~ i=1,2,3,}{\xi_i = \eta_1 + \epsilon_i, i=1,2,3,} 
\deqn{\xi_i = \eta_2 + \epsilon_i, ~ i=4,5,6,}{\xi_i = \eta_2 + \epsilon_i, i=4,5,6,} 
\deqn{\xi_i = \eta_3 + \epsilon_i, ~ i=7,8,9,}{\xi_i = \eta_3 + \epsilon_i, i=7,8,9,}
where \eqn{\epsilon_i, ~ i=1,\ldots,9}{\epsilon_i, i=1,\ldots,9} are independent and \eqn{\eta_i, ~ i=1,2,3}{\eta_i, i=1,2,3} are three hiden factors such that
\deqn{\eta_1 \sim N(0,290), ~ \eta_2 \sim N(0,300), ~ \eta_3 =  0.3 \eta_1 + 0.98 \eta_2 + \epsilon,}{\eta_1 ~ N(0,290), \eta_2 ~ N(0,300), \eta_3 = 0.3 \eta_1 + 0.98 \eta_2 + \epsilon,}
where \eqn{\epsilon \sim N(0,1)}{\epsilon ~ N(0,1)} and \eqn{\eta_1}, \eqn{\eta_2}, and \eqn{\epsilon} are independent.
}

\examples{
# load the covariance matrix
data(syn)

# the dimension of the covariance matrix
dim(syn)
}

\keyword{data set}
