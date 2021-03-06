% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bound.R
\name{bound_precision}
\alias{bound_precision}
\title{Bound Precision}
\usage{
bound_precision(vals)
}
\arguments{
\item{vals}{\code{numeric} vector of values in the interval [0, 1] to be
bounded within arbitrary machine precision. The most common use of this
functionality is to avoid indeterminate or non-finite values after the
application \code{stats::qlogis}.}
}
\value{
A \code{numeric} vector of the same length as \code{vals}, where
 the returned values are bounded to machine precision. This is intended to
 avoid numerical instability issues.
}
\description{
Bound Precision
}
\details{
Bound values in the unit interval to machine precision in order to
 avoid numerical instability issues in downstream computation.
}
