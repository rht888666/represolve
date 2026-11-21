\name{list_datasets}
\alias{list_datasets}
\title{List frozen benchmark datasets}
\description{
Reads the bundled \code{DATA_MANIFEST.csv} listing the nine frozen
replicate-half benchmark datasets and SHA-256 checksums archived at Zenodo.
}
\usage{
list_datasets()
}
\value{
A data frame with columns \code{dataset}, \code{file}, \code{bytes},
\code{sha256}, and \code{status}.
}
\examples{
list_datasets()
}
