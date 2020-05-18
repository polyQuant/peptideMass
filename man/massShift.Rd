\name{massShift}
\alias{massShift}
\title{
  Calculate the mass difference of modified peptides.
}
\description{
  Calculate the mass difference of peptides introduced by chemical modifications or heavy isotope labelling.
}
\usage{
  massShift(sequence, label = "none", aa, shift)
}
\arguments{
  \item{sequence}{
    The amino acid sequence of the peptide or protein, in one letter code. No FASTA header allowed.
    }
  \item{label}{
    Set a predefined heavy isotope label. Accepts "none", "silac_13c", "silac_13c15n" and "15n". Overwrites input in \code{aa} and \code{shift}.
    }
  \item{aa}{
    Define modified amino acids. Multiple amino acids carrying the same modification can be given as a single string. Multiple amino acids carrying different modifications can be given as a vector.
    }
  \item{shift}{
     Define the mass difference in Dalton of the amino acid(s) given in \code{aa}. If \code{aa} is a vector, \code{shift} must also be a vector of the same length.
    }
}
\references{
  For the predefined heavy isotope labels, compare:
  \item{silac_13c}{\hyperref{Unimod 188}{http://www.unimod.org/modifications_view.php?editid1=188}}
  \item{silac_13c15n}{\hyperref{Unimod 259}{http://www.unimod.org/modifications_view.php?editid1=259} and \hyperref{Unimod 267}C{http://www.unimod.org/modifications_view.php?editid1=267}}
  \item{15n}{\hyperref{Unimod 994}{http://www.unimod.org/modifications_view.php?editid1=994}, \hyperref{Unimod 995}{http://www.unimod.org/modifications_view.php?editid1=995}}, \hyperref{Unimod 996}{http://www.unimod.org/modifications_view.php?editid1=996}} and \hyperref{Unimod 897}{http://www.unimod.org/modifications_view.php?editid1=897}}
}
\examples{
  massShift("EGVNDNEEGFFSAR", label = "silac_13c")
  massShift("EGVNDNEEGFFSAR", aa = "R", shift = 6.020129)

  # To define multiple amino acids with the same modification, set \code{aa} as a single string.
  massShift("EGVNDNEEGFFSARK", aa = "RK", shift = 6.020129)

  # To define multiple amino acids with different modifications, set \code{aa} and \code{shift} as character vectors.
  massShift("EGVNDNEEGFFSAR", aa = c("RK", "C"), shift = c(6.020129, 57.021464))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetics }
\keyword{ proteins }
\keyword{ peptides }
