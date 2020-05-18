\name{proteinMass}
\alias{proteinMass}
\alias{peptideMZ}
\title{
  Mass and m/z of modified proteins and peptides.
  }
\description{
  \code{proteinMass} calculates the correct mass of modified proteins and peptides.
  \code{peptideMZ} calculates the mass over charge ratio (m/z) of modified proteins and peptides.
  This function automatically uses the monoisotopic mass.
  You can use predefined heavy isotope labels and/or define custom modifications.
  }
\usage{
  proteinMass(sequence, monoisotopic = FALSE, label = "none")
  peptideMZ(sequence, charge = 2, label = "none")
  }
\arguments{
  \item{sequence}{
    The amino acid sequence of the peptide or protein, in one letter code. No FASTA header allowed.
    }
  \item{label}{
    Set a predefined heavy isotope label. Accepts "none", "silac_13c", "silac_13c15n" and "15n".
    }
  \item{monoisotopic}{
    If FALSE, the default, the average mass of a protein is calculated. If TRUE, the monoisotopic mass is calculated.
    }
  \item{charge}{
    The charge of a peptide in a LC-MS/MS measurement (defaults to 2).
    }
}
\references{
  The unmodified protein mass is calculated using amino acid weights as defined by \hyperref{Expasy}{http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html}.
}
\seealso{
  For the exact definition of predefined heavy isotope labels and the use of custom modifications, see \code{\link{massShift}}.
}
\examples{
  proteinMass("EGVNDNEEGFFSAR", label = "silac_13c")
  peptideMass("EGVNDNEEGFFSAR", label = "silac_13c")
}
\keyword{ genetics }
\keyword{ protein }
\keyword{ peptide }
