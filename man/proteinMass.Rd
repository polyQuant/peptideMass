\name{proteinMass}
\alias{proteinMass}
\alias{peptideMZ}
\title{
  Mass and m/z of modified proteins and peptides.
  }
\description{
  \code{proteinMass} calculates the correct mass of modified proteins and peptides.
  \code{peptideMZ} calculates the mass over charge ratio (m/z) of modified proteins and peptides.
  }
\usage{
  proteinMass(sequence, label = "none", monoisotopic = FALSE)
  peptideMZ(sequence, label = "none", aa, shift, charge = 2)
  }
\arguments{
  \item{sequence}{
    The amino acid sequence of the peptide or protein, in one letter code. No FASTA header allowed.
    }
  \item{label}{
    Set a predefined heavy isotope label. Accepts "none", "silac_13c", "silac_13c15n" and "15n".
    }
  \item{monoisotopic}{
    If \code{FALSE}, the default, the average mass of a protein is calculated.
    If \code{TRUE}, the monoisotopic mass is calculated.
	\code{peptideMZ} always uses the monoisotopic mass.
    }
  \item{aa}{
	  Define modified amino acids. Multiple amino acids carrying the same modification can be given as a single string.
	  Multiple amino acids carrying different modifications can be given as a vector.
  	}
  \item{shift}{
     Define the mass difference in Dalton of the amino acid(s) given in \code{aa}.
     If \code{aa} is a vector, \code{shift} must also be a vector of the same length.
    }
  \item{charge}{
    The charge of a peptide in a LC-MS/MS measurement (defaults to 2).
    }
  \item{carbamidomethyl}{
    If \code{TRUE}, the default, a carbamidomethlyation at cysteins is added, like introduced by alkylation with iodoacetamide. The exact mass shift is 57.021464 Da at C.
  }
}
\references{
  The unmodified protein mass is calculated using amino acid weights as reported by the \href{https://proteomicsresource.washington.edu/protocols06/masses.php}{UWPR}.
}
\seealso{
  For the exact definition of predefined heavy isotope labels and the use of custom modifications, see \code{\link{massShift}}.
}
\examples{
  proteinMass("EGVNDNEEGFFSAR", label = "silac_13c")
  peptideMass("EGVNDNEEGFFSAR", label = "silac_13c")
  peptideMass("EGVNDNEEGFFSAR", aa = "R", shift = 6.020129)

  # To define multiple amino acids with the same modification, set \code{aa} as a single string.
  peptideMass("EGVNDNEEGFFSAR", aa = "RK", shift = 6.020129)

  # To define multiple amino acids with different modifications, set \code{aa} and \code{shift} as character vectors.
  peptideMass("EGVNDNEEGFFSAR", aa = c("RK", "C"), shift = c(6.020129, 57.021464), carbamidomethyl = FALSE)
}
\keyword{ genetics }
\keyword{ protein }
\keyword{ peptide }
