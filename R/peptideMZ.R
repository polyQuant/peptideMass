## Function to calculate m/z of peptides
peptideMZ <- function(sequence, charge=2, label="none"){

  # Calculate mass of uncharged peptide
  mass <- proteinMass(sequence, label, monoisotopic = TRUE)

  # Add alkylations at cysteins
  nC <- stri_count(sequence, fixed = "C")
  mass <- mass + nC * 57.021464 # For Carbamidomethylation with Iodoacetamide

  # Modify for charged peptides.
  if (charge >= 0){
    mass <- mass + charge * 1.007276 # Add weight of H+ ions.
    mass <- mass / charge # Divide by charge state.
  }

  return(mass)
}
