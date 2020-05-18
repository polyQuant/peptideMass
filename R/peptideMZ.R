peptideMZ <-
function(sequence, charge=2, label="none"){

  # Check for correct input
  if(!is.numeric(charge) | length(charge) != 1){
    stop("Charge must be given as an integer (typically between 1-4).")
  }

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
