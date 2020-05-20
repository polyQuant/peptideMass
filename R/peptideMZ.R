peptideMZ <-
function(sequence, label="none", aa, shift, charge=2, carbamidomethyl=TRUE){

  # Check for correct input
  if(!is.numeric(charge) | length(charge) != 1){
    stop("Charge must be given as an integer (typically between 1-4).")
  }

  # Calculate mass of uncharged peptide
  mass <- proteinMass(sequence, label, monoisotopic = TRUE)

  # Add Carbamidomethylation with Iodoacetamide at cysteins
  if (carbamidomethyl = TRUE){
    nC <- stri_count(sequence, fixed = "C")
    mass <- mass + nC * 57.021464
  }

  # Add manual modifications
  massShift(sequence = sequence, aa = aa, shift = shift)

  # Modify for charged peptides.
  if (charge >= 0){
    mass <- mass + charge * 1.007276 # Add weight of H+ ions.
    mass <- mass / charge # Divide by charge state.
  }

  return(mass)
}
