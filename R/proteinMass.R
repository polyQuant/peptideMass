proteinMass <-
function(sequence, label="none", monoisotopic = FALSE){

  # Check for correct input
  if(!all(strsplit(sequence, "")[[1]] %in% massTable$oneLetter)){
    stop("The input sequence contains unknown amino acids. Did you paste the sequence without Fasta header?")
  }

  # Calculate monoisotopic mass
  aaMass <- massTable$Monoisotopic
  names(aaMass) <- massTable$oneLetter
  mass <- numeric(length = 1)
  for(i in 1:length(aaMass)){
    num <- stri_count(sequence, fixed = names(aaMass)[i])
    mass <- mass + num * as.numeric(aaMass[i])
  }
  mass <- mass + 18.01056 # Add weight of water

  # Calculate average mass
  if(!monoisotopic){
    massDiff <- massTable$diffAveMono
    names(massDiff) <- massTable$oneLetter
    for(i in 1:length(massDiff)){
      num <- stri_count(sequence, fixed = names(massDiff)[i])
      mass <- mass + num * as.numeric(massDiff[i])
    }
  }

  # Add mass shifts from labels
  mass <- mass + massShift(sequence = sequence, label = label)

  return(mass)
}
