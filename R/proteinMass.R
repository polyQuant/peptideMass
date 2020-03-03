## Function to calculate mass of proteins/peptides
proteinMass <- function(sequence, monoisotopic = FALSE, label="none"){

  require(stringi)
  massTable <- read.csv("g:/programming/R_snippets/peptideMZ/AminoAcids_masses.csv", stringsAsFactors = F)
  label <- tolower(label) # Renders case-insensitive input string.

  # Calculate monoisotopic mass
  if(!(label %in% c("none", "silac_13c", "silac_15n13c", "15n"))){
    stop("Given label type unknown. Please use one of 'none', '15N', 'Silac_15N13C', or 'Silac_13C' (case-insensitive).")
  }
  aaMass <- massTable$Monoisotopic
  names(aaMass) <- massTable$oneLetters
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

  # Add mass shifts for SILAC labelling
  nK <- stri_count(sequence, fixed = "K")
  nR <- stri_count(sequence, fixed = "R")
  if(label == "silac_13c"){
    shiftK <- 6.020129
    shiftR <- 6.020129
    mass <- mass + nk * shiftK + nR * shiftR
  }else if(label == "silac_15n13c"){
    shiftK <- 8.014199
    shiftR <- 10.008269
    mass <- mass + nk * shiftK + nR * shiftR
  }

  # Add mass shifts for 15N labelling
  for (i in 1:nrow(massTable)){
    nN <- stri_count(sequence, fixed = massTable$oneLetter[i]) * massTable$AnzahlN
    mass <- mass + nN * 0.997035 # 0.997035 ist der mass shift zwischen 14N zu 15N. Da der natürliche Anteil an 15N nur 0.4% ist, kann der gleiche Unterschied für die Monoisotopic Mass und  für die Average Mass benutzt werden.
  }

  return(mass)
}

#ToDo:
# - korrekte Berechnung von 15N testen
# - evtl. noch Korrektur für avervage/monoisotopic 15N einführen
#  - evtl. noch eine Funktion aussenrum: mz aller Peptide eines PRoteins berechnen, mit Cleavage
#  - evtl. noch Hilfsfunktion (intern) schreiben, die ein bestimmtes Addukt an bestimmte AS hängt
#     - Berechnung von 15N über Mass table
#     - Addukte an einzelnen AA (C+57, KR+6/8/10) über die Hilfsfunktion
#  - better input check: allow only label "15N" and numeric charge states from 1-4.
#  - transfer into shiny app
