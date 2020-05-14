# Function to calculate additional mass by modifications on peptides
## To be called by functions 'proteinMass' and 'peptideMZ'

massShift <- function(sequence, label = "none"){

  # Check inputs
  label <- tolower(label) # Renders case-insensitive input string.
  if(!(label %in% c("none", "silac_13c", "silac_15n13c", "15n"))){
    stop("Given label type unknown. Please use one of 'none', '15N', 'Silac_15N13C', or 'Silac_13C' (case-insensitive).")
  }

  # Read mass table & dependencies
  massTable <- read.csv("g:/programming/R_snippets/peptideMZ/AminoAcids_masses.csv", stringsAsFactors = F)
  require(stringi)

  shiftTotal <- 0

  # Add mass shifts for SILAC labelling
  nK <- stri_count(sequence, fixed = "K")
  nR <- stri_count(sequence, fixed = "R")
  if(label == "silac_13c"){
    shiftK <- 6.020129
    shiftR <- 6.020129
    shiftTotal <- nK * shiftK + nR * shiftR
  }else if(label == "silac_15n13c"){
    shiftK <- 8.014199
    shiftR <- 10.008269
    shiftTotal <- nK * shiftK + nR * shiftR
  }

  # Add mass shifts for 15N labelling
  if(label == "15n"){
    for (i in 1:nrow(massTable)){
      nN <- stri_count(sequence, fixed = massTable$oneLetter[i]) * massTable$AnzahlN[i]
      shiftTotal <- shiftTotal + nN * 0.997035 # 0.997035 ist der mass shift zwischen 14N zu 15N. Da der natürliche Anteil an 15N nur 0.4% ist, kann der gleiche Unterschied für die Monoisotopic Mass und  für die Average Mass benutzt werden.
    }
  }

  return(shiftTotal)
}

## ToDo
# - Alternative input: modified aa ("aa") and mass shift ("shift")
#   - bei aa müssen mehrere angeben werden können, und auch N-terminal und C-terminal
#   - dann können auch die mass shifts für SILAC einfach über diese Funktion berechnet werden
