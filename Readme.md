# PeptideMass package

This is an R package to calculate the correct MW and the correct m/z of proteins and peptides.

It is specifically designed to handle stable isotope labels. 
Currently implemented labels are:

- SILAC (13C or 13C15N at KR)
- 15N (all 14N in all amino acids replaced by 15N)

Other modifications may be included manually, please refer to the package help files in R.

The package can be installed with the following command to be executed in R (devtools needs to be installed):

```
devtools::install_github("https://github.com/polyQuant/peptideMass")
```