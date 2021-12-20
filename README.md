# R scripts and Supplementary Material for the Manuscript "Object-mediated cultural transmission and inference"

This repository contains the Rmarkdown supplementary material and the R scripts for generating simulation outputs and figures for the manuscript:

Crema, E.R., Bortolini, E. & Lake, M.W. (2021). How cultural transmission through objects impacts inferences about cultural evolution

## File Structure

* `src.R` ... contains utility functions and simulation code for all experiments presented in the paper.
* `runscript.R` ... executes all simulations and generates figures presented in the paper.
* `figures/*.pdf` ... pdf version of all figures in the paper. Generated using scripts contained in `runscript.R`.
* `esm.Rmd` ... Suplementary material (in Rmarkdown format) containing detailed description of the models and the experiment design.

## Required R packages

```
attached base packages:                                                         
[1] parallel  stats     graphics  grDevices utils     datasets  methods         
[8] base                                                                        
                                                                                
other attached packages:                                                        
[1] RColorBrewer_1.1-2 doParallel_1.0.16  iterators_1.0.13                      
[4] foreach_1.5.1      here_1.0.1         nvimcom_0.9-115                       
                                                                                
loaded via a namespace (and not attached):                                      
[1] compiler_4.1.1   rprojroot_2.0.2  tools_4.1.1      codetools_0.2-18 
```

## Licence
CC-BY 3.0

