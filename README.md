[![DOI](https://zenodo.org/badge/424666838.svg)](https://zenodo.org/badge/latestdoi/424666838)

# R scripts and Supplementary Material for the Manuscript "How cultural transmission through objects impacts inferences about cultural evolution"

This repository contains the Rmarkdown supplementary material and the R scripts for generating simulation outputs and figures for the manuscript:

Crema, E.R., Bortolini, E. & Lake, M.W. (2023). How cultural transmission through objects impacts inferences about cultural evolution. _Journal of Archaeological Method and Theory_.

## File Structure

* `src.R` ... contains utility functions and simulation code for all experiments presented in the paper.
* `runscript.R` ... executes all simulations.
* `simres.RData` ... R image file containing the simulation outputs (generated from `runscript.R`)
* `makefigures.R` ... generates figures.
* `figures/*.pdf` ... pdf version of all figures in the paper. Generated using scripts contained in `makefigures.R`.
* `esm.Rmd` ... Suplementary material (in Rmarkdown format) containing detailed description of the models and the experiment design.
* `esm.html` ... [Rendered version](https://htmlpreview.github.io/?https://github.com/ercrema/objectmediated_transmission/blob/main/esm.html) of `esm.Rmd`.

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

