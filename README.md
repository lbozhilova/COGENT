# COGENT: COnsistency of Gene Expression NeTworks

## Overview
COGENT is an R package designed to aid the construction of gene co-expression networks without the need for annotation or other external validation data. COGENT can be used to choose between competing co-expression measures (e.g. Pearson vs. Kendall correlation coefficients), as well as to inform score cut-off choice. While designed for gene expression data, COGENT can be applied to other cases where network construction relies on similarity profiling, e.g. microbiome or synthetic lethality data.

## Installation
In order to make use of COGENT's parallel computing functionality, we recommend working in Linux or MacOS. Regardless of your operating system choice, the easiest way to install COGENT is from R, using  `devtools::install_github()`:
`install_github("lbozhilova/COGENT")`

See [the tutorial](https://lbozhilova.github.io/COGENT/tutorial/tutorial.html) for more detailed instructions on installing COGENT.

## Usage
The main COGENT functions - `cogentSingle()`, `cogentLinear()` and `cogentParallel()` - have two principal arguments. The first is a dataframe containing the gene expression data of interest, and the second is a function mapping the dataframe to a network adjacency  matrix. The output is a set of consistency measures, which quantify the robustness, or stability, of the network construction method (i.e. the function argument) as applied to the data of interest (i.e. the dataframe). The simplest COGENT function call takes the form
`cogentLinear(myData, myFunction)`.

For details on what consistency measures are available within COGENT, and how these are calculated, please refer to [the tutorial](https://lbozhilova.github.io/COGENT/tutorial/tutorial.html), which contains a fully worked example and commentary.

## Reference
If you use COGENT in your work, please cite us:

Bozhilova, L.V., Pardo-Diaz, J., Reinert, G. & Deane, C.M. [_COGENT: evaluating the consistency of gene co-expression networks_](https://www.biorxiv.org/content/10.1101/2020.06.21.163535v1) bioRxiv, 2020.

For an application of COGENT, refer to:

Pardo-Diaz, J., Bozhilova, L.V., Beguerisse-Diaz, M., Poole, P.S., Deane, C.M. & Reinert, G. [_Robust gene coexpression networks using signed distance correlation_](https://www.biorxiv.org/content/10.1101/2020.06.21.163543v1) bioRxiv, 2020.

## Updates
Nothing to see here... yet.
