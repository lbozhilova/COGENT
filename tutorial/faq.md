# Frequently asked questions

This page aims to address some implementation and use case questions. We strongly recommend first reading the [tutorial](https://lbozhilova.github.io/COGENT/tutorial/tutorial.html) for an overview of COGENT's functionality and intended use.

## How do I install and load COGENT?

COGENT can be installed directly from R, using the `devtools` package. You can then load COGENT like any other R package.
```
# Install
require("devtools")
install_github("lbozhilova/COGENT")
# Load
require("COGENT")
```

For more detailed instructions please refer to the [tutorial](https://lbozhilova.github.io/COGENT/tutorial/tutorial.html).

## Which R version do I need?

COGENT is compatible with R v.3.6.2 and above. However, we recommend always using the most recent version of R and regularly updating your packages.

## Does COGENT work on all operating systems?

COGENT has been tested on Windows, MacOS and Linux (Ubuntu and Fedora). You can access all of its functionality on MacOS and Linux. However, Windows does not support forking, and so the paralellised functions of COGENT cannot be used to improve your computation speed on Windows. We therefore recommend that either MacOS or Linux is used for computationally intensive analysis.

## What split size (or equivalently what value of `propShared`) should I use when running COGENT?

A COGENT iteration works by randomly splitting your gene expression samples in two subsets, constructing a network from each subset, and comparing the resulting networks. These sample subsets can overlap, and the overlap is controlled by the `propShared` parameter (see `?cogentSingle` for details).

The default value `propShared=0` will always result in disjoint subsets. If your dataset consists of relatively few samples (say less than 50), we recommend introducing a set overlap. The best amount of overlap will depend on your data, but we recommend setting `propShared=0.25` or `propShared=0.5` in the first instance.
 
To illustrate, given twenty samples and `propShared=0`, `cogentSingle()` will generate two disjoint subsets of ten samples each. Ten samples may not be enough to construct an informative co-expression network. Setting `propShared=0.5` will result in two overlapping subsets, each consisting of fifteen samples. Of these, ten will be shared and five will be unique to each subset.

## How many COGENT iterations (`repCount`) should I run?

Multiple independent iterations of COGENT are required to provide the most robust results. The number of iterations is controlled by the `repCount` parameter (see `?cogentLinear` for details).

We generally recommend using COGENT with the default value of `repCount=100` iterations. Please note that `cogentParallel` can be used to reduce the computational time if needed.

## When should I use density correction (`getEdgeSimilarityCorrected`)?

Density correction, also known as density adjustment, is introduced in the  `getEdgeSimilarityCorrected()` function. It can be used to compare network construction methods which result in networks of significantly different density.

You should always use density correction if you are trying to find an optimal threshold for network construction, i.e. if you are trying to determine score cut-offs using COGENT. The [tutorial](https://lbozhilova.github.io/COGENT/tutorial/tutorial.html) contains a worked example of this.

You should also use density correction when comparing two network construction methods which are not guaranteed to result in networks of similar density.

You do not need to use density correction if you are using node metric consistency as your primary consistency measure. In this case, you should ensure your node metric of choice is density-independent.

You should not use density correction for methods resulting in signed and/or weighted networks.

## Which density correction should I use?

There are two density correction methods available through the `type` parameter of the `getEdgeSimilarityCorrected()` function. These are the semi-random correction (`type="expected"`) and the fully-random correction (`type="random"`). Refer to the [Supplementary Information](https://www.biorxiv.org/content/10.1101/2020.06.21.163535v1.supplementary-material) of our paper for details on how they are computed.

By default, the semi-random correction is carried out. This is significantly faster than the fully-random correction, and we recommend using it.

## I am interested in constructing signed, possibly weighted networks. Can I use COGENT?

COGENT is designed for unsigned networks. In particular, edge similarity calculated using `getEdgeSimilarity()` and `getEdgeSimilarityCorrected()` will be entirely non-informative when calculated for signed networks.

Nevertheless, there are several ways in which you can use COGENT for signed networks.

You can convert signed networks to unsigned by rescaling or using quantiles. The `buildPearson()` function in the [tutorial](https://lbozhilova.github.io/COGENT/tutorial/tutorial.html) illustrates one possible way of doing this - Pearson correlation coefficients can be negative, but if their quantiles are used instead, thresholding and weighted comparisons can be readily performed.

You can also use node metric consistency, either through the main COGENT functions or through `getNodeSimilarity()`, to evaluate signed network construction methods. Note that in this case, your node metric of choice should produce non-negative values. This can be achieved via rescaling, for example.

## This FAQ did not answer my question. How can I get in touch?

Please feel free to contact us with any questions by emailing Charlotte Deane at <deane@stats.ox.ac.uk>.

COGENT was originally developed and is currently maintained by Lyuba. You can email her directly at <lyubabozhilova@gmail.com>.