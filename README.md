# TDMresults

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32851.svg)](http://dx.doi.org/10.5281/zenodo.32851)

Scripts and data for re-generating results from the TDM manuscript.

Use run_experiments.R to regenerate the analyses from
[Thompson et al.](https://peerj.com/preprints/1460/).

We are using Git Large File Storage to store gene expression data. You will
need to [install the client](https://git-lfs.github.com/) if you do not already
have it installed. Once you install the client, run `git lfs pull` in the
clone to retrieve the datasets. After retrieving the datasets, using
run_experiments.R will regenerate the results.

This requires the following R & Bioconductor packages be installed:

	* ggplot2
	* reshape2
	* Hmisc
	* data.table
	* scales
	* sdcMicro
	* flexclust
	* fpc
	* corrplot
	* ape
	* cluster
	* plyr
	* dplyr
	* devtools
	* quantro
	* preprocessCore
	* gridExtra
	* huge
	* caret
	* limma
	* glmnet
	* e1071
	* stringr
	* gdata
	* binr
   	* cowplot

One github package is required:

    library("devtools")
    install_github(repo = "quantroSim", username = "stephaniehicks")

By default, all input data is expected to be in this local directory and
normalized data and figure output is within this directory as well. If you'd
like to change this behavior, modify "directories.R".

Acknowledgements:
This research is funded in part by the Gordon and Betty Moore Foundation’s
Data-Driven Discovery Initiative through Grant GBMF4552 to CSG. JT is a Neukom
Graduate Fellow supported by the William H. Neukom 1964 Institute for
Computational Science. This work was supported in part by P20 GM103534,
P30 CA023108 and UL1 TR001086 from the NIH and an American Cancer Society
Research Grant, #IRG-82-003-27. The funders had no role in study design, data
collection and analysis, decision to publish, or preparation of the manuscript.
