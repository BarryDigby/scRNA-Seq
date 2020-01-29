# Kallisto | Bustools tutorial

This tutorial is designed to guide the user through a `kallisto | bustools` pipeline to generate the 3 files required for downstream analysis: `matrix.mtx`, `genes.tsv` and `barcodes.txt`. As of writing this tutorial (28th January 2020) the Pachter lab is pushing it's new tool `kb-python` which is essentially a wrapper script that executes `kallisto | bustools` for the user with some python scripts under the hood for parsing output into the correct format. That being said, the [bustools github page](https://github.com/BUStools/BUS_notebooks_python/tree/master/dataset-notebooks) is home to numerous `ipynb` files for conducting scRNA-Seq analysis using `kallisto | bustools` which are inconsistent and result in incomplete output files. Perhaps more aggregious is the [getting started](https://www.kallistobus.tools/getting_started) bustools tutorial (directed from the [bustools website](https://bustools.github.io/) which leaves the user with unannotated gene ID outputs and no counsel on how to rectify the problem. 

First things first, an outline of what the problem is. I will use the "getting started" tutorial mentioned above to illustrate it.



