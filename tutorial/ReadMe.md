# HiCtool Tutorial

The sections in this tutorial steps have to be followed in sequence. At the normalization step, choose the preferred method. Note that A/B compartment analysis can be performed only if the method from Yaffe-Tanay is used, which allows to extract the O/E contact matrix and the Pearson correlation matrix. The method using Hi-Corrector is usually preferred if you have highly deep sequenced data, and it allows both intra- and inter-chromosomal analysis.

All the scripts used in the tutorial are inside the main folder [scripts](https://github.com/Zhong-Lab-UCSD/HiCtool/tree/master/scripts). In this tutorial, all the commands are run supposing that you have downloaded the ``HiCtool-master`` folder inside your working directory, otherwise update the paths when you run the commands accordingly.

![](/tutorial/HiCtool_flowchart.png)

## [1. Data preprocessing](/tutorial/data-preprocessing.md)
## 2. Data normalization and visualization
The **explicit-factor correction model of Yaffe and Tanay** is applied to normalize and visualize only intra-chromosomal contact data and it takes into account of explicit factors such as fragment length, GC content, mappability score to normalize the data. It has to be used to perform A/B compartment analysis. 

The **matrix balancing approach of Hi-Corrector** is used to normalize and visualize globally intra- and inter-chromosomal contact maps. It does not consider specific factors for normalization and it is faster. Both methods are fine for TAD analysis.
- ### [2.1. Explicit-factor correction model of Yaffe and Tanay](/tutorial/normalization-yaffe-tanay.md)
- ### [2.2. Matrix balancing approach of Hi-Corrector](/tutorial/normalization-matrix-balancing.md)
## [3. A/B compartment analysis](/tutorial/compartment.md)
## [4. TAD analysis](/tutorial/tad-analysis.md)

***
## Supplementary information

- ### [HiCtool utility code](/tutorial/HiCtool_utility_code.md)
- ### [HiCtool contact matrix compressed format](/tutorial/HiCtool_compressed_format.md)

