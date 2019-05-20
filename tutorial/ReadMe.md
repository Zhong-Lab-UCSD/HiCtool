# HiCtool Tutorial (v2.1)

This is a full tutorial of HiCtool. The tutorial steps have to be followed in sequence. At the normalization step, choose the preferred method.

All the scripts used in the tutorial are inside the main folder [scripts](https://github.com/Zhong-Lab-UCSD/HiCtool/tree/master/scripts). In this tutorial, all the commands are run supposing that you have downloaded the ``HiCtool-master`` folder inside your working directory, otherwise update the paths when you run the commands accordingly.

![](/tutorial/HiCtool_flowchart.png)

## [1. Data preprocessing](/tutorial/data-preprocessing.md)
## 2. Data normalization and visualization
The **explicit-factor correction model of Yaffe and Tanay** is applied to normalize (and visualize) only intra-chromosomal contact data and takes into account of explicit factors such as fragment length, GC content, mappability score to normalize the data. The **matrix balancing approach of Hi-Corrector** is used to normalize and visualize globally intra- and inter-chromosomal contact maps. It does not consider specific factors for normalization and it is faster.
- ### [2.1. Explicit-factor correction model of Yaffe and Tanay](/tutorial/normalization-yaffe-tanay.md)
- ### [2.2. Matrix balancing approach of Hi-Corrector](/tutorial/normalization-matrix-balancing.md)
## [3. TAD analysis](/tutorial/tad-analysis.md)

***
## Supplementary information

- ### [HiCtool utility code](/tutorial/HiCtool_utility_code.md)
- ### [HiCtool contact matrix compressed format](/tutorial/HiCtool_compressed_format.md)

