# HiCtool compressed format

HiCtool contact map storage is documented in the Supplementary Material of [Calandrelli et al. (2018). GITAR: An open source tool for analysis and visualization of Hi-C data. *Genomics, proteomics & bioinformatics.*](https://www.sciencedirect.com/science/article/pii/S1672022918304339#s0055).

## Intra-chromosomal compressed format

Data are compressed based on the fact that contact maps are symmetric (contacts between loci i and j are the same than those between loci j and i) and usually sparse, since most of the elements are zeros, and this property is stronger with the decrease in the bin size. Given these two properties, it is not needed to save mirrored data and moreover it would be useful to “compress” the zero data within the matrices. To accomplish this, first we selected only the upper-triangular part of the contact matrices (including the diagonal) and reshaped the data by rows to form a vector. After that, we replaced all the consecutive zeros in the vector with a “0” followed by the number of zeros that are repeated consecutively; all the non-zero elements are left as they are. Finally, the data are saved in a txt file.

![](/figures/HiCtool_compression.png)

The figure shows a simplified example of the compression workflow, where the intra-chromosomal contact matrix is represented by a 4 × 4 symmetric and sparse matrix. (1) The upper-triangular part of the matrix is selected (including the diagonal); (2) data are reshaped to form a vector; (3) all the consecutive zeros are replaced with a “0” followed by the number of zeros that are repeated consecutively; (4) data are saved into a txt file.

## Inter-chromosomal compressed format

The same compression workflow applies for the inter-chromosomal contact matrix, **except for the selection of the upper triangular matrix (step 1)** which is skipped: the entire matrix data are considered (since the matrix is rectangular).



