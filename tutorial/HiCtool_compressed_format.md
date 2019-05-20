# HiCtool compressed format

HiCtool contact map storage is documented in the Supplementary Material of [Calandrelli et al. (2018). GITAR: An open source tool for analysis and visualization of Hi-C data. *Genomics, proteomics & bioinformatics.*](https://www.sciencedirect.com/science/article/pii/S1672022918304339#s0055).

## Intra-chromosomal compressed format

Data are compressed based on the fact that contact maps are symmetric (contacts between loci i and j are the same than those between loci j and i) and usually sparse, since most of the elements are zeros, and this property is stronger with the decrease in the bin size. Given these two properties, it is not needed to save mirrored data and moreover it would be useful to “compress” the zero data within the matrices. To accomplish this, first we selected only the upper-triangular part of the contact matrices (including the diagonal) and reshaped the data by rows to form a vector. After that, we replaced all the consecutive zeros in the vector with a “0” followed by the number of zeros that are repeated consecutively; all the non-zero elements are left as they are. Finally, the data are saved in a txt file.

![](/figures/HiCtool_compression.png)

The figure shows a simplified example of the compression workflow, where the intra-chromosomal contact matrix is represented by a 4 × 4 symmetric and sparse matrix. (1) The upper-triangular part of the matrix is selected (including the diagonal); (2) data are reshaped to form a vector; (3) all the consecutive zeros are replaced with a “0” followed by the number of zeros that are repeated consecutively; (4) data are saved into a txt file.

### Comparison between full and compressed data for the entire human genome (only intra-chromosomal maps)

| Data type | Bin size (kb) | Storage usage txt format (MB) | Storage usage zip format (MB) | Saving time (min:sec) | Loading time |
|-----------|---------------|-------------------------------|-------------------------------|-----------------------|--------------|
| Full      | 1000          | 6.1                           | 2.7                           | 00:01                 | 00:01        |
| Optimized | 1000          | 2.9 (48%)                     | 1.4 (52%)                     | 00:01 (100%)          | 00:01 (100%) |
| Full      | 100           | 341.4                         | 105.9                         | 01:07                 | 00:30        |
| Optimized | 100           | 110.4 (32%)                   | 49.9 (47%)                    | 00:25 (37%)           | 00:17 (57%)  |
| Full      | 40            | 1566.7                        | 213                           | 06:18                 | 03:24        |
| Optimized | 40            | 208.6 (13%)                   | 93.1 (44%)                    | 01:50 (29%)           | 01:14 (36%)  |
| Full      | 10            | 19,957.80                     | 451.3                         | 90:16              | 85:41     |
| Optimized | 10            | 393.6 (2%)                    | 177.1 (39%)                   | 26:31 (29%)           | 18:49 (22%)  |

Hardware: 2.9 GHz Intel Core i5, 16 GB of RAM. The percentage of optimization (optimized/full) at each resolution is indicated in the parentheses. kb, kilobase; MB, megabyte.

## Inter-chromosomal compressed format

The same compression workflow applies for the inter-chromosomal contact matrix, **except for the selection of the upper triangular matrix (step 1)** which is skipped: the entire matrix data are considered (since the matrix is rectangular).


