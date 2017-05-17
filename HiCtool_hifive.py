"""
Program to run HiFive functions.

1) Creating a Fend object (hifive.Fend)

A Fragment-end (Fend) object (hdf5 format) contains information about the fragments created by 
digestion of a genome by a specific restriction enzyme (RE). In this program this information is supplied 
in the form of a BED-formatted file.

2) Creating a HiCData object (hifive.HiCData)

HiC dataset created from a Fend file and mapped data in BAM format. The original BAM file have to be
splitted into two BAM files each containing the information of one of the two mapped reads, and then
they are passed to the function as elements of a list.

3) Creating a HiC project object (hifive.HiC)

The HiC project object contains links to HiCData and Fend object, information about which fends to
include in the analysis, model parameters and learned model values. This is the standard way of
working with Hi-C data in HiFive and this object will be used for learning the model, extracting
portions of data, plotting and downstream analysis.

4) Filtering HiC fends (hifive.filter_fends)

Filtering out fends which specific properties to be not taken into account for learning fend correction
parameters.

5) Finding HiC distance function (hifive.find_distance_parameters)

Estimation of the distance-dependence relationship from the data prior to normalization, in order to
avoid biases that may result due to restriction site distribution characteristics or the influence 
of distance/signal relationship. This is related to the fact that restriction sites throughout the
genome are unevenly distributed and the interaction signal is strongly inversely-related to 
inter-fragment distance.

6) Learning correction parameters using the binning algorithm (hic.find_binning_fend_corrections)

Algorithm to learn the correction values for Hi-C data. We take into account of fragment
lengths and GC content biases for the normalization.

Note!
At the end of the each step the output file is saved (.hdf5 format) allowing the maximum flexibility 
for the user to update parameters and to run then the code from a specific step of analysis.

"""

import hifive

# Creating a Fend object
fend = hifive.Fend('fend_object.hdf5', mode='w')
fend.load_fends('HindIII_hg38_gc.bed', re_name='HindIII', format='bed')
fend.save()


# Creating a HiCData object
data = hifive.HiCData('HiC_data_object.hdf5', mode='w')
data.load_data_from_bam('fend_object.hdf5',
                        ['HiCfile_pair1.bam','HiCfile_pair2.bam'], 
                        maxinsert=500)
data.save()


# Creating a HiC Project object
hic = hifive.HiC('HiC_project_object.hdf5', 'w')
hic.load_data('HiC_data_object.hdf5')
hic.save()


# Filtering HiC fends

hic = hifive.HiC('HiC_project_object.hdf5')
hic.filter_fends(mininteractions=1, mindistance=500000, maxdistance=0)
hic.save()


# Finding HiC distance function
hic = hifive.HiC('HiC_project_object.hdf5')
hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)
hic.save('HiC_distance_function.hdf5')


# Learning correction parameters using the binning algorithm
hic = hifive.HiC('HiC_distance_function.hdf5')
hic.find_binning_fend_corrections(max_iterations=1000,
                                  mindistance=500000,
                                  maxdistance=0,
                                  num_bins=[20,20],
                                  model=['len','gc'],
                                  parameters=['even','even'],
                                  usereads='cis',
                                  learning_threshold=1.0)
hic.save('HiC_norm_binning.hdf5')