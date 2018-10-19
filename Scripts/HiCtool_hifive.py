"""
Program to run HiFive functions.

1) Creating a Fend object (hifive.Fend)

A Fragment-end (Fend) object (hdf5 format) contains information about the fragments created by 
digestion of a genome by a specific restriction enzyme (RE). In this program this information is supplied 
in the form of a BED-formatted file.

2) Creating a HiCData object (hifive.HiCData)

HiC dataset created from a Fend file and mapped data in BAM format. When a data object is created PCR duplicates
are removed and reads with ends mapping to adjacent fragments on opposite strands were also excluded, 
to consider the possibility of incomplete restriction enzyme digestion and fragment circularization.

3) Creating a HiC project object (hifive.HiC)

The HiC project object contains links to HiCData and Fend object, information about which fends to
include in the analysis, model parameters and learned model values. This is the standard way of
working with Hi-C data in HiFive and this object will be used for learning the model, extracting
portions of data, plotting and downstream analysis.

4) Filtering HiC fends (hifive.filter_fends)

Filtering out fends with specific properties to be not taken into account for learning fend correction
parameters.

5) Finding HiC distance function (hifive.find_distance_parameters)

Estimation of the distance-dependence relationship from the data prior to normalization, in order to
avoid biases that may result due to restriction site distribution characteristics or the influence 
of distance/signal relationship. This is related to the fact that restriction sites throughout the
genome are unevenly distributed and the interaction signal is strongly inversely-related to 
inter-fragment distance.

6) Learning correction parameters using the binning algorithm (hic.find_binning_fend_corrections)

Algorithm to learn the correction values for Hi-C data. We take into account of fragment
lengths, inter fragments distance, GC content and mappability biases for the normalization.

Note!
At the end of the each step the output file is saved (.hdf5 format) allowing the maximum flexibility 
for the user to update parameters and to run then the code from a specific step of analysis.

"""

import hifive

fend_file = 'restrictionsites_gc_map_valid.bed'
bam_file_1 = 'HiCfile_pair1.bam'
bam_file_2 = 'HiCfile_pair2.bam'
restriction_enzyme = 'MboI'

# Creating a Fend object
fend = hifive.Fend('fend_object.hdf5', mode='w')
fend.load_fends(fend_file, re_name=restriction_enzyme, format='bed')
fend.save()

# Creating a HiCData object
data = hifive.HiCData('HiC_data_object.hdf5', mode='w')
data.load_data_from_bam('fend_object.hdf5',
                        [bam_file_1,bam_file_2], 
                        maxinsert=500,
                        skip_duplicate_filtering=False)
data.save()

# Creating a HiC Project object
hic = hifive.HiC('HiC_project_object.hdf5', 'w')
hic.load_data('HiC_data_object.hdf5')
hic.save()

# Filtering HiC fends
hic = hifive.HiC('HiC_project_object.hdf5')
hic.filter_fends(mininteractions=1, mindistance=0, maxdistance=0)
hic.save()

# Finding HiC distance function
hic = hifive.HiC('HiC_project_object.hdf5')
hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)
hic.save()

# Learning correction parameters using the binning algorithm
hic = hifive.HiC('HiC_project_object.hdf5')
hic.find_binning_fend_corrections(max_iterations=1000,
                                  mindistance=500000,
                                  maxdistance=0,
                                  num_bins=[20,20,20,20],
                                  model=['len','distance','gc','mappability'],
                                  parameters=['even','even','even','even'],
                                  usereads='cis',
                                  learning_threshold=1.0)
hic.save('HiC_norm_binning.hdf5')