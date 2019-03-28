def run_hifive(fend_file,
               bam_file_1,
               bam_file_2,
               restriction_enzyme,
               model):
    """
    Program to run HiFive functions.
    Arguments:
        fend_file: fend file from preprocessing.
        bam_file_1: bam file for the first read of the pairs.
        bam_file_2: bam file for the second read of the pairs.
        restriction_enzyme (str): restriction enzyme used in the Hi-C experiment.
        model (str): model you wish to use for the following normalization procedure: either 'Yaffe-Tanay' or 'Hi-Corrector'.
    Return: None.
    """
    import hifive
    
    if model == 'Yaffe-Tanay':
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
    
    elif model == 'Hi-Corrector':
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
