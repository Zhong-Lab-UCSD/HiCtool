# Program to run HiFive functions.

# Usage: python2.7 HiCtool_hifive.py [-h] [options]
# Options:
#  -h, --help                   show this help message and exit
#  -f    FEND_FILE              Fend file from preprocessing
#  --b1  BAM_FILE_1             BAM file 1 for the first read of the pairs
#  --b2  BAM_FILE_1             BAM file 2 for the second read of the pairs
#  -e    RESTRICTION_ENZYME     Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses a combination of MboI and Hinfl.
#  -m    MODEL                  Model you wish to use for the following normalization procedure: either 'Yaffe-Tanay' or 'Hi-Corrector'.

# Output files:
#  HiC_project_object.hdf5 and HiC_norm_binning.hdf5 if model='Yaffe-Tanay'
#  HiC_project_object.hdf5 if model='Hi-Corrector'

from optparse import OptionParser
import hifive

parameters = {'fend_file': None,
              'bam_file_1': None,
              'bam_file_2': None,
              'restriction_enzyme': None,
              'model': None}

class hi_five:
    def __init__(self, parameters):
        self.run_hifive(parameters)
    
    def run_hifive(self, parameters):
        
        fend_file = parameters['fend_file']
        bam_file_1 = parameters['bam_file_1']
        bam_file_2 = parameters['bam_file_2']
        model = parameters['model']
        
        restriction_enzymes = map(str, parameters['restriction_enzyme'].strip('[]').split(','))
        if len(restriction_enzymes) == 1:
            restriction_enzyme = restriction_enzymes[0]
        else:
            restriction_enzyme = ','.join(restriction_enzymes)
    
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
            
            
            
if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_hifive.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -f fend_file --b1 bam_file_1 --b2 bam_file_2 -e restriction_enzyme -m model')
    parser.add_option('-f', dest='fend_file', type='string', help='Fend file from preprocessing.')
    parser.add_option('--b1', dest='bam_file_1', type='string', help='BAM file 1 for the first read of the pairs.')
    parser.add_option('--b2', dest='bam_file_2', type='string', help='BAM file 2 for the second read of the pairs.')
    parser.add_option('-e', dest='restriction_enzyme', type='string', help='Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses both MboI and Hinfl.')
    parser.add_option('-m', dest='model', type='string', help='Model you wish to use for the following normalization procedure: either "Yaffe-Tanay" or "Hi-Corrector".')
    (options, args) = parser.parse_args( )

    if options.fend_file == None:
        parser.error('-h for help or provide the fend file!')
    else:
        pass
    if options.bam_file_1 == None:
        parser.error('-h for help or provide the bam file 1!')
    else:
        pass
    if options.bam_file_2 == None:
        parser.error('-h for help or provide the bam file 2!')
    else:
        pass
    if options.restriction_enzyme == None:
        parser.error('-h for help or provide the restriction enzyme(s)!')
    else:
        pass
    if options.model == None:
        parser.error('-h for help or provide the model: "Yaffe-Tanay" or "Hi-Corrector"!')
    else:
        pass

    parameters['fend_file'] = options.fend_file
    parameters['bam_file_1'] = options.bam_file_1
    parameters['bam_file_2'] = options.bam_file_2
    parameters['restriction_enzyme'] = options.restriction_enzyme
    parameters['model'] = options.model

    if parameters['model'] in ['Yaffe-Tanay', 'Hi-Corrector']:
        print "Running HiFive functions for the model " + parameters['model'] + " ..."
    else:
        parser.error("Please insert the correct model, one between Yaffe-Tanay and Hi-Corrector.")
        
    hi_five(parameters)
