# Generate the fastq file used to align the restriction enzyme(s) to generate the fend file.

# Usage: python2.7 HiCtool_build_enzyme_fastq.py [-h] -e RESTRICTION_ENZYMES
# Options:
#  -h, --help               show this help message and exit
#  -e RESTRICTION_ENZYMES   Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses a combination of MboI and Hinfl.

# Output files:
#  Fastq file of the restriction enzyme(s).

from optparse import OptionParser

parameters = {'restriction_enzymes': None}

class enzyme_fastq:
    def __init__(self, parameters):
        self.generate_re_fastq(parameters)

    def generate_re_fastq(self, parameters):
    
        re_info = dict() # rs: restriction site sequence; lj: ligation junction sequence
        re_info['HindIII']={'rs':'AAGCTT'}
        re_info['MboI']={'rs':'GATC'}
        re_info['DpnII']={'rs':'GATC'}
        re_info['Sau3AI']={'rs':'GATC'}
        re_info['BglII']={'rs':'AGATCT'}
        re_info['NcoI']={'rs':'CCATGG'}
        re_info['Hinfl']={'rs':'GANTC'}
        
        restriction_enzymes = map(str, parameters['restriction_enzymes'].strip('[]').split(','))
        
        # Check that all the restriction enzymes are available
        for i in restriction_enzymes:
            if i not in re_info.keys():
                print "ERROR! " + i + " is not among the available restriction enzymes (HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl)! Check the spelling or contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>."
                return
        
        output_file = open ('restriction_enzyme.fastq', 'w')
        for restriction_enzyme in restriction_enzymes:
            if restriction_enzyme != 'Hinfl':
                output_file.write("@" + restriction_enzyme + '\n' + re_info[restriction_enzyme]['rs'] + '\n+\n' + ''.join(['I']*len(re_info[restriction_enzyme]['rs'])) + '\n')
            else:
                for i in ['A','C','G','T']:
                    output_file.write("@" + restriction_enzyme + '\n' + re_info[restriction_enzyme]['rs'].replace('N',i) + '\n+\n' + ''.join(['I']*len(re_info[restriction_enzyme]['rs'])) + '\n')
        output_file.close()
  
                  
if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_build_enzyme_fastq.py [-h] -e RESTRICTION_ENZYMES'
    parser = OptionParser(usage = 'python2.7 %prog -e restriction_enzymes')
    parser.add_option('-e', dest='restriction_enzymes', type='string', help='Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses both MboI and Hinfl.')
    (options, args) = parser.parse_args( )

    if options.restriction_enzymes == None:
        parser.error('-h for help or provide the restriction enzyme(s)!')
    else:
        pass

    parameters['restriction_enzymes'] = options.restriction_enzymes

    enzyme_fastq(parameters)