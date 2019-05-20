# Perform pre-truncation on reads that contain potential ligation junctions. To be executed before the mapping step.

# Usage: python2.7 HiCtool_pre_truncation.py [-h] [options]
# Options:
#  -h, --help               show this help message and exit
#  -i INPUTFILES            Input fastq file or input fastq files passed like a Python list. Example: [file1.fastq,file2.fastq].
#  -e RESTRICTION_ENZYMES   Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses a combination of MboI and Hinfl.

# Output files:
#  Fastq files with pre-truncated reads.
#  Log files with pre-truncation information.
#  Plot of the distribution of truncated reads length.

from optparse import OptionParser
from collections import Counter
import os
import re
from math import ceil
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import gmtime, strftime

parameters = {'inputFiles': None,
              'restriction_enzymes': None}

class pre_truncate:
    def __init__(self, parameters):
        self.pre_truncation(parameters)

    def pre_truncation(self, parameters):
    
        re_info = dict() # rs: restriction site sequence; lj: ligation junction sequence
        re_info['HindIII']={'rs':'AAGCTT','lj':'AAGCTAGCTT'}
        re_info['MboI']={'rs':'GATC','lj':'GATCGATC'}
        re_info['DpnII']={'rs':'GATC','lj':'GATCGATC'}
        re_info['Sau3AI']={'rs':'GATC','lj':'GATCGATC'}
        re_info['BglII']={'rs':'AGATCT','lj':'AGATCGATCT'}
        re_info['NcoI']={'rs':'CCATGG','lj':'CCATGCATGG'}
        re_info['Hinfl']={'rs':'GANTC','lj':'GA[ACGT]TA[ACGT]TC'}
        
        inputFiles = map(str, parameters['inputFiles'].strip('[]').split(','))
        restriction_enzymes = map(str, parameters['restriction_enzymes'].strip('[]').split(','))
        
        # Check that all the restriction enzymes are available
        for i in restriction_enzymes:
            if i not in re_info.keys():
                print "ERROR! " + i + " is not among the available restriction enzymes (HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl)! Check the spelling or contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>."
                return
        
        for input_fastq in inputFiles:
            filename = os.path.splitext(input_fastq)[0]
            if len(restriction_enzymes) == 1:
                print "Pre-truncation of " + filename + ".fastq (restriction enzyme: " + restriction_enzymes[0] + ")..."
            else:
                print "Pre-truncation of " + filename + ".fastq (restriction enzymes: " + ", ".join(restriction_enzymes) + ")..."
            print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
            with open (input_fastq, 'r') as infile:
                lines = infile.readlines() 
                count = 0 # to count reads with potential ligation junction
                lengths = [] # to save the length of the pieces after truncation
                n_reads = len(lines)/4
                percents = {}
                for n in range(0,101,10)[1:]:
                    percents[str(n)+'%'] = int(n_reads*n*0.01)
                for i in xrange(1,len(lines),4): # iteration over lines containing the sequence
                    for key, value in percents.iteritems():
                        if i/4 == value:
                            print key + ' completed.'
                    # Search for the ligation junction
                    for j in restriction_enzymes:
                        match = re.search(re_info[j]['lj'],lines[i])
                        if match: # ligation junction present in the read sequence
                            count += 1
                            line=lines[i][0:-1] # remove the \n char at the end of the sequence
                            pieces = []
                            pieces.append(line[:match.start()])
                            pieces.append(line[match.end():])
                            max_length = max(len(x) for x in pieces)
                            if j != 'Hinfl':
                                re_seq = re_info[j]['rs']
                            else:
                                re_seq = match.group()[:4] + 'C'
                            lengths.append(max_length + len(re_seq))
                            piece = [x for x in pieces if len(x) == max_length][0]
                            start_index = re.search(piece,line).start()
                            match_index = match.start()
                            if start_index < match_index:
                                piece = piece + re_seq
                                lines[i] = piece + '\n'
                                lines[i+2] = lines[i+2][:-1][start_index:start_index+len(piece)] + '\n'
                            elif start_index > match_index:
                                piece = re_seq + piece
                                lines[i] = piece + '\n'
                                lines[i+2] = lines[i+2][:-1][start_index-len(re_seq):start_index-len(re_seq)+len(piece)] + '\n'
        
                with open (filename + '.trunc.fastq' ,'w') as fout:
                    for i in xrange(len(lines)):
                        fout.write(lines[i])
                print '100% completed.'
                
                perc_reads = ceil(float(count)/float(len(lines)/4)*10000)/100.0
                result = str(len(lines)/4) + ' reads (length = ' + str(len(lines[1])-1) + ' bp); of these:\n  ' + str(count) + ' (' + str(perc_reads) + '%) contained a potential ligation junction and have been truncated.'
                print result
        
                with open ('pre_truncation_log.txt', 'a') as fout:
                    if len(restriction_enzymes) == 1:
                        fout.write(input_fastq + ', ' + restriction_enzymes[0] + '\n' + result + '\n\n')
                    else:
                        fout.write(input_fastq + ', ' + ", ".join(restriction_enzymes) + '\n' + result + '\n\n')
                        
                rcParams.update({'figure.autolayout': True})
                cnt = Counter(lengths)
                my_lengths = [k for k, v in cnt.iteritems() if v >= 1]
                plt.close("all")
                plt.hist(lengths, bins=len(my_lengths))
                plt.title('Truncated reads (' + input_fastq + ')', fontsize=14)
                plt.xlabel('Read length', fontsize=12)
                plt.ylabel('Number of reads', fontsize=12)
                plt.xticks(fontsize=10)
                plt.yticks(fontsize=10)
                plt.savefig(input_fastq + '_truncated_reads.pdf', format = 'pdf')
                
                print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())


if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_pre_truncation.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -i inputFiles -e restriction_enzymes')
    parser.add_option('-i', dest='inputFiles', type='string', help='Input fastq file or input fastq files passed like a Python list. Example: [file1.fastq,file2.fastq].')
    parser.add_option('-e', dest='restriction_enzymes', type='string', help='Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses both MboI and Hinfl.')
    (options, args) = parser.parse_args( )

    if options.inputFiles == None:
        parser.error('-h for help or provide the input fastq file(s)!')
    else:
        pass
    if options.restriction_enzymes == None:
        parser.error('-h for help or provide the restriction enzyme(s)!')
    else:
        pass

    parameters['inputFiles'] = options.inputFiles
    parameters['restriction_enzymes'] = options.restriction_enzymes

    pre_truncate(parameters)