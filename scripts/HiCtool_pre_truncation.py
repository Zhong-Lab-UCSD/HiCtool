# Perform pre-truncation on reads that contain potential ligation junctions. To be executed before the mapping step.

# Usage: python2.7 HiCtool_pre_truncation.py [-h] [options]
# Options:
#  -h, --help               show this help message and exit
#  -i INPUTFILES            Input fastq file or input fastq files passed like a Python list. Example: [file1.fastq,file2.fastq].
#  -e RESTRICTION_ENZYMES   Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses a combination of MboI and Hinfl.
#  -p THREADS               Number of parallel threads to use.

# Output files:
#  Fastq files with pre-truncated reads.
#  Log files with pre-truncation information.

from optparse import OptionParser
import os
import re
from time import gmtime, strftime
from multiprocessing import Pool

parameters = {'inputFiles': None,
              'restriction_enzymes': None,
              'threads': None}

def pre_truncation(input_fastq):

    re_info = dict() # rs: restriction site sequence; lj: ligation junction sequence
    re_info['HindIII']={'rs':'AAGCTT','lj':'AAGCTAGCTT'}
    re_info['MboI']={'rs':'GATC','lj':'GATCGATC'}
    re_info['DpnII']={'rs':'GATC','lj':'GATCGATC'}
    re_info['Sau3AI']={'rs':'GATC','lj':'GATCGATC'}
    re_info['BglII']={'rs':'AGATCT','lj':'AGATCGATCT'}
    re_info['NcoI']={'rs':'CCATGG','lj':'CCATGCATGG'}
    re_info['Hinfl']={'rs':'GANTC','lj':'GA[ACGT]TA[ACGT]TC'}

    filename = os.path.splitext(input_fastq)[0]
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
                    print (key + ' completed - ' + input_fastq)
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
        print ('100% completed - ' + input_fastq)
        
        with open (filename + '_log.txt', 'w') as outlog:
            outlog.write(str(len(lines)/4) + "\t" + str(len(lines[1])-1) + "\t" + str(count) + "\n")


if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_pre_truncation.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -i inputFiles -e restriction_enzymes -p threads')
    parser.add_option('-i', dest='inputFiles', type='string', help='Input fastq file or input fastq files passed like a Python list. Example: [file1.fastq,file2.fastq].')
    parser.add_option('-e', dest='restriction_enzymes', type='string', help='Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses both MboI and Hinfl.')
    parser.add_option('-p', dest='threads', type='int', help='Number of parallel threads to use.')
    (options, args) = parser.parse_args()

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
    parameters['threads'] = options.threads
    
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
    threads = parameters['threads']
    
    # Check that all the restriction enzymes are available
    for i in restriction_enzymes:
        if i not in re_info.keys():
            parser.error(i + " is not among the available restriction enzymes (HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl)! Check the spelling or contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>.")
      
    if threads > 1:
        if threads > len(inputFiles):
            threads = len(inputFiles)
        if len(restriction_enzymes) == 1:
            print ("Pre-truncation in parallel (restriction enzyme: " + restriction_enzymes[0] + ") using " + str(threads) + " threads ...")
        else:
            print ("Pre-truncation in parallel (restriction enzymes: " + ", ".join(restriction_enzymes) + ") using " + str(threads) + " threads ...")
        print ("Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        pool = Pool(processes=threads)             
        pool.map(pre_truncation, inputFiles)
        print ("End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    else:
        if len(restriction_enzymes) == 1:
            print ("Pre-truncation (restriction enzyme: " + restriction_enzymes[0] + ") using a single threads ...")
        else:
            print ("Pre-truncation (restriction enzymes: " + ", ".join(restriction_enzymes) + ") using a single thread ...") 
        print ("Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        for i in inputFiles:
            pre_truncation(i)
        print ("End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime()))
