def pre_truncation(input_fastq,
                   restriction_enzyme, 
                   custom_ligation_junction='not_used'):
    """
    To perform pre-truncation on reads that contain potential ligation junctions. 
    To be executed before the mapping step.
    Arguments:
        input_fastq (str): path to the input fastq file.
        restriction_enzyme (str): restriction enzyme used in the Hi-C experiment. 
        One between 'HindIII', 'MboI', 'DpnII' and 'NcoI' (if different, you need to specify custom_ligation_junction)
        custom_ligation_junction (str): ligation junction sequence to be used if a different restriction enzyme rather than 
        the above ones is used.
    Output files:
        Fastq files with pre-truncated reads.
        Log files with pre-truncation information.
        Plot of the distribution of truncated reads length.
    Return: None.
    """
    
    from collections import Counter
    import os
    import re
    from math import ceil
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from time import gmtime, strftime
    
    if restriction_enzyme == 'HindIII':
        match = 'AAGCTAGCTT'
        re_seq = 'AAGCTT'
    elif restriction_enzyme == 'MboI' or restriction_enzyme == 'DpnII':
        match = 'GATCGATC'
        re_seq = 'GATC'
    elif restriction_enzyme == 'NcoI':
        match = 'CCATGCATGG'
        re_seq = 'CCATGG'
    else:
        match = custom_ligation_junction
    
    filename = os.path.splitext(input_fastq)[0]
    print "Pre-truncation of " + filename + ".fastq (restriction enzyme " + restriction_enzyme + ")..."
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
            if match in lines[i]:
                count += 1
                line=lines[i][0:-1] # remove the \n char at the end of the sequence
                pieces = line.split(match)
                max_length = max(len(x) for x in pieces)
                lengths.append(max_length + len(re_seq))
                piece = [x for x in pieces if len(x) == max_length][0]
                start_index = re.search(piece,line).start()
                match_index = re.search(match,line).start()
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
            fout.write(input_fastq + ', ' + restriction_enzyme + '\n' + result + '\n\n')
            
        rcParams.update({'figure.autolayout': True})
        cnt = Counter(lengths)
        my_lengths = [k for k, v in cnt.iteritems() if v >= 1]
        plt.close("all")
        plt.hist(lengths, bins=len(my_lengths))
        plt.title('Truncated reads (' + input_fastq + ')', fontsize=14)
        plt.xlabel('Read length', fontsize=14)
        plt.ylabel('Number of reads', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(input_fastq + '_truncated_reads.pdf', format = 'pdf')
        
        print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
