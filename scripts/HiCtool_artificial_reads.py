# Split the input genome into 50 bp artificial reads every 10 bp and save the output in fastq format.

# Usage: python2.7 HiCtool_artificial_reads.py [-h] [options]
# Options:
#  -h, --help               show this help message and exit
#  -g GENOME_FILE           Genome fasta file.
#  -o OUTPUT_READS_FILE     Output fastq file to save reads.

# Output file:
#  Fastq file containing all the artificial reads. A fixed and intermediate score "I" is given to each base of each read.

from optparse import OptionParser
from Bio import SeqIO
from time import gmtime, strftime

parameters = {'genome_file': None,
              'output_reads_file': None}

class artificial_reads:
    def __init__(self, parameters):
        self.generate_artificial_reads(parameters)

    def generate_artificial_reads(self, parameters):
        with open(parameters['genome_file'], "rU") as handle:
            records = list(SeqIO.parse(handle,"fasta"))
        
        for i in range(len(records)):
            sequence = records[i].seq
            reads = [str(sequence[x:x+50]) for x in range(0,len(sequence),10)]
            
            with open(parameters['output_reads_file'], "a") as fout:
                for j in range(len(reads)):
                    if len(reads[j]) == 50:
                        read = "@" + str(records[i].id) + "." + str(j) + "\n" + reads[j] + "\n" + "+" + "\n" + "I"*50 + "\n"
                        fout.write(read)


if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_artificial_reads.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -g genome_file -o output_reads_file')
    parser.add_option('-g', dest='genome_file', type='string', help='Genome fasta file.')
    parser.add_option('-o', dest='output_reads_file', type='string', help='Output fastq file to save reads.')
    (options, args) = parser.parse_args( )

    if options.genome_file == None:
        parser.error('-h for help or provide the genome fasta file!')
    else:
        pass
    if options.output_reads_file == None:
        parser.error('-h for help or provide the output fastq file!')
    else:
        pass

    parameters['genome_file'] = options.genome_file
    parameters['output_reads_file'] = options.output_reads_file
    
    print "Generating artificial reads..."
    print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
    artificial_reads(parameters)
    print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
