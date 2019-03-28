def generate_artificial_reads(genome_file,output_reads_file):
    """
        Split the input genome into 50 bp artificial reads every 10 bp and save the output in fastq format.
        Arguments:
            genome_file (str): genome fasta file.
            output_reads_file (str): output fastq file to save reads.
        Returns: None.
        Output:
            Fastq file containing all the artificial reads. A fixed and intermediate score "I" is given to each base of each read.
        """
    from Bio import SeqIO
    from time import gmtime, strftime
    
    print "Generating artificial reads..."
    print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
    
    with open(genome_file, "rU") as handle:
        records = list(SeqIO.parse(handle,"fasta"))
    
    for i in range(len(records)):
        sequence = records[i].seq
        reads = [str(sequence[x:x+50]) for x in range(0,len(sequence),10)]
        
        with open(output_reads_file, "a") as fout:
            for j in range(len(reads)):
                if len(reads[j]) == 50:
                    read = "@" + str(records[i].id) + "." + str(j) + "\n" + reads[j] + "\n" + "+" + "\n" + "I"*50 + "\n"
                    fout.write(read)
    
    print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())

