# Splits genome fasta files into individual files for every replicon,
# where every replicon has its genomic_accession as the file name, and the
# only data in the file is the dna sequence of the replicon. This is used
# as as part of an orthogonal test to ensure that the main copperon program
# is properly reading in genome_data from fasta files.

import os

fna_1 = "/home/katanga/Coding/copperon/test_assets/GCA_000152665.1_ASM15266v1/GCA_000152665.1_ASM15266v1_genomic.fna"
fna_2 = "/home/katanga/Coding/copperon/test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna"
fna_3 = "/home/katanga/Coding/copperon/test_assets/GCF_014107515.1_ASM1410751v1/GCF_014107515.1_ASM1410751v1_genomic.fna"
fna_4 = "/home/katanga/Coding/copperon/test_assets/GCF_016889785.1_ASM1688978v1/GCF_016889785.1_ASM1688978v1_genomic.fna"
output_dir = '/home/katanga/Coding/copperon/test_assets/preprocessed_test_genomes'

file_list = [fna_1, fna_2, fna_3, fna_4]

for fna in file_list:
    genomic_accessions = []
    definition_indicies = []
    replicons = []


    with open(fna, 'r') as genome_data:

        # Output logistics
        basename = os.path.basename(fna).rstrip('_genomic.fna')
        specific_output_dir = "{}/{}".format(output_dir, basename)
        os.mkdir(specific_output_dir)

        # Import and sanitize line data from file
        data = [line.strip('\n') for line in genome_data.readlines()]

        ## Define bounds for each entry in the genome
        for index, line in enumerate(data):
            if '>' in line:
                new_definition = line.split(" ")[0].strip('>')
                genomic_accessions.append(new_definition)
                definition_indicies.append(index)

        definition_indicies.append(len(data))
        sequence_bounds = []

        for start, end in zip(definition_indicies, definition_indicies[1:]):
            new_tup = (start+1, end)
            sequence_bounds.append(new_tup)

        ## Consolidate replicon sequences based on sequence bounds
        total_length = 0
        for name, (start, end) in zip(genomic_accessions, sequence_bounds):
            new_replicon = "".join(data[start:end])
            output_file = "{}/{}".format(specific_output_dir, name)

            with open(output_file, 'w') as output:
                output.writelines(new_replicon)