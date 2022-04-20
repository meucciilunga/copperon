# TITLE: Database Analysis 2 - Annotated Assembly Count and Cross Database Comparison
# AUTHOR: Katavga
# DATE: 4-Apr-2022

# PURPOSE: Of the assembly directories that were downloaded from both the RefSeq and GenBank databases,
# this script summarizes the number of available annotation statistics for each, while comparing the
# uniqueness of specific assemblies present in each of the downloaded databases. The primary summary
# released from this script is two-fold: (1) the total number of unique AND downloaded assemblies across
# both databases, PLUS (2) the number of unique, downloaded, AND annotated assemblies available across
# both databases. Primary inputs are predefined by the script: the path to root_dir of (local) GenBank
# database and the path to root_dir of (local) RefSeq database; user must define output directory of
# summary file.

import os, sys

# Count the number of assemblies with annotation data available
def annotation_set(genomes_root_dir):

    # Pull out list of all subdirectories
    subdirs = [ f.path for f in os.scandir(genomes_root_dir) if f.is_dir() ]

    # Find all assemblies (which subdirs are a proxy for) with annotation data available
    annotated_assemblies = []
    for (index, genome_dir) in enumerate(subdirs): 
        
        # Progress bar
        msg = "{}/{}".format(index+1, len(subdirs))
        print(msg)  

        genome_name = os.path.basename(genome_dir)

        # Assembly is considered annotated if its directory contains these two files
        needed_files = [
            '_genomic.gff.gz',
            '_genomic.fna.gz',
        ]

        # Verify that relevant target files exist within each assembly subdirectory
        count = 0
        for file_type in needed_files:
            target = "{}/{}{}".format(genome_dir, genome_name, file_type)

            if os.path.isfile(target):
                count = count + 1
        
        if count == len(needed_files):
            annotated_assemblies.append(genome_dir)

    # full assembly names, including parent database prefix-tag (i.e., GCF_ for RefSeq, GCA_ for GenBank)
    full_assembly_names = [os.path.basename(dir_path) for dir_path in annotated_assemblies]

    # assembly name only; EXCLUDES parent database prefix-tag (i.e., GCF_ for RefSeq, GCA_ for GenBank)
    annotated_names = [name[4:] for name in full_assembly_names]

    # Returned as set for comparison operations needed later
    return(set(annotated_names))

# INPUT FILES
genbank_dir = '/run/media/katanga/SSD1/Prokaryotic Genomes Database/Genomes/GenBank'
refseq_dir = '/run/media/katanga/SSD1/Prokaryotic Genomes Database/Genomes/RefSeq'

# Specific subsets of GenBank database based on annotation status
genbank_list = os.listdir(genbank_dir)
genbank_assemblies = [item[4:] for item in genbank_list] # removes GCA_ tag from ASM file directories, which specifies its a GenBank entry
genbank_set = set(genbank_assemblies)                    # allows use of set comparison operators
genbank_annotated = annotation_set(genbank_dir)
genbank_unannotated = genbank_set.difference(genbank_annotated)

# Specific subsets of RefSeq database based on annotation status
refseq_list = os.listdir(refseq_dir)
refseq_assemblies = [item[4:] for item in refseq_list] # removes GCF_ tag from ASM file directories, which specifies its a RefSeq entry
refseq_set = set(refseq_assemblies)                    # allows use of set comparison operators
refseq_annotated = annotation_set(refseq_dir)
refseq_unannotated = refseq_set.difference(refseq_annotated)

# Total differences b/w the two databases among ALL assemblies
genbank_minus_refseq = genbank_set.difference(refseq_set)
refseq_minus_genbank = refseq_set.difference(genbank_set)

# Total differences b/w the two databases among ONLY ANNOTATED assemblies
gb_minus_rs_annotated = genbank_annotated.difference(refseq_annotated)
rs_minus_gb_annotated = refseq_annotated.difference(genbank_annotated)

# Total differences b/w the two databases among UNANNOTATED assemblies
gb_minus_rs_unannotated = genbank_unannotated.difference(refseq_set)
rs_minus_gb_unannotated = refseq_unannotated.difference(genbank_set)

# SUMMARIZE GenBank COMPARISON
print()
genbank_summary = [
    'Total Number of Complete Chromosome/Genome GenBank Entries:\n{}\n'.format(len(genbank_set)),
    'Total Number of Complete Chromosome/Genome GenBank Entries NOT in RefSeq:\n{}\n'.format(len(genbank_minus_refseq)),
    'Total Number of Complete Chromosome/Genome GenBank Entries WITH Annotation Data:\n{}\n'.format(len(genbank_annotated)),
    'Total Number of Complete Chromosome/Genome GenBank Entries WITHOUT Annotation Data:\n{}\n'.format(len(genbank_unannotated)),
    'Total Number of Complete Chromosome/Genome GenBank Entries WITH Annotation Data and NOT in RefSeq:\n{}\n'.format(len(gb_minus_rs_annotated)),
    'Total Number of Complete Chromosome/Genome GenBank Entries WITHOUT Annotation Data and NOT in RefSeq:\n{}\n'.format(len(gb_minus_rs_unannotated)),
]
for msg in genbank_summary:
    print(msg)

print('\n')

# SUMMARIZE RefSeq COMPARISON
refseq_summary = [
    'Total Number of Complete Chromosome/Genome RefSeq Entries:\n{}\n'.format(len(refseq_set)),
    'Total Number of Complete Chromosome/Genome RefSeq Entries NOT in GenBank:\n{}\n'.format(len(refseq_minus_genbank)),
    'Total Number of Complete Chromosome/Genome RefSeq Entries WITH Annotation Data:\n{}\n'.format(len(refseq_annotated)),
    'Total Number of Complete Chromosome/Genome RefSeq Entries WITHOUT Annotation Data:\n{}\n'.format(len(refseq_unannotated)),
    'Total Number of Complete Chromosome/Genome RefSeq Entries WITH Annotation Data and NOT in GenBank:\n{}\n'.format(len(rs_minus_gb_annotated)),
    'Total Number of Complete Chromosome/Genome RefSeq Entries WITHOUT Annotation Data and NOT in GenBank:\n{}\n'.format(len(rs_minus_gb_unannotated)),
]
for msg in refseq_summary:
    print(msg)
print('\n\n')

# BUILD SETS OF ASSEMBLIES BASED ON FILTERING FOR OUTPUT
refseq_annotated_prefixed = ['GCF_' + asm for asm in refseq_annotated]
gb_only_annotated_prefixed = ['GCA_' + asm for asm in gb_minus_rs_annotated]

unique_unannotated_genomes_prefixed = ['GCA_' + asm for asm in gb_minus_rs_unannotated]
all_annotated_genomes_prefixed = refseq_annotated_prefixed + gb_only_annotated_prefixed
all_genomes_prefixed = all_annotated_genomes_prefixed + unique_unannotated_genomes_prefixed

msg = "ACTUAL: Total number of AVAILABLE and UNIQUE genomes: {}".format(len(all_genomes_prefixed))
print(msg)

msg = "ACTUAL: Total number of AVAILABLE, UNIQUE, and ANNOTATED genomes: {}".format(len(all_annotated_genomes_prefixed))
print(msg)

msg = "ACTUAL: Total number of AVAILABLE, UNIQUE, and UNANNOTATED genomes: {}".format(len(unique_unannotated_genomes_prefixed))
print(msg)

# Output filtered list of assemblies to files
output_dir = sys.argv[1]
annotated_path_list_file = "{}/list_of_all_uniq_annotated.txt".format(output_dir)
unannotated_path_list_file = "{}/list_of_all_uniq_NON-annotated.txt".format(output_dir)

with open(annotated_path_list_file, 'w') as annotation_list_file:
    lines = '\n'.join(all_annotated_genomes_prefixed)
    annotation_list_file.write(lines)

with open(unannotated_path_list_file, 'w') as annotation_list_file:
    lines = '\n'.join(unique_unannotated_genomes_prefixed)
    annotation_list_file.write(lines)