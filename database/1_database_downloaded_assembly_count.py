# TITLE: Database Analysis 1 - Available Assembly Count
# AUTHOR: Katavga
# DATE: 4-Apr-2022

# PURPOSE: This script calculates the number of assemblies containing either a complete genome/chromosome
# that are available in a given remote database (i.e., RefSeq or GenBank) based on their metadata summary
# files; it then compares that calculated number of complete genomes/chromosomes with the number of assembly-
# specific directories currently available within a specified root directory.

import os
import sys

genomes_database_dir = sys.argv[1]
genomes_database_metadata = sys.argv[2]

subdir_list = os.listdir(genomes_database_dir)

with open(genomes_database_metadata, 'r') as metadata:
    valid_genomes = []
    valid_genome_names = []
    valid_links = []
    valid_genomes_by_name = dict()

    for index, item in enumerate(metadata.readlines()):

        # Skip first line of the file
        if index == 0:
            continue
        
        # Split lines by \t marker
        entry = item.split('\t')

        # Pull specific data for each entry
        assembly_accession = entry[0]
        biosample_accession = entry[2]
        ref_seq_cat = entry[4]
        tax_id = entry[5]
        species_tax_id = entry[6]
        organism_name = entry[7]
        infraspecifc_name = entry[8]
        completion_status = entry[11].lower()
        ftp_location = entry[19]
        subdir_name = ftp_location.split('/')[-1] 
        ref_seq_status = entry[20]

        # Collapse data into a new csv entry
        new_entry = [
            tax_id,
            species_tax_id,
            organism_name,
            infraspecifc_name,
            completion_status,
            ref_seq_cat,
            ref_seq_status,
            assembly_accession,
            biosample_accession,
            subdir_name,
            "rsync" + ftp_location.lstrip("https") 
            # all database links start w/ 'https' prefix, so we remove that and replace w/ 'rsync'
        ]

        # Filter genomes by assembly level
        complete_status = (completion_status == 'complete genome')
        chromosome_status = (completion_status == 'chromosome')
        valid_ftp_link = (ftp_location != "na")

        valid_genome_cond = (complete_status or chromosome_status) and valid_ftp_link
        if valid_genome_cond:
            valid_genomes.append(new_entry)
            valid_genome_names.append(subdir_name)
            valid_links.append(new_entry[-1])
            valid_genomes_by_name[subdir_name] = new_entry

downloaded_assemblies = set(subdir_list)
known_assemblies = set(valid_genome_names)
known_links = set(valid_links)
missing_assemblies = known_assemblies.difference(downloaded_assemblies)

print("\nNumber of listed complete chromosomes/genomes:")
print(len(known_assemblies))

print("\nNumber of unique download links:")
print(len(known_links))

print("\nNumber of downloaded genomes:")
print(len(downloaded_assemblies))

print("\nNumber of missing genomes:")
print(len(missing_assemblies))



