# For quickly browsing the prokaryotic database list and 
# checking certain properties I'm interested in
import sys
import os
import csv

# Locate the open metadata files for reading
summary_data = sys.argv[1]

if 'genbank' in summary_data.lower():
    data_type = 'GenBank'
elif 'refseq' in summary_data.lower():
    data_type = 'RefSeq'
else:
    raise Exception("ERROR: must specify either 'genbank' or 'refseq' in input filename!")

summary_file = open(summary_data, 'r')
print(summary_data)

# Initialize storage structures
valid_genome_entries = []
field_names = []

# Process entries from the RefSeq/Genbank prokaryotic genome database metadata file
for index, line in enumerate(summary_file.readlines()):
    # Skip first line
    if index == 0:
        continue


    entry = line.split('\t')

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

    # Generate data sorting filters by criteria
    complete_status = (completion_status == 'complete genome')
    chromosome_status = (completion_status == 'chromosome')
    good_gen_assembly = complete_status or chromosome_status

    # Apply major filters to dataset
    if index == 1:
        new_entry.insert(-2, "Annotation Data")
        new_entry[-2] = "Assembly Directory"
        new_entry[-1] = "Rsync Path"
        field_names = new_entry
        print(field_names)

    if good_gen_assembly:
        valid_genome_entries.append(new_entry)

print(len(valid_genome_entries))

### DATABASE DOWNLOADER ######################
# Download command format
cmd_format = "rsync -r {} {}"


# Define output root directory by cmd line arg
root_dir = sys.argv[2].rstrip('/')
if not os.path.isdir(root_dir):
    os.mkdir(root_dir)

genome_subdir = root_dir + '/{}'.format(data_type)
if not os.path.isdir(genome_subdir):
    os.mkdir(genome_subdir)


# Download genome directories via Rsync
for entry in valid_genome_entries[:314]:
    rsync_cmd = cmd_format.format(entry[-1], genome_subdir)
    os.system(rsync_cmd)


# Annotation availability check function
def check_annotation(entry, all_genomes_subdir):
    # Check if genome has feature annotation data available
    annotation_file = "{}_feature_table.txt.gz".format(entry[-2])
    annotation_file = "{}/{}/{}".format(all_genomes_subdir, entry[-2], annotation_file)
    print(annotation_file)

    if os.path.isfile(annotation_file):
        entry.insert(-2, "yes")
    else:
        entry.insert(-2, "no")
    
    return(entry)


# Create summary file
summary_path = root_dir + '/{}_bacterial_genomes_summary.csv'.format(data_type)
with open(summary_path, 'w') as summary_file:
    csvwriter = csv.writer(summary_file)
    csvwriter.writerow(field_names)

    for entry in valid_genome_entries:
        entry = check_annotation(entry, genome_subdir)
        csvwriter.writerow(entry)


