# TITLE: GenBank/RefSeq Genome Assembly Individual Dataset Completeness Analyzer
# AUTHOR: Katavga
# DATE: 8-Apr-2022

desc = '''
PURPOSE: Given a root directory containing a series of GenBank/RefSeq
assembly folders, check for the existence of 9 key genome files. Any
genome directory corresponding to a genomic dataset with all 9 key files 
available is considered a 'complete genome' or a 'complete genomic dataset'
for the purposes of this project. User has ability to specifiy whether
they want the program to attempt to redownload any located incomplete
genomes via RSYNC.
'''

import os, argparse

# COMMAND LINE LOGIC
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('-i',  metavar='root_dir', type=str, required=True, action='store',
                    help='Root directory where all genome assembly folders to be processed exist as subdirectories.')
parser.add_argument('-o',  metavar='output_file_path', type=str, required=True, action='store',
                    help='File path where user would like to store the *.txt list of incomplete assemblies.')
parser.add_argument('--attempt_redownload',  metavar='attempt_redownload', nargs=2, type=str, action='store',
                    help="Whether to retry downloading assembly directories for assemblies listed as 'incomplete'; \
                        when used, user must supply two arguments: first, a path to a genome metadata file, then \
                        a path to the directory where they would like the re-attempted files to be stored.")
args = vars(parser.parse_args())

# MAJOR FUNCTIONS

# Pull list of genome assemblies that have a 'complete' datasets
def genome_dataset_completeness(genome_subdir):

    # these ten files define a 'complete' dataset
    target_basenames_list = [
        'annotation_hashes.txt',
        'md5checksums.txt'
    ]

    target_file_types = [
        'assembly_report.txt',
        'cds_from_genomic.fna.gz',
        'feature_table.txt.gz',
        'genomic.fna.gz',
        'protein.faa.gz',
        'rna_from_genomic.fna.gz',
        'translated_cds.faa.gz'
    ]

    target_file_paths = []

    # Build basenames of suffixed target files
    for suffix in target_file_types:
        assembly_name = genome_subdir.name
        file_basename = "{}_{}".format(assembly_name, suffix)
        target_basenames_list.append(file_basename)

    # Build full file path for each target file
    for basename in target_basenames_list:
        path = "{}/{}".format(genome_subdir.path, basename)
        target_file_paths.append(path)

    # Check that each of the relevant genome files exists w/i the directory
    valid_genome = True
    for path in target_file_paths:
        if not os.path.isfile(path):
            valid_genome = False
            break

    return(valid_genome)


def print_bad_asm_to_file(output_path, bad_asm_list):
    with open(output_path, 'w') as f:
        f.write('\n'.join(bad_asm_list))


def parse_genome_ftp_metadata(genomes_metadata_path):
    with open(genomes_metadata_path, 'r') as metadata:
        download_links_by_subdir_name = dict()

        for index, item in enumerate(metadata.readlines()):

            # Skip first line of the file
            if index == 0:
                continue
            
            # Split lines by \t marker
            entry = item.split('\t')

            # Pull specific data for each entry
            ftp_location = entry[19]
            subdir_name = ftp_location.split('/')[-1] 

            # Build a hashmap of subdir name to subdir download_link
            download_links_by_subdir_name[subdir_name] = "rsync" + ftp_location.lstrip("https")

    return(download_links_by_subdir_name)


def redownload_incomplete_genomes(metadata_path, output_dir, incomplete_list):
    cmd_format = "rsync --copy-links --recursive --times --verbose '{}' '{}'" # >/dev/null" # hides output from Rsync
    genomes_to_links = parse_genome_ftp_metadata(metadata_path)
    download_links = []
    
    for asm in incomplete_list:
        link = genomes_to_links[asm]
        download_links.append(link)
        
    for (index, link) in enumerate(download_links):
        progress_msg = "{}/{}".format(index+1, len(download_links))
        print(progress_msg)

        cmd = cmd_format.format(link, output_dir)
        os.system(cmd)

        print('\n\n\n')


### RUNTIME ###
def main():
    root_dir = args['i']
    output_file_path = args['o']
    redownload_flag = False

    if args['attempt_redownload'] is not None:
        metadata_path = args['attempt_redownload'][0]
        output_dir_path = args['attempt_redownload'][1]
        redownload_flag = True
        if not os.path.isdir(output_dir_path):
            raise ValueError("Output directory provided for redownloads does not exist!") 
        

    subdirs = [subdir for subdir in os.scandir(root_dir) if os.path.isdir(subdir)]

    complete_genomes = []
    incomplete_genomes = []

    for (index, genome_dir) in enumerate(subdirs):
        
        status_msg = "{}/{}".format(index+1, len(subdirs))
        print(status_msg)

        complete_dataset = genome_dataset_completeness(genome_dir)
        if complete_dataset:
            complete_genomes.append(genome_dir.name)
        else:
            incomplete_genomes.append(genome_dir.name)
    
    print_bad_asm_to_file(output_file_path, incomplete_genomes)

    ### ATTEMPT REDOWNLOAD
    if redownload_flag:
        redownload_incomplete_genomes(metadata_path, output_dir_path, incomplete_genomes)

    ### PRINT SUMMARY OF ANALYSIS
    tot = len(subdirs)
    complete = len(complete_genomes)
    incomplete = len(incomplete_genomes)
    msgs = [
        'Complete Datasets: {}/{}'.format(complete, tot),
        'Incomplete Datasets: {}/{}'.format(incomplete, tot)
    ]

    print()
    for msg in msgs:
        print(msg)
    

main()