# TITLE: GenBank/RefSeq Genome Assembly Database Integrity Validation
# AUTHOR: Katavga
# DATE: 8-Apr-2022

# INPUTS: first cmd_line argv should be the path to a root directory
# containing a series of subdirectories composed of GenBank/RefSeq
# assembly directories.

# PURPOSE: Given a root directory containing a series of GenBank/RefSeq
# assembly folders each containing a standard NCBI Genome-Assembly md5 
# checksum file, pull from md5checksums.txt the path for each known 
# md5-validated file and the value of its corresponding MD5 checksum. 
# For every path pulled from the checksum file that corresponds to a 
# valid file within that particular genome subdirectory, check the validity 
# of the current file's MD5 checksum. Notably, not all files listed in the
# checksum file need to be present in the directory for this script to be 
# used; files that are listed in a md5checksums.txt but not currently available
# in a given genome directory when the script is run will simply be ignored.
# This allows this script to validate the MD5-values of relevant genome files
# right after they've been pulled from the NCBI FTP servers AND after the
# directories have been moved and reorganized in preparation for running
# the CopY-Regulon-Search application, which generally involves another
# move of the database files to a remote server.

import hashlib, os, sys

# Open md5checksums file and pull all values into a hashmap
def pull_md5_hashes(subdir_path):
    md5_file_basename = "md5checksums.txt"
    md5_path = "{}/{}".format(subdir_path, md5_file_basename)
    md5_table = dict()

    with open(md5_path, 'r') as md5checksum_file:
        lines = md5checksum_file.readlines()

        for line in lines:
            entry = line.rstrip('\n').split('  ') # two space chars ' ' split the columns of data
            key = subdir_path + entry[1].lstrip('.')
            val = entry[0]
            md5_table[key] = val

    return(md5_table)


# Calculate md5 of local file and compare it to known md5 value
def verify_file_md5(known_md5s, input_file_path):

    # Calculate md5 of local file
    with open(input_file_path, 'rb') as f:
        data = f.read()
        md5 = hashlib.md5(data).hexdigest()

    # Compare local md5 with saved md5
    hash_cmp = (md5 == known_md5s[input_file_path])

    return(hash_cmp)


# Given a RefSeq/GeneBank directory containing a set of genome data
# files, verify the MD5 checksums of ALL available files 
def genome_md5_check(subdir):

    # Check if genome assembly directory's md5 file exists
    md5_file_basename = "md5checksums.txt"
    md5_path = "{}/{}".format(subdir, md5_file_basename)

    # If md5checksum.txt does not exist in the directory, mark genome as invalid
    if not os.path.isfile(md5_path):
        return(False)

    # Otherwise, pulls paths and checksums from file and perform md5 validation
    # for each individual, relevant file w/i the genome directory
    md5_checksums = pull_md5_hashes(subdir)
    md5_genome_file_paths = md5_checksums.keys()
    available_genome_files = [path for path in md5_genome_file_paths if os.path.isfile(path)]
    # to make this script more flexible, we're only going to validate checksums on the files in
    # the genome directory that are currently available; if a file's path is listed in the md5
    # file but not found in the current dir (likely due to reorganizing the files in preparation
    # for genome search analysis), this script won't worry about validating its md5_checksum.

    genome_valid = True
    for path in available_genome_files:
        file_validity = verify_file_md5(md5_checksums, path)
        
        # if one invalid file is found, mark whole genome as invalid
        if not file_validity:
            genome_valid = False
            break

    return(genome_valid)


#### RUNTIME #########
root_dir = sys.argv[1]
subdirs = [subdir for subdir in os.scandir(root_dir) if os.path.isdir(subdir)]
good_md5 = []
bad_md5 = []

for (index, genome_dir) in enumerate(subdirs):

    # Live progress counter
    msg = "{}/{}".format(index+1, len(subdirs))
    print(msg)

    # Perform md5 validation for every genome assembly directory in the root directory
    valid = genome_md5_check(genome_dir.path)
    if valid:
        good_md5.append(genome_dir.name)
    else:
        bad_md5.append(genome_dir.name)

msg1 = "Good Genomes: {}/{}".format(len(good_md5), len(subdirs))
msg2 = "Bad Genomes: {}/{}".format(len(bad_md5), len(subdirs))
print(msg1)
print(msg2)

for asm in bad_md5:
    print(asm)