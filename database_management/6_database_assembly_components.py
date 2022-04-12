# TITLE: GenBank/RefSeq Genome Assembly Individual Dataset Component Count
# AUTHOR: Katavga
# DATE: 11-Apr-2022

desc = '''
PURPOSE: Given a root directory containing a series of GenBank/RefSeq
assembly folders, check each for the existence of 5 key UNCOMPRESSED genome
files. After all assembly directories have been searched, output seven lists
as textfiles, each one corresponding to one of the genome file type, listing 
which assemblies contain the type of file corresponding to that specific list.
For instance, any assembly listed in the cds_from_genomic_list.txt file would
correspond to an assembly which has a *_cds_from_genomic.txt file (note that
the file being checked for is UNCOMPRESSED). These lists will be used to sort
datasets prior to their being placed into custom BLAST databases.
'''

import os, argparse

# COMMAND LINE LOGIC
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('--database_root_dir',  metavar='root_dir_path', type=str, required=True, action='store',
                    help='Root directory where genome assemblies with DECOMPRESSED component files exist.')
parser.add_argument('--output_dir',  metavar='output_dir_path', type=str, required=True, action='store',
                    help='Directory in which to output summary files.')
args = vars(parser.parse_args())


# Given an assembly subdirectory, return hashmap detailing
# which component files it does and does not have
def genome_dataset_components(genome_subdir):

    # ASSUMES FILES BEING SEARCHED FOR ARE NOT COMPRESSED!!
    target_file_types = [
        'cds_from_genomic.fna',
        'protein.faa',
        'rna_from_genomic.fna',
        'translated_cds.faa',
        'genomic.fna' # used mostly as a positive control for verifying my code works
    ]

    components_hashmap = dict()
    target_file_paths = []
    target_basenames = []

    # Build basenames of suffixed target files
    for suffix in target_file_types:
        assembly_name = genome_subdir.name
        file_basename = "{}_{}".format(assembly_name, suffix)
        target_basenames.append((suffix, file_basename))

    # Build full file path for each target file
    for (suffix, basename) in target_basenames:
        path = "{}/{}".format(genome_subdir.path, basename)
        target_file_paths.append((suffix, path))

    # Check whether each of the target genome files exists w/i the directory
    for (suffix, path) in target_file_paths:
        if os.path.isfile(path):
            components_hashmap[suffix] = True
        else:
            components_hashmap[suffix] = False
            
    return(components_hashmap)


# Builds lists of which genomes have data available for each
# of the listed target file types
def list_of_assemblies_by_component(list_of_genome_subdirs):

    output_hashmap = dict()
    cds_genomic_available = []
    protein_available = []
    rna_genomic_available = []
    cds_translated_available = []
    genome_available = []

    target_file_types = [
        'cds_from_genomic.fna',
        'protein.faa',
        'rna_from_genomic.fna',
        'translated_cds.faa',
        'genomic.fna' # used mostly as a positive control for verifying my code works
    ]
    
    for (index, assembly) in enumerate(list_of_genome_subdirs):

        msg = "{}/{}".format(index+1, len(list_of_genome_subdirs))
        print(msg)

        components_table = genome_dataset_components(assembly)

        for target in target_file_types:
            present = components_table[target]
            asm_name = assembly.name

            if present:
                if target == 'cds_from_genomic.fna':
                    cds_genomic_available.append(asm_name)
                
                elif target == 'protein.faa':
                    protein_available.append(asm_name)

                elif target == 'rna_from_genomic.fna':
                    rna_genomic_available.append(asm_name)
                
                elif target == 'translated_cds.faa':
                    cds_translated_available.append(asm_name)

                elif target == 'genomic.fna':
                    genome_available.append(asm_name)

    output_hashmap['cds_from_genomic.fna'] = cds_genomic_available
    output_hashmap['protein.faa'] = protein_available
    output_hashmap['rna_from_genomic.fna'] = rna_genomic_available
    output_hashmap['translated_cds.faa'] = cds_translated_available
    output_hashmap['genomic.fna'] = genome_available
    
    return(output_hashmap)

### RUNTIME ###
def main():
    root_dir = args['database_root_dir'].rstrip('/')
    output_dir = args['output_dir'].rstrip('/')
    subdirs = [subdir for subdir in os.scandir(root_dir) if os.path.isdir(subdir)]

    target_file_types = [
        'cds_from_genomic.fna',
        'protein.faa',
        'rna_from_genomic.fna',
        'translated_cds.faa',
        'genomic.fna' # used mostly as a positive control for verifying my code works
    ]

    components_table_list = list_of_assemblies_by_component(subdirs)

    for target_type in target_file_types:

        basename = os.path.splitext(target_type)[0]
        basename = "list_{}.txt".format(basename)
        output_path = "{}/{}".format(output_dir, basename)

        print("Number of Assemblies with *_{} Available: {}".format(target_type, len(components_table_list[target_type])))

        with open(output_path, 'w') as list_file:
            lines = '\n'.join(components_table_list[target_type])
            list_file.write(lines)


main()
