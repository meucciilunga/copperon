# TITLE: GenBank/RefSeq Genome Assembly Database Organizer
# AUTHOR: Katavga
# DATE: 11-Apr-2022

desc = '''
PURPOSE: Given a root directory containing a pair of GenBank/RefSeq
database folders and a list of assembly names, locate every assembly
directory in the list and copy it and 9 specific files over to a new
root directory. The new root directory can then be used to run the
CopY-Operon-Search program, which requires that all genome assemblies
requiring processing have their corresponding directory placed in the
same root directory.
'''

import os, shutil, argparse

# COMMAND LINE LOGIC
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('--database_root_dir',  metavar='root_dir_path', type=str, required=True, action='store',
                    help='Root directory where both GenBank and RefSeq sub-root directories exist')
parser.add_argument('--output_dir',  metavar='output_dir_path', type=str, required=True, action='store',
                    help='New root directory where user would like to place the cloned assembly directories.')
parser.add_argument('--assembly_list',  metavar='assembly_list_path', required=True, type=str, action='store',
                    help="Text file containing list of assemblies to be copied. One assembly name per line.")
args = vars(parser.parse_args())


def build_list_of_subdirs(database_root_dir, assembly_list):
    genbank_assemblies = []
    refseq_assemblies = []

    # Sort assembly names by database type
    for asm in assembly_list:

        # First 4 chars give info about which database asm is from
        genbank_status = ('GCA_' == asm[0:4])
        refseq_status = ('GCF_' == asm[0:4])

        if genbank_status:
            genbank_assemblies.append(asm)
        elif refseq_status:
            refseq_assemblies.append(asm)

    # Build paths by database type
    genbank_subdir = '{}/GenBank'.format(database_root_dir)
    genbank_paths = ['{}/{}'.format(genbank_subdir, asm) for asm in genbank_assemblies]

    refseq_subdir = '{}/RefSeq'.format(database_root_dir)
    refseq_paths = ['{}/{}'.format(refseq_subdir, asm) for asm in refseq_assemblies]

    all_assembly_paths = genbank_paths + refseq_paths

    if len(all_assembly_paths) != len(assembly_list):
        raise Exception('ERROR: assembly list contains incorrectly prefixed file!')

    return(all_assembly_paths)


def build_assembly_file_src_paths(assembly_subdir):

    target_file_paths = []

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

    # Build basenames of suffixed target files
    for suffix in target_file_types:
        assembly_name = os.path.basename(assembly_subdir)
        file_basename = "{}_{}".format(assembly_name, suffix)
        target_basenames_list.append(file_basename)

    # Build full src paths
    for basename in target_basenames_list:
        file_path = "{}/{}".format(assembly_subdir, basename)
        target_file_paths.append(file_path)

    return(target_file_paths)


def clone_assembly_dir(assembly_subdir_path, output_dir):
    asm_basename = os.path.basename(assembly_subdir_path)
    dest = "{}/{}".format(output_dir, asm_basename)

    os.makedirs(dest)

    # Copy files over from directory
    assembly_file_paths = build_assembly_file_src_paths(assembly_subdir_path)
    for src in assembly_file_paths:

        # Some files don't exist; if a particular assembly component file doesn't
        # exist, just move on to the next file
        if os.path.isfile(src):
            shutil.copy2(src, dest)
        else:
            continue

    # Locate then Decompress all *.gz files
    new_files = os.scandir(dest)
    compressed_files = []
    for new_file in new_files:
        if new_file.path.endswith('.gz'):
            compressed_files.append(new_file)

    cmd_format = "gunzip -d {}"
    for compressed in compressed_files:
        cmd = cmd_format.format(compressed.path)
        os.system(cmd)


### RUNTIME ###
def main():
    
    # Parse cmd line arguments
    database_root_dir = args['database_root_dir']
    output_dir = args['output_dir']
    assembly_list_path = args['assembly_list']

    # Parse list of assemblies from file
    with open(assembly_list_path, 'r') as asm_list_file:
        list_of_assembly_names = [name.rstrip('\n') for name in asm_list_file.readlines()]

    # Build assembly paths then copy to new output directory
    asm_subdirs = build_list_of_subdirs(database_root_dir, list_of_assembly_names)
    for (index, assembly_dir) in enumerate(asm_subdirs):
        msg = "Progress: {}/{}".format(index+1, len(asm_subdirs))
        print(msg)

        clone_assembly_dir(assembly_dir, output_dir)

main()