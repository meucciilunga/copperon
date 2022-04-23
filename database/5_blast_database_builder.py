# append all prokaryotic *.fna genome files into a single file,
# then convert that super-concatenated file into a formal BLAST database
# via makeblastdb program in BLAST+

import shutil
import os

def build_assembly_file_src_paths(assembly_subroot_dir):

    # Storage logistics
    target_file_paths = []

    # List of assembly directories
    for subdir in os.scandir(assembly_subroot_dir):
        file_path = "{}/{}_genomic.fna".format(subdir.path, subdir.name)
        target_file_paths.append(file_path)

    return(target_file_paths)

### RUNTIME #######
def main():
    annotated_db = '/run/media/katanga/SSD2/organized_prokaryotic_database/annotated'
    unannotated_db = '/run/media/katanga/SSD2/organized_prokaryotic_database/unannotated'
    output_file_path = '/run/media/katanga/SSD2/organized_prokaryotic_database/blast_db/prokaryotic_genome_database.fna'

    full_list_of_src_paths = build_assembly_file_src_paths(annotated_db)
    full_list_of_src_paths.extend(build_assembly_file_src_paths(unannotated_db))

    total_size = 0
    with open(output_file_path, 'wb') as output_file:

        for index, path in enumerate(full_list_of_src_paths):
            progress_msg = "{}/{}".format(index+1, len(full_list_of_src_paths))
            total_size = total_size + os.path.getsize(path)
            
            with open(path, 'rb') as fd:
                shutil.copyfileobj(fd, output_file)
            print(progress_msg)
    
    print("TOTAL SIZE: {} bytes".format(total_size))

main()

# Command used for putting super-concatenated file into a formal blastdb
#/home/katanga/ncbi-blast-2.13.0+-src/c++/ReleaseMT/bin/makeblastdb -in '/run/media/katanga/SSD2/organized_prokaryotic_database/blast_db/prokaryotic_genome_database.fna'\
#-title "Prokaryotic_Genome_Database" -dbtype 'nucl' -out 'prok_genome_db' -max_file_sz 4GB