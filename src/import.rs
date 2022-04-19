#![allow(dead_code)]
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

fn import_database_summary(summary_file_path: PathBuf) -> HashMap<String, AssemblyMetadata> {

    // Storage logistics
    let mut database_metadata_table: HashMap<String, AssemblyMetadata> = HashMap::new(); 

    // Read-in Metadata file
    let mut file = File::open(summary_file_path).expect("ERROR: could not open metadata summary file!");
    let mut contents = String::new();
    file.read_to_string(&mut contents).expect("ERROR: could not read string data from metadata summary file!");
    
    // Parse lines from Metadata file
    let metadata_content = contents.split('\n');
    let mut metadata_entries: Vec<String> = vec![];
    for (index, item) in metadata_content.enumerate() {

        // Skip first line of metadata file!
        if index == 0 {
            continue;
        }

        // Store entries
        metadata_entries.push(item.to_string());
    }

    // Filter off last entry, since it is just an empty string; aka, remove parsing artifact
    let modified_metadata_entries = metadata_entries.into_iter()
                                                    .filter(|s| !s.eq(&"".to_string()))
                                                    .collect::<Vec<String>>();

    // Parse individual data entries from file
    for line in modified_metadata_entries {
        let entries: Vec<String> = line.split('\t')
                                       .map(|s| s.chars().collect::<String>())
                                       .collect();

        let assembly_entry = AssemblyMetadata::new(entries);
        database_metadata_table.insert(assembly_entry.subdir_name.clone(), assembly_entry);
    }
    
    // Map assembly metadata (val) against assembly directory name (key)
    database_metadata_table
}

fn parse_genome_sequence(genome_name: String, genome_fasta_file: PathBuf) -> GenomeSequence {
    
    // Read-in Metadata file
    let mut file = File::open(genome_fasta_file).expect("ERROR: could not open metadata summary file!");
    let mut contents = String::new();
    file.read_to_string(&mut contents).expect("ERROR: could not read string data from metadata summary file!");

    // Parse fasta data into individual files
    let lines_of_sequence = contents.split('\n');
    let mut parsed_sequence: Vec<String> = vec![];
    let mut definition_line_indicies: Vec<usize> = vec![];
    let mut sequence_bounds: Vec<(usize, usize)> = vec![];
    let mut replicon_data: Vec<RepliconSequence> = vec![];

    // Determine replicon sequence bounds by determining location of 
    // every definition line by looking for '>' character
    for (index, line) in lines_of_sequence.enumerate() {

        // Final bound is defined by line with no characters at the end of every file
        if line.contains('>') || line.len() == 0 {
            definition_line_indicies.push(index);
        }
        parsed_sequence.push(line.to_string());
    }

    // Define the replicon sequence boundaries
    let start_indicies = definition_line_indicies.iter();
    let end_indicies = definition_line_indicies[1..].iter();
    for (&start, &end) in start_indicies.zip(end_indicies) {
        let new_bound = (start+1, end);
        sequence_bounds.push(new_bound);
    }

    // Parse out individual subsequences from FASTA data via bound indicies
    for (start, end) in sequence_bounds {
        let replicon_sequence = (&parsed_sequence[start..end].join("")).clone();

        // Parse out genomic accession identifier from definition
        let replicon_accession = parsed_sequence[start-1].clone()
                                                         .split(' ')
                                                         .map(|s| s.to_string())
                                                         .collect::<Vec<String>>()[0]
                                                         .clone()
                                                         .chars()
                                                         .filter(|&c| c != '>')
                                                         .collect::<String>();

        let new_replicon = RepliconSequence {
            replicon_accession,
            replicon_sequence,
        };

        replicon_data.push(new_replicon);
    }

    // Wrap data into a genome sequence
    GenomeSequence {
        genome_name,
        genomic_elements: replicon_data
    }
}

#[derive(PartialEq, Debug, Clone)]
struct AssemblyMetadata {
    assembly_accession:     String,
    biosample_accession:    String,
    ref_seq_cat:            String,
    tax_id:                 String,
    species_tax_id:         String,
    organism_name:          String,
    infraspecifc_name:      String,
    assembly_level:         String,
    assembly_name:          String,
    ftp_location:           String,
    subdir_name:            String,
    ref_seq_status:         String,
}

impl AssemblyMetadata {
    fn new(entry_data: Vec<String>) -> AssemblyMetadata {

        let assembly_accession = entry_data[0].to_string();
        let assembly_name = entry_data[15].to_string();
        let subdir_name = format!("{}_{}", assembly_accession, assembly_name);

        AssemblyMetadata {
            assembly_accession,
            biosample_accession:    entry_data[2].to_string(),
            ref_seq_cat:            entry_data[4].to_string(),
            tax_id:                 entry_data[5].to_string(),
            species_tax_id:         entry_data[6].to_string(),
            organism_name:          entry_data[7].to_string(),
            infraspecifc_name:      entry_data[8].to_string(),
            assembly_level:         entry_data[11].to_string().to_lowercase(),
            assembly_name,
            ftp_location:           entry_data[19].to_string(),
            subdir_name,
            ref_seq_status:         entry_data[20].to_string()
        }
    }
}

struct GenomeSequence {
    genome_name: String,
    genomic_elements: Vec<RepliconSequence>
}

#[derive(PartialEq, Debug, Clone)]
struct RepliconSequence {
    replicon_accession: String,
    replicon_sequence: String,
}

struct GeneralAnnotationData {}

struct CopYAnnoationData {}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use super::*;

    #[test]
    fn test_database_metadata_import_1() {
        let test_file = "database_management/data_assets/assembly_summary_refseq.txt";
        let test_file = PathBuf::from(test_file);
        let database_metadata = import_database_summary(test_file);

        let expected_answer = AssemblyMetadata { 
            assembly_accession:     "GCF_001735525.1".to_string(), 
            biosample_accession:    "SAMN05384437".to_string(), 
            ref_seq_cat:            "representative genome".to_string(), 
            tax_id:                 "23".to_string(), 
            species_tax_id:         "23".to_string(),
            organism_name:          "Shewanella colwelliana".to_string(),
            infraspecifc_name:      "strain=CSB03KR".to_string(),
            assembly_level:         "scaffold".to_string(),
            assembly_name:          "ASM173552v1".to_string(),
            ftp_location:           "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/735/525/GCF_001735525.1_ASM173552v1".to_string(),
            subdir_name:            "GCF_001735525.1_ASM173552v1".to_string(),
            ref_seq_status:         "".to_string()
        };

        let test_key = "GCF_001735525.1_ASM173552v1".to_string();
        let actual_answer = (*database_metadata.get(&test_key).unwrap()).clone();
        println!("{:?}", expected_answer);

        assert_eq!(expected_answer, actual_answer)
    }

    #[test]
    fn test_database_metadata_import_2() {
        let test_file = "database_management/data_assets/assembly_summary_genbank.txt";
        let test_file = PathBuf::from(test_file);
        let database_metadata = import_database_summary(test_file);

        let expected_answer = AssemblyMetadata { 
            assembly_accession:     "GCA_001465115.1".to_string(), 
            biosample_accession:    "SAMN04296138".to_string(), 
            ref_seq_cat:            "na".to_string(), 
            tax_id:                 "332949".to_string(), 
            species_tax_id:         "332949".to_string(),
            organism_name:          "Enterococcus silesiacus".to_string(),
            infraspecifc_name:      "strain=LMG 23085".to_string(),
            assembly_level:         "complete genome".to_string(),
            assembly_name:          "ASM146511v1".to_string(),
            ftp_location:           "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/465/115/GCA_001465115.1_ASM146511v1".to_string(),
            subdir_name:            "GCA_001465115.1_ASM146511v1".to_string(),
            ref_seq_status:         "missing tRNA genes".to_string()
        };

        let test_key = "GCA_001465115.1_ASM146511v1".to_string();
        let actual_answer = (*database_metadata.get(&test_key).unwrap()).clone();
        println!("{:?}", expected_answer);

        assert_eq!(expected_answer, actual_answer)
    }

    #[test]
    fn parse_genome_sequence_1() {
        let test_file = "test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let genome = parse_genome_sequence("GCF_000009725.1_ASM972v1".to_string(), test_file);
    }
}