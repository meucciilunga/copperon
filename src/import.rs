#![allow(dead_code)]
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, BufRead};
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

fn parse_genome_sequence(assembly_name: String, genome_fasta_file: PathBuf) -> GenomeSequence {
    
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
        assembly_name,
        genomic_elements: replicon_data
    }
}

fn parse_genome_annotation(genome_annotation_file: PathBuf) -> Vec<AnnotationEntry> {

    // Storage logistics
    let mut annotation_data: Vec<AnnotationEntry> = Vec::new();
      
    // Read-in Metadata file
    let file = File::open(genome_annotation_file).expect("ERROR: could not open annotation data file!");
    let file_lines = io::BufReader::new(file).lines();

    // Parse each line in the file into AnnotationEntry
    for line in file_lines {

        // Skip any line that starts with '#'
        let parse_error = "ERROR: could not unwrap line data.";
        if line.as_ref().expect(parse_error).chars().rev().last() == Some('#') {
            continue;
        }

        // Parse every data entry line
        if let Ok(data) = line {
            let entries = data.split('\t').map(|s| s.to_string()).collect::<Vec<String>>();
            let new_annotation = AnnotationEntry::from_split_line_vec(entries);
            annotation_data.push(new_annotation)
        }
    }

    annotation_data
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

#[derive(PartialEq, Debug, Clone)]
struct GenomeSequence {
    assembly_name: String,
    genomic_elements: Vec<RepliconSequence>
}

#[derive(PartialEq, Debug, Clone)]
struct RepliconSequence {
    replicon_accession: String,
    replicon_sequence: String,
}

#[derive(PartialEq, Debug, Clone)]
struct AnnotationEntry {
    genomic_accession: String,
    source: String,
    feature_type: String,
    ord_start_index: usize,
    ord_end_index: usize,
    replicon_strand: String,
    attributes: HashMap<String, String>,
}

impl AnnotationEntry {

    // Parse annotation entry from an input vector deriving from a split operation
    fn from_split_line_vec(input_vec: Vec<String>) -> AnnotationEntry {

        let start_index_err = "ERROR: could not parse annotation start index";
        let end_index_err = "ERROR: could not parse annotation start index";

        let mut attributes: HashMap<String, String> = HashMap::new();
        let raw_attributes = input_vec[8].split(';');

        for item in raw_attributes {
            let new_pair = item.split('=')
                .map(|s| s.to_string())
                .collect::<Vec<String>>();

            // Update attributes hash_table; assumes no duplicate entries
            attributes.insert(new_pair[0].clone(), new_pair[1].clone());
        }

        AnnotationEntry {
            genomic_accession:   input_vec[0].clone(),
            source:             input_vec[1].clone(),
            feature_type:       input_vec[2].clone(),
            ord_start_index:    input_vec[3].parse::<usize>().expect(start_index_err),
            ord_end_index:      input_vec[4].parse::<usize>().expect(end_index_err),
            replicon_strand:    input_vec[6].clone(),
            attributes,
        }
    }
}

struct BlastDerivedAnnotationData {}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::fs::{self, File};
    use std::io::Read;
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

    // For parsing genome sequence from preprocessed test files
    fn parse_confirmation_genome_dir(dir_path: PathBuf) -> GenomeSequence {
        let assembly_name = dir_path.file_name().unwrap().to_str().unwrap().to_string();
        let paths = fs::read_dir(dir_path).expect("ERROR: could not open genomic testing file directory.");
        let mut replicons: Vec<RepliconSequence> = Vec::new();

        for path in paths {
            let file_path = path.unwrap().path();
            let replicon_accession = file_path.file_name().unwrap().to_str().unwrap().to_string();

            let mut file = File::open(file_path).expect("ERROR: could not open replicon sequence file!");
            let mut replicon_sequence = String::new();

            file.read_to_string(&mut replicon_sequence).expect("ERROR: could not read sequence data from replicon file!");

            let new_replicon = RepliconSequence {
                replicon_accession,
                replicon_sequence,
            };

            replicons.push(new_replicon);
        }

        GenomeSequence {
            assembly_name,
            genomic_elements: replicons,
        }
    }

    #[test]
    fn test_parse_genome_sequence_1() {
        let test_file = "test_assets/GCA_000152665.1_ASM15266v1/GCA_000152665.1_ASM15266v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence("GCA_000152665.1_ASM15266v1".to_string(), test_file);

        let confirmation_dir = PathBuf::from("test_assets/preprocessed_test_genomes/GCA_000152665.1_ASM15266v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_eq!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_genome_sequence_2() {
        let test_file = "test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence("GCF_000009725.1_ASM972v1".to_string(), test_file);

        let confirmation_dir = PathBuf::from("test_assets/preprocessed_test_genomes/GCF_000009725.1_ASM972v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_eq!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_genome_sequence_3() {
        let test_file = "test_assets/GCF_014107515.1_ASM1410751v1/GCF_014107515.1_ASM1410751v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence("GCF_014107515.1_ASM1410751v1".to_string(), test_file);

        let confirmation_dir = PathBuf::from("test_assets/preprocessed_test_genomes/GCF_014107515.1_ASM1410751v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_eq!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_genome_sequence_4() {
        let test_file = "test_assets/GCF_016889785.1_ASM1688978v1/GCF_016889785.1_ASM1688978v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence("GCF_016889785.1_ASM1688978v1".to_string(), test_file);

        let confirmation_dir = PathBuf::from("test_assets/preprocessed_test_genomes/GCF_016889785.1_ASM1688978v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_eq!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_genome_sequence_5() {
        let test_file = "test_assets/GCF_016889785.1_ASM1688978v1/GCF_016889785.1_ASM1688978v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence("GCF_016889785.1_ASM1688978v1".to_string(), test_file);

        let confirmation_dir = PathBuf::from("test_assets/preprocessed_test_genomes/GCF_000009725.1_ASM972v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_ne!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_genome_sequence_6() {
        let test_file = "test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence("GCF_000009725.1_ASM972v1_missing_files".to_string(), test_file);

        let confirmation_dir = PathBuf::from("test_assets/preprocessed_test_genomes/GCF_000009725.1_ASM972v1_missing_files");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_ne!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_annotation_1() {
        let test_file = "test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.gff";
        let test_file = PathBuf::from(test_file);
        let parsed_annotation_entries = parse_genome_annotation(test_file);

        // Manually defined database entries
        let expected_entry_1 = AnnotationEntry {
            genomic_accession: "NC_000911.1".to_string(),
            source: "Protein Homology".to_string(),
            feature_type: "CDS".to_string(),
            ord_start_index: 3065617,
            ord_end_index: 3066207,
            replicon_strand: "-".to_string(),
            attributes: {
                let mut test_hashmap: HashMap<String, String> = HashMap::new();
                let test_keys = vec!["ID", "Parent", "Dbxref", "Name", "gbkey", 
                                 "inference", "locus_tag", "product", "protein_id", "transl_table"];
                let test_vals = vec!["cds-WP_010873937.1", "gene-SGL_RS16115", "Genbank:WP_010873937.1",
                                 "WP_010873937.1", "CDS", "COORDINATES: similar to AA sequence:RefSeq:WP_011153858.1",
                                 "SGL_RS16115", "DUF305 domain-containing protein", "WP_010873937.1", "11"];

                let test_keys = test_keys.iter();
                let test_vals = test_vals.iter();

                for (key, val) in test_keys.zip(test_vals) {
                    test_hashmap.insert(key.to_string(), val.to_string());
                }

                test_hashmap
            },
        };

        let expected_entry_2 = AnnotationEntry {
            genomic_accession: "NC_005232.1".to_string(),
            source: "RefSeq".to_string(),
            feature_type: "gene".to_string(),
            ord_start_index: 9520,
            ord_end_index: 10167,
            replicon_strand: "-".to_string(),
            attributes: {
                let mut test_hashmap: HashMap<String, String> = HashMap::new();

                let test_keys = vec!["ID", "Name", "gbkey", "gene_biotype",
                                     "locus_tag", "old_locus_tag"];
                let test_vals = vec!["gene-SGL_RS01415", "SGL_RS01415", "Gene", "protein_coding",
                                     "SGL_RS01415", "sll6010"];

                let test_keys = test_keys.iter();
                let test_vals = test_vals.iter();

                for (key, val) in test_keys.zip(test_vals) {
                    test_hashmap.insert(key.to_string(), val.to_string());
                }

                test_hashmap
            },
        };

        // Pull corresponding entries taken from automatic parser
        let mut entry_comparison: HashMap<usize, AnnotationEntry> = HashMap::new();
        for item in parsed_annotation_entries {

            let entry_1 = item.ord_start_index == 3065617 && item.genomic_accession.eq("NC_000911.1") && item.feature_type.eq("CDS");
            let entry_2 = item.ord_start_index == 9520 && item.genomic_accession.eq("NC_005232.1") && item.feature_type.eq("gene");

            if entry_1 || entry_2 {
                entry_comparison.insert(item.ord_start_index, item);
            }
        }

        // Compare parsed entries against manual entries
        assert_eq!(expected_entry_1, *entry_comparison.get(&expected_entry_1.ord_start_index).unwrap());
        assert_eq!(expected_entry_2, *entry_comparison.get(&expected_entry_2.ord_start_index).unwrap());
    }
}