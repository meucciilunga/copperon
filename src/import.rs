use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, BufRead};
use std::path::PathBuf;

pub fn import_database_summary(summary_file_path: PathBuf) -> HashMap<String, AssemblyMetadata> {

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

pub fn parse_genome_sequence(assembly_name: &String, genome_fasta_file: &PathBuf) -> ProtoGenome {
    
    // Read-in genome sequence file
    let file = File::open(genome_fasta_file).expect("ERROR: could not open genome sequence file!");
    let file_lines = io::BufReader::new(file).lines();
    
    // Storage Logistics
    let mut parsed_sequence: Vec<String> = vec![];
    let mut definition_line_indicies: Vec<usize> = vec![];
    let mut sequence_bounds: Vec<(usize, usize)> = vec![];
    let mut replicon_data: Vec<ProtoReplicon> = vec![];

    // Determine replicon sequence bounds by determining location of 
    // every definition line by looking for '>' character
    for (index, line) in file_lines.enumerate() {

        // Specially target lines that start with '>'
        let parse_error = "ERROR: could not unwrap line data. (GENOME PARSER)";
        if line.as_ref().expect(parse_error).chars().rev().last() == Some('>') {
            definition_line_indicies.push(index);
        }

        if let Ok(data) = line {
            // Only definition lines should contain this character
            if data.contains('>') {
                // Push definition lines as they are
                parsed_sequence.push(data.to_string());
            } else {
                // Ensure lines containing sequence data 
                // are entirely composed of uppercase chars 
                // before adding them to the storage vector
                parsed_sequence.push(data.to_uppercase());
            }
        }

        // NOTE: the final step collects the lines of the FASTA file
        // into the storage vector, while ensuring both that all FASTA sequence  
        // data is in a capital letter format AND that case-sensitive text data
        // in the FASTA definition lines isn't altered, since data from the
        // definition line data needs to be used as is
    }

    // Define the final bound of FASTA file, given by the number of lines it has
    definition_line_indicies.push(parsed_sequence.len());

    // Formally define the sequence boundaries for every replicon in FASTA based on the indicies of each 
    // defn line; NOTE: REMEMBER THAT THE END BOUND (in this implementation) IS EXCLUSIVE AND NOT INCLUSIVE
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

        let new_replicon = ProtoReplicon {
            replicon_accession,
            replicon_sequence,
            proto_replicon_type: ProtoRepliconType::Plasmid,
            // Default all new replicons to 'plasmid' type, will
            // correct for this declaration a few lines down
        };

        replicon_data.push(new_replicon);
    }

    // Find length of the longest sequence in FASTA
    let mut max_replicon_length = 0;
    for replicon in replicon_data.iter() {
        if replicon.replicon_sequence.len() > max_replicon_length {
            max_replicon_length = replicon.replicon_sequence.len();
        }
    };

    // Adjust replicon type based on length of longest replicon
    for replicon in &mut replicon_data {
        if replicon.replicon_sequence.len() == max_replicon_length {
            replicon.proto_replicon_type = ProtoRepliconType::Chromosome;
        } else {
            replicon.proto_replicon_type = ProtoRepliconType::Plasmid;
        }
    }

    // Wrap all derived data into a genome sequence struct
    ProtoGenome {
        assembly_name: assembly_name.clone(),
        proto_replicons: replicon_data
    }
}

pub fn parse_genome_annotation(genome_annotation_file: &PathBuf) -> Vec<AnnotationEntry> {

    // Storage logistics
    let mut annotation_data: Vec<AnnotationEntry> = Vec::new();
      
    // Read-in Metadata file
    let file = File::open(genome_annotation_file).expect("ERROR: could not open annotation data file!");
    let file_lines = io::BufReader::new(file).lines();

    // Parse each line in the file into AnnotationEntry
    for line in file_lines {

        // Skip any line that starts with '#'
        let parse_error = "ERROR: could not unwrap line data. (GENOME PARSER)";
        if line.as_ref().expect(parse_error).chars().rev().last() == Some('#') {
            continue;
        }

        // Parse every data entry line
        if let Ok(data) = line {
            let entries = data.split('\t')
                              .map(|s| s.to_string())
                              .collect::<Vec<String>>();
            let new_annotation = AnnotationEntry::from_split_line_vec(entries);
            annotation_data.push(new_annotation)
        }
    }

    annotation_data
}

pub fn parse_annotations_from_blast_results(blast_result_file: PathBuf) -> Vec<BlastDerivedAnnotation> {
    
    // Storage logistics
    let mut blast_annotation_data: Vec<BlastDerivedAnnotation> = Vec::new();
    
    // Read-in Metadata file
    let file = File::open(blast_result_file).expect("ERROR: could not open BLAST results file!");
    let file_lines = io::BufReader::new(file).lines();

    // Parse each line in the file into AnnotationEntry
    for line in file_lines {

        // Skip any line that starts with '#'
        let parse_error = "ERROR: could not unwrap line data. (BLAST PARSER)";
        if line.as_ref().expect(parse_error).chars().rev().last() == Some('#') {
            continue;
        }

        // Parse every data entry line
        if let Ok(data) = line {
            let entries = data.split('\t').map(|s| s.to_string()).collect::<Vec<String>>();
            let new_blast_annotation = BlastDerivedAnnotation::from_split_line_vec(entries);
            blast_annotation_data.push(new_blast_annotation)
        }
    }

    blast_annotation_data
}

#[derive(PartialEq, Debug, Clone)]
pub enum ProtoRepliconType {
    Chromosome,
    Plasmid,
}

#[derive(PartialEq, Debug, Clone)]
pub struct AssemblyMetadata {
    pub assembly_accession:     String,
    pub biosample_accession:    String,
    pub ref_seq_cat:            String,
    pub tax_id:                 String,
    pub species_tax_id:         String,
    pub organism_name:          String,
    pub infraspecifc_name:      String,
    pub assembly_level:         String,
    pub assembly_name:          String,
    pub ftp_location:           String,
    pub subdir_name:            String,
    pub ref_seq_status:         String,
}

impl AssemblyMetadata {
    fn new(entry_data: Vec<String>) -> AssemblyMetadata {

        let assembly_accession = entry_data[0].to_string();
        let assembly_name = entry_data[15].to_string()
                                          .replace("  ", " ")
                                          .replace(" ", "_")
                                          .replace("#", "_")
                                          .replace("/", "_")
                                          .replace("(", "_")
                                          .replace(")", "_");
        // RefSeq/GenBank metadata files occasionally don't match assembly
        // directory names exactly due to use of special characters in assembly
        // name, so have to make these replacements in order to compensate for that ^
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
pub struct ProtoGenome {
    pub assembly_name: String,
    pub proto_replicons: Vec<ProtoReplicon>
}

#[derive(PartialEq, Debug, Clone)]
pub struct ProtoReplicon {
    pub replicon_accession: String,
    pub replicon_sequence: String,
    pub proto_replicon_type: ProtoRepliconType,
}

#[derive(PartialEq, Debug, Clone)]
pub struct AnnotationEntry {
    pub genomic_accession:  String,
    pub source_database:    String,
    pub feature_type:       String,
    pub start_index_ord:    usize,
    pub end_index_ord:      usize,
    pub replicon_strand:    char,
    pub attributes:         HashMap<String, String>,
}

impl AnnotationEntry {

    // Parse annotation entry from an input vector deriving from a split operation
    fn from_split_line_vec(input_vec: Vec<String>) -> AnnotationEntry {

        let start_index_err =   "ERROR: could not parse annotation start index";
        let end_index_err =     "ERROR: could not parse annotation end index";
        let strand_err =        "ERROR: could not parse strand information";

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
            genomic_accession:  input_vec[0].clone(),
            source_database:    input_vec[1].clone(),
            feature_type:       input_vec[2].clone(),
            start_index_ord:    input_vec[3].parse::<usize>().expect(start_index_err),
            end_index_ord:     input_vec[4].parse::<usize>().expect(end_index_err),
            replicon_strand:    input_vec[6].parse::<char>().expect(strand_err),
            attributes,
        }
    }
}

#[derive(PartialEq, Debug)]
pub struct BlastDerivedAnnotation {
    pub sseqid: String,
    pub sstart: usize,
    pub send: usize,
    pub match_length: usize,
    pub sframe: i32,
    pub pident: f64,
    pub evalue: f64,
    pub qseqid: String,
    pub qstart: usize,
    pub qend: usize,
}

impl BlastDerivedAnnotation {
    fn from_split_line_vec(input_vec: Vec<String>) -> BlastDerivedAnnotation {
        let start_index_err =   "ERROR: could not parse BLAST annotation start index";
        let end_index_err =     "ERROR: could not parse BLAST annotation start index";
        let length_err =        "ERROR: could not parse BLAST annotation 'match length' (column 4)";
        let reading_frame_err = "ERROR: could not parse BLAST annotation 'reading frame' (column 5)";
        let pident_err =        "ERROR: could not parse BLAST annotation 'pident' (column 6)";
        let eval_err =          "ERROR: could not parse BLAST annotation 'evalue' (column 7)";

        BlastDerivedAnnotation {
            sseqid:         input_vec[0].clone(),
            sstart:         input_vec[1].parse::<usize>().expect(start_index_err),
            send:           input_vec[2].parse::<usize>().expect(end_index_err),
            match_length:   input_vec[3].parse::<usize>().expect(length_err),
            sframe:         input_vec[4].parse::<i32>().expect(reading_frame_err),
            pident:         input_vec[5].parse::<f64>().expect(pident_err),
            evalue:         input_vec[6].parse::<f64>().expect(eval_err),
            qseqid:         input_vec[7].clone(),
            qstart:         input_vec[8].parse::<usize>().expect(start_index_err),
            qend:           input_vec[9].parse::<usize>().expect(end_index_err),
        }
    }
}















// UNIT TESTS
#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::fs::{self, File};
    use std::io::Read;
    use super::*;

    #[test]
    fn test_database_metadata_import_1() {
        let test_file = "database/data_assets/assembly_summary_refseq.txt";
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
        let test_file = "database/data_assets/assembly_summary_genbank.txt";
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
    fn parse_confirmation_genome_dir(dir_path: PathBuf) -> ProtoGenome {
        let assembly_name = dir_path.file_name().unwrap().to_str().unwrap().to_string();
        let paths = fs::read_dir(dir_path).expect("ERROR: could not open genomic testing file directory.");
        let mut replicons: Vec<ProtoReplicon> = Vec::new();

        for path in paths {
            let file_path = path.unwrap().path();
            let replicon_accession = file_path.file_name().unwrap().to_str().unwrap().to_string();

            let mut file = File::open(file_path).expect("ERROR: could not open replicon sequence file!");
            let mut replicon_sequence = String::new();

            file.read_to_string(&mut replicon_sequence).expect("ERROR: could not read sequence data from replicon file!");

            // temporarily set all replicon types to chromosome
            let proto_replicon_type = ProtoRepliconType::Plasmid;

            let new_replicon = ProtoReplicon {
                replicon_accession,
                replicon_sequence,
                proto_replicon_type,
            };

            replicons.push(new_replicon);
        }

        // Find length of longest replicon in genome
        let mut max_len = 0_usize;
        for replicon in replicons.iter() {
            if replicon.replicon_sequence.len() > max_len {
                max_len = replicon.replicon_sequence.len()
            }
        }

        // Adjust replicon_types based on maximum replicon length
        for replicon in &mut replicons {
            if replicon.replicon_sequence.len() == max_len {
                replicon.proto_replicon_type = ProtoRepliconType::Chromosome;
            } else {
                replicon.proto_replicon_type = ProtoRepliconType::Plasmid;
            }
        }

        ProtoGenome {
            assembly_name,
            proto_replicons: replicons,
        }
    }

    #[test]
    fn test_parse_genome_sequence_1() {
        let test_file = "tests/test_assets/GCA_000152665.1_ASM15266v1/GCA_000152665.1_ASM15266v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let test_genome = parse_genome_sequence(&"GCA_000152665.1_ASM15266v1".to_string(), &test_file);

        let confirmation_dir = PathBuf::from("tests/test_assets/preprocessed_test_genomes/GCA_000152665.1_ASM15266v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        // Validate that both methods parse to identical objects
        assert_eq!(test_genome, actual_genome);
        
        // Validate that the correct number of replicons are parsed
        assert_eq!(test_genome.proto_replicons.len(), 1);

        // Validate total length of genome against known value
        assert_eq!(test_genome.proto_replicons[0].replicon_sequence.len(), 4_109_689);
    }

    #[test]
    fn test_parse_genome_sequence_2() {
        let test_file = "tests/test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let test_genome = parse_genome_sequence(&"GCF_000009725.1_ASM972v1".to_string(), &test_file);

        let confirmation_dir = PathBuf::from("tests/test_assets/preprocessed_test_genomes/GCF_000009725.1_ASM972v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        let mut test_replicons: HashMap<String, ProtoReplicon> = HashMap::new();
        let mut actual_replicons: HashMap<String, ProtoReplicon> = HashMap::new();

        for item in test_genome.proto_replicons.into_iter() {
            test_replicons.insert(item.replicon_accession.clone(), item);
        }

        for item in actual_genome.proto_replicons.into_iter() {
            actual_replicons.insert(item.replicon_accession.clone(), item);
        }

        // This test is a little weird; basically, the program passes the test on manual review, 
        // but fails it when undergoing the automatic review in cargo; the reason for the 
        // auto-fail is that the order of the replicons in the test_genome vector differ slightly
        // from the order of the replicons in the actual_genome vector; 
        // this is due in part to differences in how they're imported: importing the actual_genome
        // via the file system puts the replicons in alpha-numerical order; importing via
        // the fasta uses the same order for the replicons as is listed in the fasta file. 
        // These may not be the same lol. To get around this for this test, and so I don't
        // have to butcher the rest of the codebase for such a minor problem, basically turn the 
        // replicon vector into a accession -> replicon struct hash_table, then run the comparison 
        // against the hash_table and against the assembly names seperately :) the result is the
        // same: we check that parsing a genome via the fasta function gives the proper result 
        assert_eq!(test_genome.assembly_name, actual_genome.assembly_name);
        assert_eq!(test_replicons, actual_replicons);

        // Validate total length of genome against known value
        let mut total_len = 0_usize;
        for genomic_element in test_replicons.values() {
            total_len += genomic_element.replicon_sequence.len()
        }
        assert_eq!(total_len, 3_947_019);
    }

    #[test]
    fn test_parse_genome_sequence_3() {
        let test_file = "tests/test_assets/GCF_014107515.1_ASM1410751v1/GCF_014107515.1_ASM1410751v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let test_genome = parse_genome_sequence(&"GCF_014107515.1_ASM1410751v1".to_string(), &test_file);

        let confirmation_dir = PathBuf::from("tests/test_assets/preprocessed_test_genomes/GCF_014107515.1_ASM1410751v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        // Validate that both methods parse to identical objects
        assert_eq!(test_genome, actual_genome);

        // Validate that the correct number of replicons are parsed
        assert_eq!(test_genome.proto_replicons.len(), 3);

        // Validate total length of genome against known value
        let mut total_len = 0_usize;
        for genomic_element in &test_genome.proto_replicons {
            total_len += genomic_element.replicon_sequence.len()
        }
        assert_eq!(total_len, 5_356_494);
    }

    #[test]
    fn test_parse_genome_sequence_4() {
        // NOTE: TEST INPUT IS ENTIRELY IN LOWERCASE CHARS, SO ALSO CHECKS IF LOWERCASE -> UPPERCASE ADJUSTMENT WORKS
        // SINCE CONFIRMATION INPUT IS ALL UPPERCASE
        let test_file = "tests/test_assets/GCF_016889785.1_ASM1688978v1/GCF_016889785.1_ASM1688978v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let test_genome = parse_genome_sequence(&"GCF_016889785.1_ASM1688978v1".to_string(), &test_file);

        let confirmation_dir = PathBuf::from("tests/test_assets/preprocessed_test_genomes/GCF_016889785.1_ASM1688978v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        // Validate that both methods parse to identical objects
        assert_eq!(test_genome, actual_genome);

        // Validate that the correct number of replicons are parsed
        assert_eq!(test_genome.proto_replicons.len(), 1);

        // Validate total length of genome against known value
        let mut total_len = 0_usize;
        for genomic_element in &test_genome.proto_replicons {
            total_len += genomic_element.replicon_sequence.len()
        }
        assert_eq!(total_len, 2_886_678);
    }

    #[test]
    fn test_parse_genome_sequence_5() {
        let test_file = "tests/test_assets/GCF_016889785.1_ASM1688978v1/GCF_016889785.1_ASM1688978v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let test_genome = parse_genome_sequence(&"GCF_016889785.1_ASM1688978v1".to_string(), &test_file);

        let confirmation_dir = PathBuf::from("tests/test_assets/preprocessed_test_genomes/GCF_000009725.1_ASM972v1");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_ne!(test_genome, actual_genome);
    }

    #[test]
    fn test_parse_genome_sequence_6() {
        let test_file = "tests/test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna";
        let test_file = PathBuf::from(test_file);
        let expected_genome = parse_genome_sequence(&"GCF_000009725.1_ASM972v1_missing_files".to_string(), &test_file);

        let confirmation_dir = PathBuf::from("tests/test_assets/preprocessed_test_genomes/GCF_000009725.1_ASM972v1_missing_files");
        let actual_genome = parse_confirmation_genome_dir(confirmation_dir);

        assert_ne!(expected_genome, actual_genome);
    }

    #[test]
    fn test_parse_annotation_1() {
        let test_file = "tests/test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.gff";
        let test_file = PathBuf::from(test_file);
        let parsed_annotation_entries = parse_genome_annotation(&test_file);

        // Manually defined database entries
        let expected_entry_1 = AnnotationEntry {
            genomic_accession: "NC_000911.1".to_string(),
            source_database: "Protein Homology".to_string(),
            feature_type: "CDS".to_string(),
            start_index_ord: 3065617,
            end_index_ord: 3066207,
            replicon_strand: '-',
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
            source_database: "RefSeq".to_string(),
            feature_type: "gene".to_string(),
            start_index_ord: 9520,
            end_index_ord: 10167,
            replicon_strand: '-',
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

            let entry_1 = item.start_index_ord == 3065617 && item.genomic_accession.eq("NC_000911.1") && item.feature_type.eq("CDS");
            let entry_2 = item.start_index_ord == 9520 && item.genomic_accession.eq("NC_005232.1") && item.feature_type.eq("gene");

            if entry_1 || entry_2 {
                entry_comparison.insert(item.start_index_ord, item);
            }
        }

        // Compare parsed entries against manual entries
        assert_eq!(expected_entry_1, *entry_comparison.get(&expected_entry_1.start_index_ord).unwrap());
        assert_eq!(expected_entry_2, *entry_comparison.get(&expected_entry_2.start_index_ord).unwrap());
    }

    #[test]
    fn test_parse_blast_1() {
        let test_file = "tests/test_assets/CopZ-TcrZ_blast_result.txt";
        let test_file = PathBuf::from(test_file);
        let parsed_blast_entries = parse_annotations_from_blast_results(test_file);

        let expected_1 = BlastDerivedAnnotation {
            sseqid: "NZ_OD940440.2".to_string(),
            sstart: 261732,
            send: 261932,
            match_length: 67,
            sframe: 3,
            pident: 44.776,
            evalue: 8.68e-10,
            qseqid: "WP_002294571.1".to_string(),
            qstart: 1,
            qend: 67,
        };

        let expected_2 = BlastDerivedAnnotation {
            sseqid: "NZ_CP064343.1".to_string(),
            sstart: 1870968,
            send: 1870762,
            match_length: 69,
            sframe: -3,
            pident: 66.667,
            evalue: 8.45e-15,
            qseqid: "WP_010718490.1".to_string(),
            qstart: 1,
            qend: 69,
        };

        
        assert_eq!(parsed_blast_entries[637], expected_1);
        assert_eq!(parsed_blast_entries[158], expected_2);
    }

    #[test]
    #[should_panic]
    fn test_parse_blast_2() {
        let test_file = "blast/results/CopZ-TcrZ_blast_result.txt";
        let test_file = PathBuf::from(test_file);
        let parsed_blast_entries = parse_annotations_from_blast_results(test_file);

        let expected_1 = BlastDerivedAnnotation {
            sseqid: "NZ_OD940440.2".to_string(),
            sstart: 261732,
            send: 261932,
            match_length: 67,
            sframe: 3,
            pident: 44.776,
            evalue: 8.68e-10,
            qseqid: "WP_002294571.1".to_string(),
            qstart: 1,
            qend: 67,
        };

        let expected_2 = BlastDerivedAnnotation {
            sseqid: "NZ_CP064343.1".to_string(),
            sstart: 1870968,
            send: 1870762,
            match_length: 69,
            sframe: -3,
            pident: 66.667,
            evalue: 8.45e-15,
            qseqid: "WP_010718490.1".to_string(),
            qstart: 1,
            qend: 69,
        };

        
        assert_eq!(parsed_blast_entries[637], expected_2);
        assert_eq!(parsed_blast_entries[158], expected_1);
    }

    use crate::cop_specific_analysis;
    #[test]
    fn test_assembly_names_in_table_match_directory_names() {
        
        // Define primary file paths
        let path_to_database = PathBuf::from("/run/media/katanga/SSD2/organized_prokaryotic_database/annotated");
        let refseq = PathBuf::from("database/data_assets/assembly_summary_refseq_updated.txt");
        let genbank = PathBuf::from("database/data_assets/assembly_summary_genbank.txt");
        let blast_data = cop_specific_analysis::build_list_of_blast_result_files("blast/results");
        
        // Load shared resources
        let (metadata, _, _) = cop_specific_analysis::load_resources(refseq, genbank, Some(blast_data));

        // Build task list from list of genome assembly directories
        let genome_dirs = cop_specific_analysis::build_list_of_genome_directories(path_to_database);

        for item in genome_dirs.iter() {
            let asm_name = item.file_name().and_then(|x| x.to_str()).map(|x| x.to_string()).unwrap();
            let attempt = metadata.get(&asm_name);
            attempt.unwrap();
        }
    }

}