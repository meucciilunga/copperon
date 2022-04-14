#![allow(dead_code)]

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::error::Error;

fn parse_assembly_report(assembly_report_path: PathBuf) -> Result<AssemblyReportMetadata, Box<dyn Error>> {
    let mut file = File::open(assembly_report_path)?;
    let mut contents = String::new();

    file.read_to_string(&mut contents)?;
    let metadata_content = contents.split('\n');

    // Pull target metadata entries
    let mut metadata_tmp: Vec<String> = Vec::with_capacity(6);
    let mut replicon_tmp: Vec<AssemblyReportRepliconMetadata> = Vec::new();
    let mut replicon_flag: bool = false;
    let metadata_tags = [
        "# Assembly name:  ",
        "# Organism name:  ",
        "# Taxid:          ",
        "# Assembly level: ",
        "# GenBank assembly accession: ",
        "# RefSeq assembly accession: ",
    ];

    // Interrogate file line-by-line
    for line in metadata_content {

        // ignore lines that are blank
        if line.len() == 0 {
            continue;
        }

        // Parse top-level metadata
        if !replicon_flag { // if replicon flag hasn't been triggered
            for tag in metadata_tags {
                if line.contains(tag) {
                    let mut new_line = line.split(tag).collect::<Vec<&str>>()[1].to_string();
                    new_line.pop().unwrap(); // remove newline character
                    metadata_tmp.push(new_line.to_string());
                    break;
                }
            }
        }

        // Anything after this line would correspond to individual replicon metadata
        if line.contains("# Sequence-Name	Sequence-Role",) {
            replicon_flag = true;
            continue;
        }

        // Parse individual replicon metadata
        if replicon_flag {
            let cleaned_line = line.strip_suffix('\n');
            let new_line = match cleaned_line {
                None => line,
                Some(x) => x,
            };

            let replicon_parse: Vec<&str> = new_line.split('\t').collect();
            
            let replicon_meta = AssemblyReportRepliconMetadata {
                seq_name: replicon_parse[0].to_string(),
                location: replicon_parse[3].to_string(),
                genbank_accession: replicon_parse[4].to_string(),
                refseq_accession: replicon_parse[6].to_string(),
                sequence_len: str::parse(replicon_parse[8]).unwrap(),
            };

            replicon_tmp.push(replicon_meta);
        }
    }

    // Consolidate parsed data
    let output = AssemblyReportMetadata {
        assembly_name: metadata_tmp[0].clone(),
        organism_name: metadata_tmp[1].clone(),
        taxid: str::parse(&metadata_tmp[2]).unwrap(),
        assembly_level: metadata_tmp[3].clone(),
        genbank_accesion: metadata_tmp[4].clone(),
        refseq_accesion: metadata_tmp[5].clone(),
        replicon_metadata: replicon_tmp,
    };

    Ok(output)
}

fn parse_genome_sequence() {}


#[derive(PartialEq, Debug)]
struct AssemblyReportMetadata {
    taxid: usize,
    assembly_name: String,
    organism_name: String,
    assembly_level: String,
    genbank_accesion: String,
    refseq_accesion: String,
    replicon_metadata: Vec<AssemblyReportRepliconMetadata>,
}

#[derive(PartialEq, Debug)]
struct AssemblyReportRepliconMetadata {
    seq_name: String,
    location: String,
    genbank_accession: String,
    refseq_accession: String,
    sequence_len: usize,
}

struct GeneralAnnotationData {}

struct CopYAnnoationData {}

struct GenomeSequenceData {}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use super::*;
    
    #[test]
    fn assembly_report_import_test_1() {
        
        // Import test asset
        let test_file = "test_assets/GCA_000152665.1_ASM15266v1/GCA_000152665.1_ASM15266v1_assembly_report.txt";
        let test_file = PathBuf::from(test_file);
        
        // Define correct answer manually
        let replicon_meta = AssemblyReportRepliconMetadata {
            seq_name: "ANONYMOUS".to_string(),
            location: "Chromosome".to_string(),
            genbank_accession: "CM000287.4".to_string(),
            refseq_accession: "NZ_CM000287.1".to_string(),
            sequence_len: 4109689_usize,
        };
        
        let expected_assembly_meta = AssemblyReportMetadata {
            taxid: 367459_usize,
            assembly_name: "ASM15266v1".to_string(),
            organism_name: "Clostridioides difficile QCD-32g58 (firmicutes)".to_string(),
            assembly_level: "Chromosome".to_string(),
            genbank_accesion: "GCA_000152665.1".to_string(),
            refseq_accesion: "GCF_000152665.1".to_string(),
            replicon_metadata: vec![replicon_meta],
        };

        let actual_assembly_meta = parse_assembly_report(test_file).unwrap();
        assert_eq!(expected_assembly_meta, actual_assembly_meta);
    }

    #[test]
    fn assembly_report_import_test_2() {
        
        // Import test asset
        let test_file = "test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_assembly_report.txt";
        let test_file = PathBuf::from(test_file);
        
        // Define correct answer manually
        let replicon_meta_1 = AssemblyReportRepliconMetadata {
            seq_name: "ANONYMOUS".to_string(),
            location: "Chromosome".to_string(),
            genbank_accession: "BA000022.2".to_string(),
            refseq_accession: "NC_000911.1".to_string(),
            sequence_len: 3573470_usize,
        };

        let replicon_meta_2 = AssemblyReportRepliconMetadata {
            seq_name: "pSYSA".to_string(),
            location: "Plasmid".to_string(),
            genbank_accession: "AP004311.1".to_string(),
            refseq_accession: "NC_005230.1".to_string(),
            sequence_len: 103307_usize,
        };

        let replicon_meta_3 = AssemblyReportRepliconMetadata {
            seq_name: "pSYSG".to_string(),
            location: "Plasmid".to_string(),
            genbank_accession: "AP004312.1".to_string(),
            refseq_accession: "NC_005231.1".to_string(),
            sequence_len: 44343_usize,
        };

        let replicon_meta_4 = AssemblyReportRepliconMetadata {
            seq_name: "pSYSM".to_string(),
            location: "Plasmid".to_string(),
            genbank_accession: "AP004310.1".to_string(),
            refseq_accession: "NC_005229.1".to_string(),
            sequence_len: 119895_usize,
        };

        let replicon_meta_5 = AssemblyReportRepliconMetadata {
            seq_name: "pSYSX".to_string(),
            location: "Plasmid".to_string(),
            genbank_accession: "AP006585.1".to_string(),
            refseq_accession: "NC_005232.1".to_string(),
            sequence_len: 106004_usize,
        };
        
        let expected_assembly_meta = AssemblyReportMetadata {
            taxid: 1148,
            assembly_name: "ASM972v1".to_string(),
            organism_name: "Synechocystis sp. PCC 6803 (cyanobacteria)".to_string(),
            assembly_level: "Complete Genome".to_string(),
            genbank_accesion: "GCA_000009725.1".to_string(),
            refseq_accesion: "GCF_000009725.1".to_string(),
            replicon_metadata: vec![
                replicon_meta_1,
                replicon_meta_2,
                replicon_meta_3,
                replicon_meta_4,
                replicon_meta_5,

            ],
        };

        let actual_assembly_meta = parse_assembly_report(test_file).unwrap();
        assert_eq!(expected_assembly_meta, actual_assembly_meta);
    }
}