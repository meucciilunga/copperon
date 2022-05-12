mod permutations;
mod import;
mod genome;
mod search;
mod analysis;
mod logging;
pub mod multithreader;


pub mod cop_specific_analysis {
    use std::io::{Write, BufWriter};
    use std::{collections::HashMap, sync::Arc, sync::Mutex, path::PathBuf};
    use crate::genome::BlastAssociationType;
    use crate::logging::SummaryLog;
    use crate::{analysis::*, import::AssemblyMetadata, permutations::SequencePermutations, genome::BlastHitsTable};
    use crate::multithreader::ThreadPoolTask;
    use std::fs;
    pub use crate::analysis::prepare_common_resources as load_resources;

    // Build set of all CopY operator permutations
    pub fn build_cop_permutation_table() -> (String, HashMap<char, Vec<char>>) {
        let operator_seq = "RNYKACANNYGTMRNY".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
    
            let character_flags = ['R', 'N', 'Y', 'K', 'M'];
            let substituents = [
                vec!['A', 'G'],
                vec!['A', 'C', 'G', 'T'],
                vec!['C', 'T'],
                vec!['G', 'T'],
                vec!['A', 'C'],
            ];
    
            for (key, val) in character_flags.into_iter().zip(substituents) {
                tmp_table.insert(key, val);
            }
            tmp_table
        };
    
        return(operator_seq, permutation_table)
    }

    // Build a closure for a worker thread to run
    pub fn build_genome_analysis_task(task_num: usize,
                                      genome_dir: PathBuf,
                                      metadata: Arc<HashMap<String, AssemblyMetadata>>,
                                      consensus: Arc<SequencePermutations>,
                                      blast: Arc<Option<Vec<BlastHitsTable>>>,
                                      genomes_summary_file: Arc<Mutex<BufWriter<fs::File>>>,
                                      orphans_file: Arc<Mutex<BufWriter<fs::File>>>,
                                      operators_file: Arc<Mutex<BufWriter<fs::File>>>,
                                      operons_file: Arc<Mutex<BufWriter<fs::File>>>
                                      ) -> ThreadPoolTask {
        
        let task = move || {
            println!("[{}] {}", task_num, genome_dir.file_name().unwrap().to_str().unwrap());
            let raw_genome = build_genome_from_dir(&genome_dir, &metadata);
            let processed_genome = search(&raw_genome, &consensus, &blast);
            
            processed_genome.log_table_entry(&genomes_summary_file);
            processed_genome.log_orphans(&orphans_file);
            processed_genome.log_operator_summaries(&operators_file);
            processed_genome.log_operon_summaries(&operons_file)
        };

        ThreadPoolTask::new(Box::new(task))
    }

    // Build paths to BLAST results data from a root directory
    pub fn build_list_of_blast_result_files(input_dir: &str) -> Vec<(BlastAssociationType, PathBuf)> {
        let dir_files = fs::read_dir(input_dir).unwrap();
        let mut output: Vec<(BlastAssociationType, PathBuf)> = Vec::new();
        
        for entry in dir_files {
            let path = entry.unwrap().path();
            let assoc_type = match path.file_name().and_then(|s| s.to_str()).unwrap() {
                "CopA" => BlastAssociationType::CopA,
                "CupA" => BlastAssociationType::CupA,
                "CopY" => BlastAssociationType::CopY,
                "CopZ" => BlastAssociationType::CopZ,
                  _ => panic!("INVALID BLAST DATA."),
            };

            let blast_datum = (assoc_type, path);
            output.push(blast_datum);
        }

        output
    }

    // Build paths to genome assembly directories from a root directory
    pub fn build_list_of_genome_directories(genomes_root_dir: &PathBuf) -> Vec<PathBuf> {
        let dir_files = fs::read_dir(genomes_root_dir).unwrap();
        dir_files.into_iter().map(|x| x.unwrap().path()).collect::<Vec<PathBuf>>()
    }

    // Build a shared log file for storing analysis summary data derived from multiple threads
    pub fn generate_shared_logging_file(log_file_path: &PathBuf) -> Arc<Mutex<BufWriter<fs::File>>> {

        let file = fs::OpenOptions::new().append(true)
                                         .create_new(true)
                                         .open(log_file_path)
                                         .expect("ERROR: could not create shared log file!");
        let buffered_file = BufWriter::new(file);
        let guarded_file = Mutex::new(buffered_file);
        
        Arc::new(guarded_file)
    }
}