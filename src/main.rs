use num_cpus;
use std::sync::Arc;
use std::path::PathBuf;
use copperon::multithreader;
use copperon::cop_specific_analysis;
use std::time::Instant;

fn main() {

    // Start threadpool
    let num_cpus: usize = num_cpus::get();
    let threadpool = multithreader::Manager::new(2 * num_cpus);

    // Define primary file paths
    let path_to_database = PathBuf::from("/run/media/katanga/SSD2/organized_prokaryotic_database/annotated");
    let refseq = PathBuf::from("database/data_assets/assembly_summary_refseq_updated.txt");
    let genbank = PathBuf::from("database/data_assets/assembly_summary_genbank.txt");
    let blast_data = cop_specific_analysis::build_list_of_blast_result_files("blast/results");
    
    // Load shared resources
    let (metadata, blast, consensus) = cop_specific_analysis::load_resources(refseq, genbank, Some(blast_data));

    // Load paths to all genome directories
    let genome_dirs = cop_specific_analysis::build_list_of_genome_directories(path_to_database);
    
    // Build task list from list of genome assembly directories
    let mut tasks: Vec<multithreader::ThreadPoolTask> = Vec::new();
    for (index, dir) in genome_dirs.into_iter().enumerate() {
        let metadata = Arc::clone(&metadata);
        let consensus = Arc::clone(&consensus);
        let blast = Arc::clone(&blast);
        let new_task = cop_specific_analysis::build_genome_analysis_task(index+1, dir, metadata, consensus, blast);

        tasks.push(new_task);
    }

    let now = Instant::now();
    let num_to_analyze = tasks.len();
    // Run analysis
    threadpool.execute_queue(tasks);

    let total_time = now.elapsed().as_millis();
    let average_time = total_time / num_to_analyze as u128;
    println!("Total Analysis Time: {}s", total_time / 1000);
    println!("Average Time: {}s", average_time);
}