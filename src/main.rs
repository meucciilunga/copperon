use num_cpus;
use std::io::Write;
use std::sync::Arc;
use std::path::PathBuf;
use copperon::multithreader;
use copperon::cop_specific_analysis;
use std::time::Instant;

fn main() {

    // Start threadpool
    let num_cpus: usize = num_cpus::get();
    let threadpool = multithreader::Manager::new(2*num_cpus);

    // Define primary file paths
    let path_to_database = PathBuf::from("/run/media/katanga/SSD2/organized_prokaryotic_database/annotated");
    let refseq = PathBuf::from("database/data_assets/assembly_summary_refseq_updated.txt");
    let genbank = PathBuf::from("database/data_assets/assembly_summary_genbank.txt");
    let blast_data = cop_specific_analysis::build_list_of_blast_result_files("blast/results");

    // Generate shared log files for storage
    let orphans_file = cop_specific_analysis::generate_shared_logging_file(&PathBuf::from("output/orphans.log"));
    let genomes_summary_file = cop_specific_analysis::generate_shared_logging_file(&PathBuf::from("output/genomes.summary"));
    let operators_summary_file = cop_specific_analysis::generate_shared_logging_file(&PathBuf::from("output/operators.summary"));
    let operon_summary_file = cop_specific_analysis::generate_shared_logging_file(&PathBuf::from("output/operons.summary"));
    
    // Load shared resources
    let (metadata, blast, consensus) = cop_specific_analysis::load_resources(refseq, genbank, Some(blast_data));

    // Load paths to all genome directories
    let genome_dirs = cop_specific_analysis::build_list_of_genome_directories(&path_to_database);
    
    // Build task list from list of genome assembly directories
    let mut tasks: Vec<multithreader::ThreadPoolTask> = Vec::new();
    for (index, dir) in genome_dirs.iter().enumerate() {
        let metadata = Arc::clone(&metadata);
        let consensus = Arc::clone(&consensus);
        let blast = Arc::clone(&blast);
        let genomes_log = Arc::clone(&genomes_summary_file);
        let orphans_log = Arc::clone(&orphans_file);
        let operators_log = Arc::clone(&operators_summary_file);
        let operons_log = Arc::clone(&operon_summary_file);

        let new_task = cop_specific_analysis::build_genome_analysis_task(index+1, dir.clone(), metadata, consensus, blast, genomes_log, orphans_log, operators_log, operons_log);
        tasks.push(new_task);
    }

    let now = Instant::now();
    let num_to_analyze = tasks.len();

    // Run analysis
    threadpool.execute_queue(tasks);
    let total_time = now.elapsed().as_millis();

    // Flush logs
    orphans_file.lock().expect("ERROR: COULD NOT ACQUIRE LOCK TO FLUSH ORPHANS").flush().expect("ERROR: COULD NOT FLUSH ORPHANS LOG");
    genomes_summary_file.lock().expect("ERROR: COULD NOT ACQUIRE LOCK TO FLUSH GENOMES").flush().expect("ERROR: COULD NOT FLUSH GENOMES LOG");

    
    let average_time = total_time / num_to_analyze as u128;
    println!("Total Analysis Time: {}ms", total_time / num_to_analyze as u128);
    println!("Average Time: {}ms", average_time);
}