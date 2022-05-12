use crate::search::{OperatorType, SearchGenome, SearchGene, Operator, Operon};
use crate::genome::{self, RepliconType, BlastAssociationType};
use std::fs;
use std::io::{BufWriter, Write};
use std::sync::Mutex;
use std::fmt::Display;

pub trait SummaryLog {
    fn log_table_entry(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>);
}

#[allow(non_snake_case)]
impl<'genome, 'blast, 'query> SummaryLog for SearchGenome<'genome, 'blast, 'query> {

    fn log_table_entry(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>) {

        let operators_expectation = match self.operator_query {
            Some(x) => x.expectation,
            None => 0.0,
        };

        let actual_num_operators = match &self.operators {
            None => 0,
            Some(x) => x.len(),
        };

        let num_interior_operators = match &self.operators {
            None => 0_usize,
            Some(x) => x.iter().filter(|x| x.operator_type == OperatorType::Interior).count()
        }; 

        let num_exterior_operators = match &self.operators {
            None => 0_usize,
            Some(x) => x.iter().filter(|x| x.operator_type == OperatorType::Exterior).count()
        };

        let num_boundary_operators = match &self.operators {
            None => 0_usize,
            Some(x) => x.iter().filter(|x| x.operator_type == OperatorType::Boundary).count()
        }; 

        let (num_operons, num_operon_genes) = match &self.operons {
            None => (0,0),
            Some(x) => {
                let a = x.len();
                let b = x.iter().fold(0_usize, |acc, x| acc + x.genes.len());

                (a,b)
            }
        };

        let num_replicons = self.genome.replicons.len();

        let (gap_len, genes_len, total_genome_len) = match self.total_intergenic_space() {
            Some(x) => x,
            None => (0_usize, 0_usize, 0_usize),
        };

        let actual_operator_density_per_1Mbp = actual_num_operators as f64 / (total_genome_len as f64 / 1_000_000.0);

        let num_orphans = match self.CopY_orphans() {
            Some(x) => x.len(),
            None => 0,
        };

        let expected_num_operators = operators_expectation * total_genome_len as f64;
        let expected_operator_density_per_1Mbp = operators_expectation * 1_000_000.0;

        let gap_frac = (gap_len as f64) / (total_genome_len as f64);
        let gene_frac = (genes_len as f64) / (total_genome_len as f64);

        let assoc_types = [
            (genome::BlastAssociationType::CopA, 'A'),
            (genome::BlastAssociationType::CupA, 'a'),
            (genome::BlastAssociationType::CopY, 'Y'),
            (genome::BlastAssociationType::CopZ, 'Z'),
        ];

        let assoc_code = match &self.blast_associations() {
            Some(table) => {
                let mut code = "[".to_string();
                for (bat, c) in assoc_types.iter() {
                    if table.contains(bat) {
                        code.push(*c);
                    }
                }
                code.push(']');

                code
            }
            
            None => String::from(""),
        };

        let output = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                                self.genome.metadata.assembly_accession,
                                self.genome.metadata.species_tax_id,
                                self.genome.metadata.tax_id,
                                self.genome.metadata.organism_name,
                                self.genome.metadata.infraspecifc_name,
                                num_replicons,
                                total_genome_len,
                                gap_frac,
                                gene_frac,
                                expected_num_operators,
                                actual_num_operators,
                                expected_operator_density_per_1Mbp,
                                actual_operator_density_per_1Mbp,
                                num_exterior_operators,
                                num_interior_operators,
                                num_boundary_operators,
                                num_operons,
                                num_operon_genes,
                                num_orphans,
                                assoc_code);

        // Append log entry to file
        let mut file = shared_log_file_stream.lock().expect("ERROR: COULD NOT ACQUIRE MUTEX LOCK FOR WRITING TO SHARED LOG FILE.");
        file.write_all(output.as_bytes()).expect("ERROR: COULD NOT WRITE LOG ENTRY TO BUFFER.");
    }
}

impl<'genome, 'blast, 'query> SearchGenome<'genome, 'blast, 'query> {

    // Log a genome analysis's list of CopY orphans
    pub fn log_orphans(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>) {
        
        let orphans = match self.CopY_orphans() {
            Some(x) => x,
            None => return,
        };

        let mut output = "".to_string();
        for gene in orphans {
            let tmp = format!("{}\n", gene);
            output.push_str(&tmp);
        }

        let mut file = shared_log_file_stream.lock().expect("ERROR: could not lock 'CopY orphans' log file!");
        file.write_all(output.as_bytes()).expect("ERROR: COULD NOT WRITE 'ORPHANS' ENTRY TO BUFFER.");
    }

    pub fn log_operon_summaries(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>) {
        if let Some(x) = &self.operons {
            for operon in x {
                operon.log_table_entry(shared_log_file_stream)
            }
        }
    }

    pub fn log_operator_summaries(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>) {
        if let Some(x) = &self.operators {
            for operator in x {
                operator.log_table_entry(shared_log_file_stream)
            }
        }
    }
}

impl<'genome> SummaryLog for Operator<'genome> {

    fn log_table_entry(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>) {
        let operator_seq = &self.seq;
        let replicon = &self.linear_location.replicon.accession_id;
        let replicon_type = &self.linear_location.replicon.replicon_type;

        let num_operon_genes = match &self.operon {
            Some(x) => x.genes.len(),
            None => 0,
        };
        
        let postional_status = &self.operator_type;
        let class = &self.operator_class;
        let palindrome = &self.dimension;

        let assoc_types = [
            (BlastAssociationType::CopA, 'A'),
            (BlastAssociationType::CupA, 'a'),
            (BlastAssociationType::CopY, 'Y'),
            (BlastAssociationType::CopZ, 'Z'),
        ];

        let assoc_code = match &self.blast_association {
            Some(table) => {
                let mut code = "[".to_string();
                for (bat, c) in assoc_types.iter() {
                    if table.contains(bat) {
                        code.push(*c);
                    }
                }
                code.push(']');

                code
            }
            
            None => String::from(""),
        };

        let output = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                                operator_seq,
                                replicon,
                                replicon_type,
                                self.linear_location.start_bound+1,
                                self.linear_location.strand,
                                num_operon_genes,
                                postional_status,
                                class,
                                palindrome,
                                assoc_code);

        let mut file = shared_log_file_stream.lock().expect("ERROR: could not lock 'Operator Summary' log file!");
        file.write_all(output.as_bytes()).expect("ERROR: COULD NOT WRITE 'OPERATOR SUMMARY' ENTRY TO BUFFER.");
    }
}



impl<'genome> SummaryLog for Operon<'genome> {

    fn log_table_entry(&self, shared_log_file_stream: &Mutex<BufWriter<fs::File>>) {
        let replicon = &self.linear_location.replicon.accession_id;
        let replicon_type = &self.linear_location.replicon.replicon_type;
        
        let num_operators = match &self.operators {
            Some(x) => x.len(),
            None => 0,
        };

        let num_genes = self.genes.len();

        let output = format!("{}\t{}\t{}\t{}\t{}\t{}\n", 
                                replicon,
                                replicon_type,
                                self.linear_location.start_bound+1,
                                self.linear_location.strand,
                                num_operators,
                                num_genes);

        let mut file = shared_log_file_stream.lock().expect("ERROR: could not lock 'Operator Summary' log file!");
        file.write_all(output.as_bytes()).expect("ERROR: COULD NOT WRITE 'OPERATOR SUMMARY' ENTRY TO BUFFER.");
    }
}

impl<'genome> Display for SearchGene<'genome> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let name = match self.gene.attributes.get("Name") {
            None => match self.gene.attributes.get("ID") {
                Some(x) => x,
                None => "NO NAME",
            },
            Some(x) => x,
        };

        let pdt = match self.gene.attributes.get("product") {
            None => match self.gene.attributes.get("ID") {
                Some(x) => x,
                None => "UNKNOWN PRODUCT",
            },
            Some(x) => x,
        };

        write!(f, "{}\t{}\t{}\t{}\t{}\n", name, pdt, 
                    self.linear_location.replicon.accession_id, 
                    self.linear_location.start_bound+1, 
                    self.linear_location.end_bound)
    }
}

impl Display for RepliconType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match self {
            RepliconType::Chromosome => "Chromosome",
            RepliconType::Plasmid => "Plasmid",
        };

        write!(f, "{}", text)
    }
}