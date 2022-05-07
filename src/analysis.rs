use std::fs;
use std::sync::Arc;
use std::hash::Hash;
use std::path::PathBuf;
use std::collections::HashMap;
use crate::genome::{self, BlastHitsTable, StrandSense};
use crate::import::{self, AssemblyMetadata};
use crate::permutations::{self, SequencePermutations};
use crate::search::{self, Operator, NeighborType, OverlapType, Relationship,
                    SearchGene, SearchGenome, SpatialRelationship};
use crate::cop_specific_analysis::build_cop_permutation_table;

#[derive(PartialEq, Eq, Hash)]
enum GenomeFileType {
    Annotation,
    FASTA,
}

#[derive(PartialEq, Eq, Hash, Debug)]
pub enum OperatorType {
    Exterior,
    Interior,
}

// Load datasets that are to be shared across threads
pub fn prepare_common_resources(refseq: PathBuf, genbank: PathBuf, path_to_blasts: Option<Vec<(genome::BlastAssociationType, PathBuf)>>) 
-> (Arc<HashMap<String, AssemblyMetadata>>, Arc<Option<Vec<BlastHitsTable>>>, Arc<SequencePermutations>) {

    // Import respective RefSeq/GenBank Database Metadata Files
    let refseq_meta = import::import_database_summary(refseq).into_iter();
    let genbank_meta = import::import_database_summary(genbank).into_iter();
    let full_metadata: HashMap<String, import::AssemblyMetadata> = refseq_meta.chain(genbank_meta).collect();
    let metadata = Arc::new(full_metadata);

    // Load BLAST results into a replicon-searchable table
    let blast_tables = match path_to_blasts {

        Some(tables) => {
            let mut storage: Vec<genome::BlastHitsTable> = Vec::new();

            for (id, path) in tables {
                let new_table = genome::BlastHitsTable::build_table_from_file(id, path);
                storage.push(new_table);
            }

            Some(storage)
        },

        None => None
    };
    let blast = Arc::new(blast_tables);

    // Load CopY Operators
    let (operator_seq, table) = build_cop_permutation_table();
    let operators = permutations::SequencePermutations::new("CopY-family Consensus Operator".to_string(), operator_seq, table);
    let operators = Arc::new(operators);

    (metadata, blast, operators)
}

// Build Genome from a directory containing its core files
pub fn build_genome_from_dir(dir: PathBuf, database_metadata: &HashMap<String, AssemblyMetadata>) -> genome::Genome {

    // Define name of assembly by its directory name
    let assembly_name = dir.file_name().and_then(|s| s.to_str()).and_then(|s| Some(s.to_string())).unwrap();

    // Hash paths by extension
    let mut genome_files: HashMap<GenomeFileType, PathBuf> = HashMap::new();

    // Open directory and sort files inside
    if dir.is_dir() {
        let dir_files = fs::read_dir(dir).unwrap();

        for entry in dir_files {
            let path = entry.unwrap().path();
            let extension = path.extension()
                                .expect("ERROR: file in directory is missing an extesion.")
                                .to_str()
                                .unwrap();

            match extension {
                "fna" => {
                    let key = GenomeFileType::FASTA;
                    let val = path;
                    genome_files.insert(key, val);
                },
                "gff" => {
                    let key =  GenomeFileType::Annotation;
                    let val = path;
                    genome_files.insert(key, val);
                },
                    _  => {}
            }
        }

    } else {
        panic!("PATH IS NOT A VALID DIRECTORY.");
    }

    // Pull file paths from the hash table; pull assembly name from name of the directory
    let genome_seq_file = genome_files.get(&GenomeFileType::FASTA).unwrap();
    let genome_annotation_file = genome_files.get(&GenomeFileType::Annotation);

    // Import data relevant to TIGR4
    let protogenome = import::parse_genome_sequence(&assembly_name, genome_seq_file);
    let genes = match genome_annotation_file {
        Some(annotation_file) => {
            let protogenes = import::parse_genome_annotation(&annotation_file);
            Some(genome::Gene::convert_annotation_entry_list(&protogenes))
        }
        None => None,
    };

    // Build genome
    let asm_pull_error = "ERROR: could not find assembly name in metadata database!";
    let metadata = database_metadata.get(&assembly_name).expect(asm_pull_error);
    
    genome::Genome::from_proto(protogenome, metadata, genes)
}

pub fn search<'genome, 'blast>(target: &'genome genome::Genome, 
                                    consensus: &SequencePermutations, 
                                    tables: &'blast Option<Vec<BlastHitsTable>>) -> search::SearchGenome<'genome, 'blast> {
    search::SearchGenome::new(target, consensus, tables)
}

pub fn predict_operator_genes<'genome>(input: &'genome Operator, input_genome: &'genome SearchGenome) -> Option<Operon<'genome>> {
    
    // Can't search genome for annotated genes if it has no gene annotation data
    let neighbors = match &input_genome.genes {
        Some(list_of_genes) => list_of_genes.clone(),
        None => return None,
    };

    // Filter off genes (which are generally duplicates of CDS) and the start region
    let gene_region_filter = |x: &SearchGene| -> bool {
        (x.gene.feature_type != genome::FeatureType::Region) && (x.gene.feature_type != genome::FeatureType::Gene)
    };

    let neighbors = neighbors.into_iter()
                                 .filter(|x| gene_region_filter(x))
                                 .collect::<Vec<SearchGene>>();

    // Ensure each feature region has only one representative annotation (i.e., if a pseudogene and CDS represent
    // the same genome region in the annotation table, only one of them should be included in the neighbors list)
    let mut unique_neighbors: HashMap<(usize, usize), SearchGene> = HashMap::with_capacity(neighbors.len());
    for neighbor in neighbors.into_iter() {
        let indicies = (neighbor.linear_location.start_bound, neighbor.linear_location.end_bound);
        unique_neighbors.insert(indicies, neighbor);
    }
    let mut neighbors = unique_neighbors.into_values().collect::<Vec<SearchGene>>();

    // Sort by finding a distance of Element A to input gene
    let sort_routine = |a: &search::SearchGene, b: &search::SearchGene| {
        
        // Find distance of a to input gene
        let a_dist = match a.relative_to(input) {
            search::SpatialRelationship::Overlap(_) => 0,
            search::SpatialRelationship::None => usize::MAX,
            search::SpatialRelationship::Neighbor(variant) => {
                match variant {
                    NeighborType::FivePrime(x) => x,
                    NeighborType::ThreePrime(x) => x,
                }
            }
        };

        // Find distance of b to input gene
        let b_dist = match b.relative_to(input) {
            search::SpatialRelationship::Overlap(_) => 0,
            search::SpatialRelationship::None => usize::MAX,
            search::SpatialRelationship::Neighbor(variant) => {
                match variant {
                    NeighborType::FivePrime(x) => x,
                    NeighborType::ThreePrime(x) => x,
                }
            }
        };

        a_dist.cmp(&b_dist)
    };

    // Sort neighbors to input by distance from left or right bound of feature (whichever is less)
    neighbors.sort_by(sort_routine);

    // Filtering Functions
    let five_prime_filter = |relationship| -> bool {
        match relationship {
            SpatialRelationship::Neighbor(NeighborType::FivePrime(_)) => true,
            _ => false,
        }
    };

    let three_prime_filter = |relationship| -> bool {
        match relationship {
            SpatialRelationship::Neighbor(NeighborType::ThreePrime(_)) => true,
            _ => false,
        }
    };

    let overlap_filter = |relationship| -> bool {
        match relationship {
            SpatialRelationship::Overlap(OverlapType::ContainerOf) => true,
            SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary) => true,
            SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary) => true,
            _ => false,
        }
    };

    const MAX_NEIGHBOR_DIST: usize = 40_000;
    let overlap_neighbors = neighbors.iter()
                                     .filter(|&g| overlap_filter(g.relative_to(input)))
                                     .map(|x| x.clone())
                                     .collect::<Vec<SearchGene>>();

    let mut five_prime_neighbors = neighbors.iter()
                                            .filter(|&g| five_prime_filter(g.relative_to(input)))
                                            .filter(|&g| g.relative_to(input).distance() <= MAX_NEIGHBOR_DIST)
                                            .map(|x| x.clone())
                                            .collect::<Vec<SearchGene>>();

    let mut three_prime_neighbors = neighbors.iter()
                                             .filter(|&g| three_prime_filter(g.relative_to(input)))
                                             .filter(|&g| g.relative_to(input).distance() <= MAX_NEIGHBOR_DIST)
                                             .map(|x| x.clone())
                                             .collect::<Vec<SearchGene>>(); 

    // Neighbors who are overlapping the input always go to the top of the neighbors list
    // in whose direction they match the first element of; if the overlapping neighbor matches
    // neither of the surrounding directions, we know that the overlapping neighbor (and 
    // its corresponding direction) must define the entire 'pseudo-operon'
    let five_prime_init_strand = match five_prime_neighbors.get(0) {
        Some(x) => x.linear_location.strand.clone(),
        None => StrandSense::Other,
    };

    let three_prime_init_strand = match three_prime_neighbors.get(0) {
        Some(x) => x.linear_location.strand.clone(),
        None => StrandSense::Other,
    };
    
    let mut engulfing_neighbors: Vec<SearchGene<'genome>> = vec![];

    // Determine operator type now before consuming the overlap_neighbors vector
    let operator_type = if overlap_neighbors.len() == 0 {
        OperatorType::Exterior
    } else {
        OperatorType::Interior
    };

    for neighbor in overlap_neighbors.into_iter().rev() {

        let five_prime_cond = neighbor.linear_location.strand == five_prime_init_strand;
        let three_prime_cond = neighbor.linear_location.strand == three_prime_init_strand;

        if five_prime_cond && three_prime_cond {
            three_prime_neighbors.insert(0, neighbor.clone());
            five_prime_neighbors.insert(0, neighbor);
        } else if five_prime_cond {
            five_prime_neighbors.insert(0, neighbor);
        } else if three_prime_cond {
            three_prime_neighbors.insert(0, neighbor);
        } else {
            engulfing_neighbors.push(neighbor);
        }
    }

    // Assemble proto-operons, which is the list of contiguous genes on either 5' or 3' reasonably spaced away from operator
    const MAX_INITIAL_OPERON_DIST: usize = 100;
    const MAX_INTER_OPERON_DIST: usize = 25;

    let proto_operon_builder = |list_of_neighbors_by_dist: Vec<SearchGene<'genome>>| -> Vec<SearchGene<'genome>> {

        // If input list is empty, return empty output list
        if list_of_neighbors_by_dist.len() == 0 {
            let tmp: Vec<SearchGene<'genome>> = Vec::new();
            return tmp;
        }

        let mut tmp_storage: Vec<(SearchGene, usize)> = Vec::new();

        let init_dist = list_of_neighbors_by_dist[0].relative_to(input).distance();
        let init_strand = list_of_neighbors_by_dist[0].linear_location.strand.clone();
        if init_dist > MAX_INITIAL_OPERON_DIST {
            let tmp: Vec<SearchGene<'genome>> = Vec::new();
            return tmp;
        }
        tmp_storage.push((list_of_neighbors_by_dist[0].clone(), init_dist));

        let mut offset_neighbors: Vec<SearchGene> = Vec::new();
        offset_neighbors.extend_from_slice(&list_of_neighbors_by_dist[1..]);

        for (first, second) in list_of_neighbors_by_dist.into_iter().zip(offset_neighbors.into_iter()) {
            let dist = first.relative_to(&second).distance();
            let new_neighbor = (second, dist);
            tmp_storage.push(new_neighbor);
        }

        // Because of init_dist test above, if the code reaches this far, the first neighbor is confirmed under the operon;
        // this line adds it to the confirmed operator list, while removing it from the tmp_storage list
        let mut output: Vec<SearchGene> = vec![tmp_storage.remove(0).0];

        // Must validate every neighbor after the first
        for (neighbor, dist) in tmp_storage {

            let break_condition = (dist > MAX_INTER_OPERON_DIST) || (neighbor.linear_location.strand != init_strand);

            if break_condition {
                break;
            }

            output.push(neighbor);
        }

        output
    };


    let mut storage_vec: Vec<SearchGene> = Vec::new();
    let mut five_prime_genes = proto_operon_builder(five_prime_neighbors);
    let mut three_prime_genes = proto_operon_builder(three_prime_neighbors);
    let mut engulfing_genes = proto_operon_builder(engulfing_neighbors);

    // If operator is on the 3'-end of a block of neighbors, operator is likely not
    // a valid regulator of that block since it occurs at the end of their run; this
    // check is obviously only needed if the operator has any neighbors within range
    let end_of_five_prime_neighborhood_block = if five_prime_genes.len() != 0 {
        match input.relative_to(&five_prime_genes[0]) {
            SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary) => true,
            SpatialRelationship::Neighbor(NeighborType::ThreePrime(_)) => true,
            _ => false,
        }
    } else {
        false
    };

    let end_of_three_prime_neighborhood_block = if three_prime_genes.len() != 0 {
        match input.relative_to(&three_prime_genes[0]) {
        SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary) => true,
        SpatialRelationship::Neighbor(NeighborType::ThreePrime(_)) => true,
        _ => false,
        }
    } else {
        false
    };

    let five_prime_invalid = end_of_five_prime_neighborhood_block && (three_prime_genes.len() == 0);
    let three_prime_invalid = end_of_three_prime_neighborhood_block && (five_prime_genes.len() == 0);
    let bidirectional_invalid = end_of_five_prime_neighborhood_block && end_of_three_prime_neighborhood_block;

    // Engulfing genes are always valid operon contenders; genes to which the operator is a three prime
    // neighbor of any variety are not very likely to be valid operon contenders, so are excluded
    storage_vec.append(&mut engulfing_genes);

    if !(five_prime_invalid || bidirectional_invalid) {
        storage_vec.append(&mut five_prime_genes);
    }

    if !(three_prime_invalid || bidirectional_invalid) {
        storage_vec.append(&mut three_prime_genes);
    }
    
    // Similar to above; remove duplicates
    let mut unique_operon_elements: HashMap<(usize, usize), SearchGene> = HashMap::with_capacity(storage_vec.len());
    for neighbor in storage_vec.into_iter() {
        let indicies = (neighbor.linear_location.start_bound, neighbor.linear_location.end_bound);
        unique_operon_elements.insert(indicies, neighbor);
    }
    let regulated_genes = unique_operon_elements.into_values().collect::<Vec<SearchGene>>();


    match regulated_genes.len() {
        0 => None,
        _ => {
            let output = Operon {
                operator: input,
                operon: regulated_genes,
                operator_type,
            };

            Some(output)
        }
    }
}

// Simple function for proposing an operon for an operator based on clusters of nearby genes
pub fn find_operons<'genome>(input: &'genome SearchGenome) -> Option<Vec<Operon<'genome>>> {

    let operators = match &input.operators {
        None => return None,
        Some(x) => x,
    };

    let mut storage: Vec<Operon<'genome>> = Vec::new();
    for operator in operators {
        match predict_operator_genes(operator, input) {
            None => continue,
            Some(x) => storage.push(x),
        }
    }

    match storage.len() {
        0 => None,
        _ => Some(storage),
    }
}

pub struct Operon<'genome> {
    pub operator: &'genome Operator<'genome>,
    pub operon: Vec<SearchGene<'genome>>,
    pub operator_type: OperatorType,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_operons() {
        let refseq = PathBuf::from("tests/test_assets/assembly_summary_refseq.txt");
        let genbank = PathBuf::from("tests/test_assets/assembly_summary_genbank.txt");
        let t4_dir = PathBuf::from("tests/test_assets/GCF_000006885.1_ASM688v1");
        let (metadata, blast, consensus) = prepare_common_resources(refseq, genbank, None);
        let t4_genome = build_genome_from_dir(t4_dir, &metadata);
        let t4_search = search(&t4_genome, &consensus, &blast);

        match find_operons(&t4_search) {
            None => {},
            Some(x) => {
                assert_eq!(14, x.len());
                
                println!("T4 has {} Operons:\n", x.len());
                for operon in x {
                    println!("{}", operon.operator);
                    println!("{:?}", operon.operator_type);
                    println!("{}", operon.operon.len());
                    for gene in operon.operon {
                        println!("{:?}", gene.gene.feature_type);
                    }
                    println!("\n\n\n");
                }
            }
        }
    }

}