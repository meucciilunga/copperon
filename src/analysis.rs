use std::fs;
use std::sync::Arc;
use std::hash::Hash;
use std::path::PathBuf;
use std::collections::{HashMap, HashSet};
use crate::genome::{self, BlastHitsTable, StrandSense, BlastAssociationType};
use crate::import::{self, AssemblyMetadata};
use crate::permutations::{self, SequencePermutations};
use crate::search::{self, Operator, NeighborType, OverlapType, Relationship,
                    SearchGene, SearchGenome, SpatialRelationship, Operon, OperatorType,
                    LinearGenomeLocation, LocateOnLinearGenome, CircularGenomeLocation};
use crate::cop_specific_analysis;

#[derive(PartialEq, Eq, Hash)]
enum GenomeFileType {
    Annotation,
    FASTA,
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
    let (consensus_operator_seq, table) = cop_specific_analysis::build_cop_permutation_table();
    let operators = permutations::SequencePermutations::new(consensus_operator_seq, table);
    let operators = Arc::new(operators);

    (metadata, blast, operators)
}

// Build Genome from a directory containing its core files
pub fn build_genome_from_dir(dir: &PathBuf, database_metadata: &HashMap<String, AssemblyMetadata>) -> genome::Genome {

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

pub fn search<'genome, 'blast, 'query>( target:    &'genome genome::Genome, 
                                        consensus: &'query SequencePermutations, 
                                        tables:    &'blast Option<Vec<BlastHitsTable>>) -> search::SearchGenome<'genome, 'blast, 'query> {

    // Searches for operator sequences
    let mut output = search::SearchGenome::new(target, consensus, tables);
    
    // Locates operons and updates information about operators based on their operons
    output.locate_operons();
    
    output
}

fn predict_operon<'genome, 'blast>(input: &Operator<'genome>, input_genes: &Option<Vec<SearchGene<'genome>>>) -> Option<Vec<SearchGene<'genome>>> {
    
    // Can't search genome for annotated genes if it has no gene annotation data
    let neighbors = match input_genes {
        Some(list_of_genes) => list_of_genes.clone(),
        None => return None,
    };

    // Filter off 'Gene' annotations (which are generally duplicates of CDS) and the start 'Regions'
    let gene_region_filter = |x: &SearchGene| -> bool {
        (x.gene.feature_type != genome::FeatureType::Region) && (x.gene.feature_type != genome::FeatureType::Gene)
    };

    let neighbors = neighbors.into_iter()
                             .filter(|x| gene_region_filter(x))
                             .collect::<Vec<SearchGene>>();

    // Ensure each feature region has only one representative annotation (i.e., if a pseudogene and CDS represent
    // the same genome region in the annotation table, only one of them should be included in the neighbors list)
    let mut unique_neighbors: HashMap<(usize, usize, String), SearchGene> = HashMap::with_capacity(neighbors.len());
    for neighbor in neighbors.into_iter() {
        let loc = (neighbor.linear_location.start_bound, 
                   neighbor.linear_location.end_bound, 
                   neighbor.linear_location.replicon.accession_id.clone());
        unique_neighbors.insert(loc, neighbor);
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
            SpatialRelationship::Overlap(_) => true,
            _ => false,
        }
    };

    const MAX_NEIGHBOR_DIST: usize = 40_000;
    let overlap_neighbors = neighbors.iter()
                                     .filter(|&g| overlap_filter(g.relative_to(input)))
                                     .cloned()
                                     .collect::<Vec<SearchGene>>();

    let mut five_prime_neighbors = neighbors.iter()
                                            .filter(|&g| five_prime_filter(g.relative_to(input)))
                                            .filter(|&g| g.relative_to(input).distance() <= MAX_NEIGHBOR_DIST)
                                            .cloned()
                                            .collect::<Vec<SearchGene>>();

    let mut three_prime_neighbors = neighbors.iter()
                                             .filter(|&g| three_prime_filter(g.relative_to(input)))
                                             .filter(|&g| g.relative_to(input).distance() <= MAX_NEIGHBOR_DIST)
                                             .cloned()
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
    const MAX_INTER_OPERON_DIST: usize = 40;

    let proto_operon_builder = |list_of_neighbors_by_smallest_dist: Vec<SearchGene<'genome>>| -> Vec<SearchGene<'genome>> {

        // If input list is empty, return empty output list
        if list_of_neighbors_by_smallest_dist.len() == 0 {
            let tmp: Vec<SearchGene<'genome>> = Vec::new();
            return tmp;
        }

        let mut tmp_storage: Vec<(SearchGene, usize)> = Vec::new();

        let init_dist = list_of_neighbors_by_smallest_dist[0].relative_to(input).distance();
        let init_strand = list_of_neighbors_by_smallest_dist[0].linear_location.strand.clone();

        // If closest neighbor to operator is too far away, return empty list   
        if init_dist > MAX_INITIAL_OPERON_DIST {
            let tmp: Vec<SearchGene<'genome>> = Vec::new();
            return tmp;
        }
        tmp_storage.push((list_of_neighbors_by_smallest_dist[0].clone(), init_dist));

        let mut offset_neighbors: Vec<SearchGene> = Vec::new();
        offset_neighbors.extend_from_slice(&list_of_neighbors_by_smallest_dist[1..]);

        for (first, second) in list_of_neighbors_by_smallest_dist.into_iter().zip(offset_neighbors.into_iter()) {
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

    let neighborhood_gap_too_big = if (five_prime_genes.len() != 0) && (three_prime_genes.len() != 0) {
        if five_prime_genes[0].relative_to(&three_prime_genes[0]).distance() > MAX_INTER_OPERON_DIST {
            true
        } else {
            false
        }
    } else {
        true
    };

    let five_prime_invalid = end_of_five_prime_neighborhood_block && (three_prime_genes.len() == 0 || neighborhood_gap_too_big);
    let three_prime_invalid = end_of_three_prime_neighborhood_block && (five_prime_genes.len() == 0 || neighborhood_gap_too_big);
    let bidirectional_invalid = end_of_five_prime_neighborhood_block && end_of_three_prime_neighborhood_block;

    // Engulfing genes are always valid operon contenders; genes to which the operator is a three prime
    // neighbor are not very likely to be valid operon contenders (unless occuring in larger operon), so are excluded
    storage_vec.append(&mut engulfing_genes);

    if !(five_prime_invalid || bidirectional_invalid) {
        storage_vec.append(&mut five_prime_genes);
    }

    if !(three_prime_invalid || bidirectional_invalid) {
        storage_vec.append(&mut three_prime_genes);
    }
    
    // Similar to code above; removes duplicates via hashmapping against sequence location
    let mut unique_operon_elements: HashMap<(usize, usize, String), SearchGene> = HashMap::with_capacity(storage_vec.len());
    for neighbor in storage_vec.into_iter() {
        let loc = (neighbor.linear_location.start_bound, 
                   neighbor.linear_location.end_bound,
                   neighbor.linear_location.replicon.accession_id.clone()
                );
        unique_operon_elements.insert(loc, neighbor);
    }
    let regulated_genes = unique_operon_elements.into_values().collect::<Vec<SearchGene>>();

    match regulated_genes.len() {
        0 => None,
        _ => {
            Some(regulated_genes)
        }
    }
}

impl<'genome, 'blast, 'query> SearchGenome<'genome, 'blast, 'query> {
    
    // After linking every operator to its operon, this function will
    // use that data to determine the OperatorType and OperatorClass
    fn update_operator_class_and_type(&mut self) {

        // Exit call if genome has no operons
        if let None = &self.operons {
            return
        }

        // Exit call if genome has no operators
        let operators = match &mut self.operators {
            Some(x) => x,
            None => return
        };

        // Update every operator in the SearchGenome
        for operator in operators {

            // Sort by finding a distance of Element A to operator
            let sort_routine = |a: &search::SearchGene, b: &search::SearchGene| {
                
                // Find distance of a to input gene
                let a_dist = match a.relative_to(operator) {
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
                let b_dist = match b.relative_to(operator) {
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

            // Pull full operon structure plus just a vector of its genes
            let (full_operon, mut operon_genes) = match &operator.operon {
                Some(operon) => {
                    (operon, operon.genes.clone())
                },

                // If operator has no operon, set defaults and move to next operator for analysis
                None => {   
                            operator.operator_type = search::OperatorType::Exterior;
                            operator.operator_class = search::OperatorClass::D;
                            continue
                        },
            };

            // Sort genes in operon by their distance from operator 
            operon_genes.sort_by(sort_routine);

            // Update operator's OperatorType based on it's spatial relationship
            // with the gene closest to it on its corresponding the operon
            let operator_type = match operator.relative_to(&operon_genes[0]) {
                SpatialRelationship::Overlap(x) => match x {
                    OverlapType::FivePrimeBoundary => search::OperatorType::Boundary,
                    OverlapType::ThreePrimeBoundary => search::OperatorType::Boundary,
                    _ => search::OperatorType::Interior,
                },
                SpatialRelationship::Neighbor(_) => search::OperatorType::Exterior,
                SpatialRelationship::None => search::OperatorType::Unknown,
            };

            // Update operator's OperatorClass based on it's spatial relationship
            // with the closest gene on its operon

            // If operator is not five_prime to entire operon, OperatorClass depends on 
            // operator's relationship to the gene in the operon closest to the operator
            let non_five_prime_scenario = || -> search::OperatorClass {

                match full_operon.linear_location.strand {
                    StrandSense::Forward | StrandSense::Reverse => {
                        let operon_center = CircularGenomeLocation::new(&full_operon.linear_location).center;
                        let operator_center = CircularGenomeLocation::new(&operator.linear_location).center;
                        let mut orientation = operon_center.cross(&operator_center);
                        const ZERO: f64 = 0.0;
    
                        orientation = match full_operon.linear_location.strand {
                            StrandSense::Forward => 1.0 * orientation,
                            StrandSense::Reverse => -1.0 * orientation,
                            StrandSense::Other => ZERO,
                        };
    
                        if orientation > ZERO {
                            if operon_genes.len() == 1 {
                                search::OperatorClass::C
                            } else {
                                search::OperatorClass::CMulti
                            }
                        } else {
                            if operon_genes.len() == 1 {
                                search::OperatorClass::B
                            } else {
                                search::OperatorClass::BMulti
                            }
                        }
                    },
                    
                    // 
                    StrandSense::Other => match operator.relative_to(&operon_genes[0]) {
                        SpatialRelationship::Neighbor(NeighborType::FivePrime(dist)) => search::OperatorClass::A(dist),
                        SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary) => search::OperatorClass::A(0),
                        _ => search::OperatorClass::Unknown,
                    },
                }
            };

            let operator_class = match operator.relative_to(full_operon) {
                SpatialRelationship::Overlap(x) => match x {
                    OverlapType::FivePrimeBoundary => search::OperatorClass::A(0),
                    _ => non_five_prime_scenario(),
                },

                SpatialRelationship::Neighbor(x) => match x {
                    NeighborType::FivePrime(dist) => search::OperatorClass::A(dist),
                    NeighborType::ThreePrime(_) => non_five_prime_scenario(),
                },

                SpatialRelationship::None => search::OperatorClass::Unknown,
            };

            // Merge operators association table with association tables of its operon genes
            let mut association_table: HashSet<BlastAssociationType> = HashSet::new(); 
            for gene in operon_genes {

                // Update operator association table based on operon association table
                match gene.blast_association {
                    Some(x) => {
                        association_table = association_table.union(&x).cloned().collect();
                    },
                    None => continue,
                }
            }

            operator.operator_class = operator_class;
            operator.operator_type = operator_type;
            operator.blast_association = match association_table.len() {
                0 => None,
                _ => Some(association_table),
            }
        }
    }

    // For every operator in a genome, locate its corresponding operon
    fn locate_operons(&mut self) {

        // Take read-only reference to all operators in the genome; if no operators, return
        let operators = match &mut self.operators {
            None => return,
            Some(x) => x,
        };

        // Calculate and store operon analysis result for each operator in genome
        let mut storage: Vec<Operon> = Vec::with_capacity(operators.len());
        let proto_operons = operators.iter()
                                     .map(|operator| predict_operon(operator, &self.genes))
                                     .collect::<Vec<Option<Vec<SearchGene>>>>();

        // Convert each list of potential operon genes into a formal operon struct
        for (proto_operon, operator) in proto_operons.iter().zip(operators.into_iter()) {
            match proto_operon {
                Some(x) => {
                    // The broader match statement containing this arm ensures that
                    // the vector entering into this proto-operon has at least 1 element;
                    // because of that guarantee, we are clear to use an '.unwrap()' call
                    let operon_genes = x.clone();
                    let operon = Operon::new(operon_genes).unwrap();
                    operator.operon = Some(operon.clone());
                    storage.push(operon);
                },
                None => continue,
            };
        }

        // Remove duplicates by hashing operon against its linear_genome_location; two operons
        // sharing the exact same linear genome location will only be counted once
        let mut tmp: HashMap<LinearGenomeLocation, Operon> = HashMap::with_capacity(storage.len());
        for mut operon in storage {
            tmp.insert(operon.linear_location.clone(), operon);
        } 

        // Collect hashed set of operons into a vector then sort operons by their start index (low -> high)
        let mut storage = tmp.into_values().collect::<Vec<Operon>>();
        let start_bound_sort = |a: &Operon, b: &Operon| {
            let a_start = a.linear_location.start_bound;
            let b_start = b.linear_location.start_bound;
            a_start.cmp(&b_start)
        };
        storage.sort_by(start_bound_sort);

        // Store list of operons in SearchGenome struct
        self.operons = Some(storage);

        // After assigning list of operons to each operator, update their Class/Type statuses
        self.update_operator_class_and_type();
    }

    // Search genome for every protein that has an association with CopY, but no parent operator
    #[allow(non_snake_case)]
    pub fn CopY_orphans(&self) -> Option<Vec<SearchGene>> {
        
        let CopY_genes = match self.list_all_CopY() {
            None => return None,
            Some(x) => x,
        };

        // Closure filter that keeps a gene based on whether annotated feature is NOT near an operator
        let not_operator_associated = |input_gene: &SearchGene| -> bool {

            const SEARCH_RADIUS: usize = 100;

            // Return list of genome operators; if genome has no operators, 
            // given gene cannot be associated with an operator
            let operators = match &self.operators {
                Some(x) => x,
                None => return true
            };

            // Return list of nearby operators; if gene has no nearby operators, 
            // given gene cannot be associated with an operator
            let nearby = match search::find_nearby_operators(input_gene, operators, SEARCH_RADIUS) {
                Some(x) => x,
                None => return true
            };

            // Filter by directionality; operators are only valid if NOT on 3'-prime end of a gene
            for operator in nearby {
                match operator.relative_to(input_gene) {
                    SpatialRelationship::Neighbor(NeighborType::FivePrime(_)) => return false,
                    SpatialRelationship::Neighbor(NeighborType::ThreePrime(_)) => continue,
                    SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary) => continue,
                    SpatialRelationship::None => continue,
                    SpatialRelationship::Overlap(_) => return false,
                }
            }

            // If arrive at end of the filter, there is no operator associated with gene
            true
        };

        // Check if potential orphan belongs to an existing CopY operon; if search gene 
        // overlaps a proposed CopY operon, then SearchGene is likely not an orphan
        let not_operon_associated = |input_gene: &SearchGene| -> bool {
            let mut not_associated = true;
            if let Some(x) = &self.operons {
                for operon in x {
                    if input_gene.relative_to(operon).distance() == 0 {
                        not_associated = false;
                        break;
                    }
                }
            };

            not_associated
        };


        let CopY_orphans = CopY_genes.into_iter()
                                     .filter(|x| not_operon_associated(x))
                                     .filter(|x| not_operator_associated(x))
                                     .collect::<Vec<SearchGene>>();

        match CopY_orphans.len() {
            0 => None,
            _ => Some(CopY_orphans)
        }
    }

    #[allow(non_snake_case)]
    fn list_all_CopY(&self) -> Option<Vec<SearchGene>> {
        
        let genes = match &self.genes {
            Some(x) => x.clone(),
            None => return None
        };

        fn wash(mut input: Vec<SearchGene>) -> Vec<SearchGene> {

            // Storage logistics
            let mut washer: HashMap<(usize, usize, String), SearchGene> = HashMap::new();
            
            // Filter off 'Gene' feature annotations (which are generally duplicates 
            // of the more useful CDS feature) plus the start genome region annotations
            let gene_region_filter = |x: &SearchGene| -> bool {
                (x.gene.feature_type != genome::FeatureType::Region) && 
                (x.gene.feature_type != genome::FeatureType::Gene)
            };

            // Filter off features corresponding to Gene or Region-type annotations
            input = input.into_iter().filter(|x| gene_region_filter(x)).collect::<Vec<SearchGene>>();

            // Filter duplicates by start/end bound; every possible valid (start, end) pair
            // describing a feature in the genome should be represented by just one annotation
            for gene in input {
                let start = gene.linear_location.start_bound;
                let end = gene.linear_location.end_bound;
                let replicon = gene.linear_location.replicon.accession_id.clone();
                let key = (start, end, replicon);
                washer.insert(key, gene);
            }

            washer.into_values().collect::<Vec<SearchGene>>()
        }

        // Closure filter that keeps only genes that are CopY-associated
        let CopY_associated = |x: &SearchGene| -> bool {
            match &x.blast_association {
                Some(x) => {
                    let CopY_fragment = BlastAssociationType::CopY;
                    if x.contains(&CopY_fragment) {
                        true
                    } else {
                        false
                    }
                },
                None => false,
            }
        };
       
        let CopY_genes = wash(genes).into_iter().filter(|x| CopY_associated(x)).collect::<Vec<SearchGene>>();
        match CopY_genes.len() {
            0 => None,
            _ => Some(CopY_genes),
        }
    }

    pub fn blast_associations(&self) -> Option<HashSet<BlastAssociationType>> {
        
        let genes = match &self.genes {
            Some(x) => x.clone(),
            None => return None
        };

        let mut storage: HashSet<BlastAssociationType> = HashSet::new();
        
        for gene in genes.iter() {
            let assoc = match &gene.blast_association {
                None => continue,
                Some(x) => x,
            };

            storage = storage.union(&assoc).cloned().collect();
        }

        match storage.len() {
            0 => None,
            _ => Some(storage)
        }
    }

    // Calculate total amount of space in a genome that lies b/w genes
    pub fn total_intergenic_space(&self) -> Option<(usize, usize, usize)> {
    
        let genes = match &self.genes {
            Some(x) => x,
            None => return None,
        };

        // Filter off gene annotations (which are generally duplicates of CDS) plus the start region
        let gene_region_filter = |x: &SearchGene| -> bool {
            (x.gene.feature_type != genome::FeatureType::Region) && 
            (x.gene.feature_type != genome::FeatureType::Gene) &&
            (x.gene.feature_type != genome::FeatureType::Pseudogene)
        };

        let all_genes = genes.into_iter()
                             .filter(|x| gene_region_filter(x))
                             .cloned()
                             .collect::<Vec<SearchGene>>();

        // Ensure each feature region has only one representative annotation (i.e., if a pseudogene and CDS represent
        // the same genome region in the annotation table, only one of them should be included in the neighbors list)
        let mut unique_features: HashMap<(usize, usize, String), SearchGene> = HashMap::with_capacity(all_genes.len());
        for feature in all_genes.into_iter() {
            let loc = (feature.linear_location.start_bound, 
                       feature.linear_location.end_bound,
                       feature.linear_location.replicon.accession_id.clone());
            unique_features.insert(loc, feature);
        }
        let all_relevant_unique_genes = unique_features.into_values().collect::<Vec<SearchGene>>();

        // Sort all genes into seperate vectors based on the replicon they're derived from
        let mut genes_by_replicon: Vec<(String, Vec<SearchGene>)> = Vec::new();
        for replicon_id in self.genome.replicons.keys() {
            let replicon_genes = all_relevant_unique_genes.iter()
                                                          .filter(|&x| x.gene.location.replicon_accession.eq(replicon_id))
                                                          .cloned()
                                                          .collect::<Vec<SearchGene>>();

            let tmp = (replicon_id.clone(), replicon_genes);
            genes_by_replicon.push(tmp);
        }

        // Sum up total amount of inter/intra-genic space for each replicon
        let mut net_gap_length = 0_usize;
        let mut net_genes_length_upper_bound = 0_usize;
        let mut total_relevant_genome_len = 0_usize;

        for (id, mut replicon_genes) in genes_by_replicon.into_iter() {

            // For some reason or another, certain replicons can have no genes associated with them,
            // so don't factor these replicons into the spatial calculations when encountered
            if replicon_genes.len() == 0 {
                continue;
            }

            // Sort list of remaining replicon features by start bound index
            let start_bound_sort = |a: &SearchGene, b: &SearchGene| {
                let a_start = a.linear_location.start_bound;
                let b_start = b.linear_location.start_bound;
                a_start.cmp(&b_start)
            };
            replicon_genes.sort_by(start_bound_sort);

            // Calculate space between n and n+1 replicon features for all n-1 features
            for (left, right) in replicon_genes.iter().zip(replicon_genes[1..].iter()) {
                net_gap_length += left.relative_to(right).distance();
            }

            // Calculate maximum possible length of intragenic space across genome;
            // maximum bound b/c calculation assumes genes cannot overlap
            let coding_length = replicon_genes.iter()
                                              .fold(0.0, |acc, x| acc + search::CircularGenomeLocation::new(x.get_linear_location()).arc_length);
            net_genes_length_upper_bound += coding_length.round() as usize;

            // Calculate total genome length for replicons considered
            total_relevant_genome_len += self.genome.replicons.get(&id).unwrap().fwd_strand.len();
        }
        
        // Calculate total length of genome across all replicons CONSIDERED FOR THE INTEGENIC ANALYSIS
        let output = (net_gap_length, net_genes_length_upper_bound, total_relevant_genome_len);

        match net_gap_length {
            0 => None,
            _ => Some(output)
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use cop_specific_analysis::build_list_of_genome_directories;

    #[test]
    fn test_find_operons() {

        // Import and search T4 Genome
        let refseq = PathBuf::from("tests/test_assets/assembly_summary_refseq.txt");
        let genbank = PathBuf::from("tests/test_assets/assembly_summary_genbank.txt");
        let target_dir = PathBuf::from("tests/test_assets/GCF_000006885.1_ASM688v1");
        //let target_dir = PathBuf::from("tests/test_assets/GCF_016127235.1_ASM1612723v1");
        let blast_data = Some(cop_specific_analysis::build_list_of_blast_result_files("blast/results"));
        let (metadata, blast, consensus) = prepare_common_resources(refseq, genbank, blast_data);
        let genome = build_genome_from_dir(&target_dir, &metadata);
        let search = search(&genome, &consensus, &blast);

        match &search.operators {
            Some(operators) => {
                println!("Num operators: {}", operators.len());
                let mut count = 0;

                for operator in operators.iter() {
                    count += 1;
                    println!("[{}] {}", count, operator);
                }
            },
            None => {}
        }

        match &search.operons {
            Some(operons) => {
                println!("\n\nNum operons: {}\n", operons.len());
                let mut count = 0;

                for operon in operons.iter() {
                    count += 1;
                    println!("[{}]", count);

                    for gene in &operon.genes {
                        println!("{}", gene.linear_location);
                    }

                    match &operon.operators {
                        Some(x) => {
                            for operator in x {
                                println!("{}", operator.linear_location);
                            }
                        }
                        None => continue,
                    }

                }
            },
            None => {}
        }

    }

    #[test]
    fn test_intergenic_space() {

        // Import and search T4 Genome
        let refseq = PathBuf::from("tests/test_assets/assembly_summary_refseq.txt");
        let genbank = PathBuf::from("tests/test_assets/assembly_summary_genbank.txt");
        let genome_dir = PathBuf::from("tests/test_assets/GCA_000008885.1_ASM888v1");
        let blast_data = Some(cop_specific_analysis::build_list_of_blast_result_files("blast/results"));
        let (metadata, blast, consensus) = prepare_common_resources(refseq, genbank, blast_data);
        let genome = build_genome_from_dir(&genome_dir, &metadata);
        let search = search(&genome, &consensus, &blast);

        search.total_intergenic_space();
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_CopY_orphans() {
        
        // Import and search T4 Genome
        let refseq = PathBuf::from("tests/test_assets/assembly_summary_refseq.txt");
        let genbank = PathBuf::from("tests/test_assets/assembly_summary_genbank.txt");
        let blast_data = Some(cop_specific_analysis::build_list_of_blast_result_files("blast/results"));
        let (metadata, blast, consensus) = prepare_common_resources(refseq, genbank, blast_data);
        let genomes_root_dir = PathBuf::from("/run/media/katanga/SSD2/organized_prokaryotic_database/annotated");

        let list_of_genomes = build_list_of_genome_directories(&genomes_root_dir);

        let now = std::time::Instant::now();
        const NUM: usize = 10;

        for (index, directory) in list_of_genomes[..NUM].iter().enumerate() {
            println!("[{}]", index+1);

            let genome = build_genome_from_dir(directory, &metadata);
            let search = search(&genome, &consensus, &blast);

            match search.CopY_orphans() {
                Some(x) => {
                    for orphan in x {
                        println!("{}", orphan.linear_location)
                    }
                },
                None => println!("No orphans."),
            }
            println!("\n")
        }

        println!("{}", now.elapsed().as_secs_f64() / (NUM as f64))
    }
}