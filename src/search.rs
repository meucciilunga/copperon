#![allow(dead_code)]
//#![allow(unused_imports)]

use crate::genome::{self, Gene, Genome, GenomeRegion, StrandSense, BlastFragment,
                    Replicon, BlastHitsTable};
use crate::permutations::SequencePermutations;
use std::f64::consts::PI;
use std::ops::{Add, Sub};
use std::collections::{HashSet, HashMap};

#[derive(Debug, Clone, PartialEq)]
enum SpatialRelationship {
    Neighbor(NeighborType),
    Overlap(OverlapType),
    None,
}

#[derive(Debug, Clone, PartialEq)]
enum NeighborType {
    FivePrime(usize),
    ThreePrime(usize),
}

#[derive(Debug, Clone, PartialEq)]
enum OverlapType {
    FivePrimeBoundary,
    ThreePrimeBoundary,
    EngulfedBy,
    ContainerOf,
    PerfectEclipse,
    Crescent,
    CrescentGap,
}

#[derive(Debug, Clone, PartialEq)]
enum GenomeObject<'blast, 'genome> {
    Gene(SearchGene<'genome>),
    Operator(Operator<'genome>),
    Blast(SearchBlastFragment<'blast, 'genome>),
}

#[derive(Debug, Clone, PartialEq)]
enum OperatorDimension {
    Double,
    Single,
}

#[derive(Debug, Clone, PartialEq)]
struct GenomeVector(f64, f64);

impl Add for &GenomeVector {
    type Output = GenomeVector;

    fn add(self, rhs: Self) -> GenomeVector {
        let x = self.0 + rhs.0;
        let y = self.1 + rhs.1;

        GenomeVector(x,y)
    }
}

impl Sub for &GenomeVector {
    type Output = GenomeVector;

    fn sub(self, rhs: Self) -> GenomeVector {
        let x = self.0 - rhs.0;
        let y = self.1 - rhs.1;

        GenomeVector(x,y)
    }
}

impl GenomeVector {

    // Converts a genome position in linear space to circular space, based on the size of the circular genome / replicon;
    // because bounds can only fall on real numbers, the std version of this method operates only on usize, since all ordinal
    // indicies are usize, and the bound indicies are derived directly from those
    fn new(bound_index: usize, replicon_size: usize) -> GenomeVector {

        let n = (bound_index % replicon_size) as f64;
        let theta = (2.0 * PI * n) / (replicon_size as f64);

        let x = theta.cos();
        let y = theta.sin();

        GenomeVector(x, y)
    }

    // Calculates the length of a vector
    fn magnitude(&self) -> f64 {
        (self.0).hypot(self.1)
    }

    // Stretches a vector by a scalar factor
    fn scalar(&self, scalar: f64) -> GenomeVector {
        let (x,y) = (scalar * self.0, scalar * self.1);
        GenomeVector(x,y)
    }

    // Given a vector, returns a vector with same direction but of length 1
    fn normalize(&self) -> GenomeVector {
        let factor = self.magnitude().recip();
        self.scalar(factor)
    }

    // Returns cross product of two vectors where AxB == A.cross(B)
    fn cross(&self, other: &Self) -> f64 {
        let ad = self.0 * other.1;
        let bc = other.0 * self.1;
        let pdt = ad - bc;
    
        pdt
    }

    // Returns dot product of two vectors where A⋅B == A.dot(B)
    fn dot(&self, other: &Self) -> f64 {
        let a1b1 = self.0 * other.0;
        let a2b2 = self.1 * other.1;
        a1b1 + a2b2
    }

    // Returns shortest angle (in radians) between two vectors
    fn angle(&self, other: &Self) -> f64 {

        // OLD ALGORITHM, based on simple linear algebra -- Keeping in codebase, commented out, for future reference;
        // I compared the results I got from this old algorithm with new one below for calculating the values of arc lengths,
        // which depends on accurately knowing the angle between two vectors; the new algorithm has a delta-value (i.e., 
        // the magnitude of the difference between real value and calculate value) significantly lower than this old algorithm
        
        // let theta = self.dot(&other) / (self.magnitude() * other.magnitude());
        // theta.acos()

        // NEW ALGORITHM, taken from: Miscalculating Area and Angles of a Needle-like Triangle (2014) by Prof. W. Kahan:
        // http://http.cs.berkeley.edu/~wkahan/Triangle.pdf
        // Algorithm was also discussed on this SciComp stackexchange thread, which is where I originally found it:
        // https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors
        // While writing unit tests, I found that while the old algorithm was 'good enough', this one produced much smaller
        // arc-length delta-values, so I obviously elected to use this algorithm instead of the old one.
        
        let a = self.magnitude();
        let b = other.magnitude();
        let c = (self - other).magnitude();

        let mu = if (b >= c) && (c >= 0.0) {
            c - (a - b)
        } else if (c > b) && (b >= 0.0) {
            b - (a - c)
        } else {
            panic!("ERROR: Invalid triangle!")
        };

        let top = ((a - b) + c) * mu;
        let bottom = (a + (b + c)) * ((a - c) + b);
        let theta = (top / bottom).sqrt();
        let angle = 2.0 * theta.atan();

        angle
    }   
}

// Defines a location on a linear genome by its BOUNDARY INDICIES
// NOT BY ITS ORDINAL INDICIES; BOUNDARY INDICIES ∈ ℝ while ORDINAL INDICIES ∈ ℕ;
//          START BOUND == START INDEX - 1; END BOUND == END INDEX
// NOTE: implementation of the reference structure of this struct is such that 
// you need to have a genome instance before you can have a LinearGenomeLocation; 
// this makes some practical and intuitive sense as a GenomeLocation (whether circular 
// or linear) only really makes sense when tied to a specific genome. It wouldn't  
// make sense to have a location pointing to a place on an object that doesn't exist
#[derive(Debug, Clone, PartialEq)]
struct LinearGenomeLocation<'genome> {
    replicon: &'genome Replicon,
    strand: StrandSense,
    start_bound: usize,
    end_bound: usize,
}

impl<'genome> LinearGenomeLocation<'genome> {
    fn new(input: &GenomeRegion, parent: &'genome Genome) -> LinearGenomeLocation<'genome> {
        let genome_err = "ERROR: could not find replicon denoted in GenomeRegion in the parent Genome supplied!";
        let replicon = parent.replicons.get(&input.replicon_accession).expect(genome_err);
        let strand = input.replicon_strand.clone();

        let left = input.start_index_ord - 1;
        let right = input.end_index_ord;

        LinearGenomeLocation {
            replicon,
            strand,
            start_bound: left,
            end_bound: right,
        }
    }
}

// Simply holds a linear genome location converted to its circular form
#[derive(Debug, Clone, PartialEq)]
struct CircularGenomeLocation<'genome> {
    replicon: &'genome Replicon,
    strand: StrandSense,
    start_unit_vec: GenomeVector,
    end_unit_vec: GenomeVector,
    center: GenomeVector,
    arc_length: f64
}

impl<'genome> CircularGenomeLocation<'genome> {

    // Simply converts a linear genome location to its circular equivalent
    fn new(input: &LinearGenomeLocation<'genome>) -> CircularGenomeLocation<'genome> {

        // Replicon logistics
        let replicon: &'genome Replicon = input.replicon;
        let replicon_len = replicon.fwd_strand.len();
        let replicon_radius = (replicon_len as f64) / (2.0 * PI);
        let strand = input.strand.clone();

        // Circular positioning logistics
        let start_unit_vec = GenomeVector::new(input.start_bound, replicon_len);
        let end_unit_vec = GenomeVector::new(input.end_bound, replicon_len);
        let arc_orientation = start_unit_vec.cross(&end_unit_vec);
        let center_unit_vec: GenomeVector;
        let arc_length: f64;
        const ZERO: f64 = 0.0;
        
        if arc_orientation > ZERO {
            center_unit_vec = (&start_unit_vec + &end_unit_vec).normalize();
            arc_length = start_unit_vec.angle(&end_unit_vec) * replicon_radius;
        } else if arc_orientation < ZERO {
            center_unit_vec = (&start_unit_vec + &end_unit_vec).normalize().scalar(-1.0);
            arc_length = (2.0 * PI - start_unit_vec.angle(&end_unit_vec)) * replicon_radius;
        } else {
            if input.start_bound == input.end_bound {
                center_unit_vec = (&start_unit_vec + &end_unit_vec).normalize();
                arc_length = start_unit_vec.angle(&end_unit_vec) * replicon_radius;
            } else {
                center_unit_vec = (&start_unit_vec + &end_unit_vec).normalize().scalar(-1.0);
                arc_length = (2.0 * PI - start_unit_vec.angle(&end_unit_vec)) * replicon_radius;
            }
        }

        CircularGenomeLocation {
            replicon,
            strand,
            start_unit_vec,
            end_unit_vec,
            center: center_unit_vec,
            arc_length
        }
    }
}

// Implement light wrappers around elements from Genome module to facilitate easier implementation 
// of LocateOnLinearGenome trait
struct SearchGenome<'genome, 'blast> {
    genome:     &'genome Genome,
    genes:      Option<Vec<SearchGene<'genome>>>,
    operators:  Option<Vec<Operator<'genome>>>,
    fragments:  Option<Vec<SearchBlastFragment<'blast, 'genome>>>
}

impl<'genome, 'blast, 'seq> SearchGenome<'genome, 'blast> {
    fn new(genome: &'genome Genome, permutations: &'seq SequencePermutations, blast_tables: Option<&'blast Vec<BlastHitsTable>>) -> SearchGenome<'genome, 'blast> {
        
        // If genome has gene annotation data available...
        let genes = match &genome.genes {
            Some(genes) => {
                let mut list_of_searchable_genes: Vec<SearchGene> = Vec::with_capacity(genes.len());

                for gene in genes.iter() {
                    let new_search_gene = SearchGene {
                        gene,
                        linear_location: LinearGenomeLocation::new(&gene.location, &genome),
                        blast_association: None,
                    };
                    list_of_searchable_genes.push(new_search_gene);
                }

                Some(list_of_searchable_genes)
            }
            None => None,
        };

        // If BLAST result tables have been provided...
        let fragments = if let Some(tables) = blast_tables {

            let mut genome_proto_fragments: Vec<(&BlastFragment, genome::BlastAssociationType)> = Vec::new();
            let replicon_ids: Vec<String> = genome.replicons.keys().map(|x| x.clone()).collect();

            // Check each table to see if it contains a list of fragments for a given replicon
            // in the genome; repeat this process for every replicons in the genome
            for id in replicon_ids.iter() {
                for table in tables.iter() {
                    match table.table.get(id) {
                        Some(list_of_fragments) => {
                            // Collect (references to) all blast fragments relevant to this particular
                            // replicon paired with the type of BLAST search they came from
                            let mut tmp = list_of_fragments.into_iter()
                                                           .map(|x| (x, table.association.clone()))
                                                           .collect::<Vec<(&BlastFragment, genome::BlastAssociationType)>>();
                            genome_proto_fragments.append(&mut tmp);
                        },
                        None => continue,
                    }
                }
            }

            // Convert the protofragments into a list of SearchBlastFragments
            match genome_proto_fragments.len() {
                0 => None,
                _ => {

                    let mut storage: Vec<SearchBlastFragment> = Vec::new();

                    for proto_fragment in genome_proto_fragments.iter() {
                        
                        let region = &proto_fragment.0.location;
                        let linear_location = LinearGenomeLocation::new(region, genome);
                        let blast_type = proto_fragment.1.clone();

                        let new_fragment = SearchBlastFragment {
                            blast: proto_fragment.0,
                            linear_location,
                            blast_type,
                        };

                        storage.push(new_fragment);
                    }

                    Some(storage)
                }
            }

        } else {
            None
        };

        // Assemble enough of genome to complete BLAST linking
        let mut proto_search_genome = SearchGenome {
            genome,
            genes,
            fragments,
            operators: None,
        };

        // Link BLAST results to genes
        blast_results_linker(&mut proto_search_genome);

        // Search genome for operators
        find_genome_operators(&mut proto_search_genome, permutations);

        proto_search_genome
    }
}

#[derive(Debug, Clone, PartialEq)]
struct SearchGene<'genome> {
    gene: &'genome Gene,
    linear_location: LinearGenomeLocation<'genome>,
    blast_association: Option<HashSet<genome::BlastAssociationType>>,
}

#[derive(Debug, Clone, PartialEq)]
struct SearchBlastFragment<'blast, 'genome> {
    blast: &'blast BlastFragment,
    linear_location: LinearGenomeLocation<'genome>,
    blast_type: genome::BlastAssociationType,
}

// Operator Data Structre + Related Methods
#[derive(Debug, Clone, PartialEq)]
struct Operator<'genome> {
    linear_location: LinearGenomeLocation<'genome>,
    seq: String,
    dimension: OperatorDimension,
}

impl<'genome> genome::GetSequence for Operator<'genome> {
    fn get_sequence(&self) -> String {
        self.seq.clone()
    }
}

impl<'genome> genome::ReverseComplement for Operator<'genome> {}

// Implement locate on LinearGenomeTraits for all searchable elements
trait LocateOnLinearGenome {
    fn get_linear_location(&self) -> &LinearGenomeLocation;
}

impl<'genome> LocateOnLinearGenome for SearchGene<'genome> {
    fn get_linear_location(&self) -> &LinearGenomeLocation {
        &self.linear_location
    }
}

impl<'genome> LocateOnLinearGenome for Operator<'genome> {
    fn get_linear_location(&self) -> &LinearGenomeLocation {
        &self.linear_location
    }
}

impl<'blast, 'genome> LocateOnLinearGenome for SearchBlastFragment<'blast, 'genome> {
    fn get_linear_location(&self) -> &LinearGenomeLocation {
        &self.linear_location
    }
}

// How object 'self' is related to a given origin; turns spatial relationship into a method call on objects
trait Relationship: LocateOnLinearGenome {
    fn relationship_with<T: LocateOnLinearGenome>(&self, origin: &T) -> SpatialRelationship {
        let origin = origin.get_linear_location();
        let target = &self.get_linear_location();
        spatial_relationship(origin, target)
    }
}
impl<'genome> Relationship for Operator<'genome> {}
impl<'genome> Relationship for SearchGene<'genome> {}
impl<'blast, 'genome> Relationship for SearchBlastFragment<'blast, 'genome> {}

// Given a genome and all the sequence permutations of an (operator) consensus sequence,
// locate all places in the genome where some variant of the consensus sequence is found.
fn find_genome_operators<'genome>(target_genome: &mut SearchGenome<'genome, '_>, query: &SequencePermutations) {

    // Define chunk length based on size of first element in sequence permutation list
    let chunk_len: usize = query.sequences[0].len();

    // Convert list of sequences into a hashset
    let permutations: HashSet<&str> = query.sequences.iter().map(|x| &x[..]).collect();

    // Private function for performing search operation against a single sequence
    fn locate_matching_subsequences(chunk_len: usize, queries: &HashSet<&str>, subject_seq: String, replicon_len: usize, strand: StrandSense) -> Option<Vec<(String, (usize, usize))>> { 

        let search_len = subject_seq.len() - chunk_len;
        let mut results: Vec<(String, (usize, usize))> = Vec::new();

        for start_bound_index in 0..=search_len {
            let start = start_bound_index;
            let end = start + chunk_len;

            // Check if subject-derived 'word' (16-mer chunk) is in the database of known permutations 
            let valid = queries.get(&subject_seq[start..end]);

            match valid {
                Some(&seq) => {
                    let mut ord_start = start_bound_index + 1;
                    let mut ord_end = ord_start + chunk_len - 1;
                    let ord_indicies;

                    match strand {
                        StrandSense::Forward | StrandSense::Other => {},
                        StrandSense::Reverse => {
                            let ord_start_tmp = (-1 * ord_start as i64).rem_euclid(replicon_len as i64) as usize;
                            let ord_end_tmp = (-1 * ord_end as i64).rem_euclid(replicon_len as i64) as usize;

                            ord_start = ord_end_tmp + 1;
                            ord_end = ord_start_tmp + 1;
                        },
                    }

                    ord_indicies = (ord_start, ord_end);
                    let new_entry = (seq.to_string(), ord_indicies);
                    results.push(new_entry)
                },
                None => continue,
            };
        }

        match results.len() {
            0 => None,
            _ => Some(results)
        }
    }

    fn process_proto_operator_data<'genome>(finds: Vec<(String, (usize, usize))>, direction: StrandSense, parent_genome: &SearchGenome<'genome, '_>, parent_replicon: &Replicon) -> Vec<Operator<'genome>> {

        let mut storage: Vec<Operator> = Vec::new();

        for proto_operator in finds {
            let tmp_genome_region = GenomeRegion {
                replicon_accession: parent_replicon.accession_id.clone(),
                replicon_strand: direction.clone(),
                start_index_ord: proto_operator.1.0,
                end_index_ord: proto_operator.1.1,
            };
            
            let new_operator = Operator {
                linear_location: LinearGenomeLocation::new(&tmp_genome_region, parent_genome.genome),
                seq: proto_operator.0,
                dimension: OperatorDimension::Single,
            };

            storage.push(new_operator);
        }

        storage
    }

    // Storage Logistics
    let mut operator_list: Vec<Operator<'genome>> = Vec::new();

    // Perform operator search on each replicon, parse every find into an Operator struct, then store Operators in operator_list
    for replicon in target_genome.genome.replicons.values() {
        
        // Append n basepairs (where n = chunk-length) from front of sequence
        // to end to allow for searching over circular genome discontinuity, 
        // then search for operators
        let mut fwd_modified = replicon.fwd_strand.clone();
        let fwd_wrap_around_extension = fwd_modified[0..chunk_len].to_string();
        fwd_modified.push_str(&fwd_wrap_around_extension[..]);
        let fwd_matches = locate_matching_subsequences(chunk_len, &permutations, fwd_modified, replicon.fwd_strand.len(), StrandSense::Forward);

        // Repeat process for reverse strand
        let mut rev_modified = replicon.rev_strand.clone();
        let rev_wrap_around_extension = rev_modified[0..chunk_len].to_string();
        rev_modified.push_str(&rev_wrap_around_extension[..]);
        let rev_matches = locate_matching_subsequences(chunk_len, &permutations, rev_modified, replicon.rev_strand.len(), StrandSense::Reverse);

        if let Some(finds) = fwd_matches {
            let direction = StrandSense::Forward;
            operator_list.append(&mut process_proto_operator_data(finds, direction, target_genome, replicon));
        }

        if let Some(finds) = rev_matches {
            let direction = StrandSense::Reverse;
            // Reverse order in which reverse strand operators are listed before appending
            operator_list.append(&mut process_proto_operator_data(finds, direction, target_genome, replicon).into_iter().rev().collect());
        }
    }

    // Mark any operators in the list as 'two-dimensional' if their start index appears more than once
    let operator_start_indicies: Vec<usize> = operator_list.iter().map(|x| x.linear_location.start_bound).collect();
    for operator in operator_list.iter_mut() {
        let start_bound = operator.linear_location.start_bound;
        let count = operator_start_indicies.iter().filter(|&n| *n == start_bound).count();
        if count > 1 {
            operator.dimension = OperatorDimension::Double
        }
    }

    // Remove any duplicate entries, defined as two operator who have the same forward strand start position
    let mut indices_hash_map: HashMap<(usize, usize), Operator> = HashMap::new();
    for operator in operator_list.into_iter().rev() { 
        // Reverse ordering in order to process reverse strand operators first so
        // that they get replaced preferentially by forward strand operators
        let key = (operator.linear_location.start_bound, operator.linear_location.end_bound);
        indices_hash_map.insert(key, operator);
    }

    // Sort operator list by start_index of the Operators
    let mut operator_list: Vec<Operator> = indices_hash_map.into_values().collect();
    operator_list.sort_by(|a,b| a.linear_location.start_bound.cmp(&b.linear_location.start_bound));

    match operator_list.len() {
        0 => target_genome.operators = None,
        _ => target_genome.operators = Some(operator_list),
    }
}

// Determine what spatial relationship two genomic elements A and B have with each other;
// this function is non-commutative, as the spatial relationship returned is always given
// as the positioning of B relative to A
#[allow(non_snake_case)]
#[allow(unused_parens)]
fn spatial_relationship(location_A: &LinearGenomeLocation, location_B: &LinearGenomeLocation) -> SpatialRelationship {   

    // Define circular genome location of each element
    let A = CircularGenomeLocation::new(location_A);
    let B = CircularGenomeLocation::new(location_B);
    let replicon_radius = (A.replicon.len as f64) / (2.0 * PI);
    // NOTE: replicon.len == circumference of the genome; 2*pi*radius = circumference

    // Determine if B is to the 'left' or 'right' of A
    // For an element A on the forward strand:
    // |AxB| > 0 => B is a 3' neighbor
    // |AxB| < 0 => B is a 5' neighbor
    // For Fwd strand: CCW orientation == 5' -> 3'
    // For Fwd strand: CCW orientation == 3' -> 5'
    let A_hat = A.center;
    let B_hat = B.center;
    let orientation = match A.strand {
        StrandSense::Forward | StrandSense::Other => A_hat.cross(&B_hat),
        StrandSense::Reverse => -A_hat.cross(&B_hat),
    }; // Invert orientation result based on which strand element A is on;
    // This inversion of the sign of the orientation X-pdt depending on  
    // which strand a given gene is on allows the spatial relationship
    // algorithm to process relationship results for both Fwd/Rev strand
    // genes without having to duplicate very similar code with only minor
    // changes accounting for strand-specific Left/Right determination

    let A_radius = A.arc_length / 2.0;
    let B_radius = B.arc_length / 2.0;
    
    // Define length of arc intervening between center of A and center of B
    let Δ_angle = A_hat.angle(&B_hat);            // Angle between centers of element A and B
    let Δ_arc_length = replicon_radius * Δ_angle; // Length of an arc = theta*radius, where theta is the angle of arc in radians

    // Calculate values of the algorithm
    const ZERO: f64 = 0.0;
    let different_replicons = (A.replicon.accession_id != B.replicon.accession_id);
    let no_overlap = Δ_arc_length >= (A_radius + B_radius);
    let engulfing_overlap = f64::max(A_radius, B_radius) > (Δ_arc_length + f64::min(A_radius, B_radius));
    let crescent = (2.0 * PI * replicon_radius - f64::min(2.0*A_radius, 2.0*B_radius)) < f64::max(2.0*A_radius, 2.0*B_radius);

    // Apply rules of the algorithm
    let output = if different_replicons {
        SpatialRelationship::None
    } else if no_overlap {

        // Arc distance should be a whole integer, so round before casting; need to round
        // because casting to usize automatically applies a floor rounding, which means
        // a number like 18361.999999999796 gets rounded to 18361 instead of 18362. Applying
        // a normal integer rounding operation before casting should stop that from happening.
        let arc_distance = (Δ_arc_length - A_radius - B_radius).round() as usize;

        if orientation > ZERO { 
            // AxB > 0 => center of B is to the LEFT of center of A (for fwd strand)
            SpatialRelationship::Neighbor(NeighborType::ThreePrime(arc_distance))
        } else { 
            // AxB < 0 => center of B is to the RIGHT of center of A
            SpatialRelationship::Neighbor(NeighborType::FivePrime(arc_distance))
        }
    } else if crescent {
        if A_radius > B_radius {
            SpatialRelationship::Overlap(OverlapType::CrescentGap)
        } else {
            SpatialRelationship::Overlap(OverlapType::Crescent)
        }
    } else if engulfing_overlap {
        if A_radius > B_radius {
            SpatialRelationship::Overlap(OverlapType::EngulfedBy)
        } else {
            SpatialRelationship::Overlap(OverlapType::ContainerOf)
        }
    } else {
        if orientation == ZERO {
            SpatialRelationship::Overlap(OverlapType::PerfectEclipse)
        } else if orientation > ZERO {
            SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary)
        } else {
            SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary)
        }
    };

    output
}

fn search_bubble<'genome, T: LocateOnLinearGenome>(query: &'genome T, search_radius: usize) -> LinearGenomeLocation {
    let mut proto_bubble = query.get_linear_location().clone();
    let proto_bubble_len = CircularGenomeLocation::new(&proto_bubble).arc_length.round() as usize;

    if !(proto_bubble_len + 2 * search_radius >= proto_bubble.replicon.len) {
        let start_i64: i64 = (proto_bubble.start_bound as i64) - (search_radius as i64);
        
        proto_bubble.start_bound = start_i64.rem_euclid(proto_bubble.replicon.len as i64) as usize;
        proto_bubble.end_bound += search_radius;
    }

    proto_bubble
}

// Given a genomic element and a list of all other elements (to be considered) on that genome,
// return a list of all genome elements within a distance 'search_radius' of the query's center
fn find_nearby_elements<'blast, 'genome, R>(query: &R, possible_elements: &'genome Vec<GenomeObject<'blast, 'genome>>, search_radius: usize)
-> Option<Vec<(&'genome GenomeObject<'blast, 'genome>, SpatialRelationship)>>
where R: LocateOnLinearGenome,
{
    // Storage logistics
    let mut nearby_list: Vec<(&'genome GenomeObject, SpatialRelationship)> = Vec::new(); 

    // Determine location of the query and its corresponding search bubble
    let query_loc = query.get_linear_location();
    let bubble_loc = search_bubble(query, search_radius);

    for element in possible_elements {

        // Pull the location of a given element
        let element_loc = match element {
            GenomeObject::Gene(gene) => gene.get_linear_location(),
            GenomeObject::Operator(operator) => operator.get_linear_location(),
            GenomeObject::Blast(blast) => blast.get_linear_location(),
        };

        // Perform relational search against the bubble, looking for elements overlapping the bubble
        match spatial_relationship(&bubble_loc, element_loc) {

            // We only care about the elements that overlap our search bubble in some way
            SpatialRelationship::Overlap(_) => {
                let relationship = spatial_relationship(query_loc, element_loc);
                let find = (element, relationship);
                nearby_list.push(find);
            },
            _ => continue, 
            // We don't care if an element is neighboring the SEARCH BUBBLE; 
            // we care ONLY about how it neighbors the query itself; we also don't
            // care about elements that are on completely different replicons
        }

    }

    match nearby_list.len() {
        0 => None,
        _ => Some(nearby_list),
    }
}

// Given a genome with genes and a BLAST table, updates all the genes in the genome 
// instance to reflect whether each is associated with any known BLAST hits
fn blast_results_linker(genome: &mut SearchGenome) {

    // Check if genome has associated blast hits; exit call if genome has no associated blasts
    let blast_hits = match &genome.fragments {
        Some(hits) => hits,
        None => return,
    };

    // Check if genome has associated gene annotations; exit call if genome has no associated genes
    let genes = match &mut genome.genes {
        Some(genes) => genes,
        None => return,
    };

    // Check the relationship every gene has with every known BLAST hit
    for gene in genes {

        for hit in blast_hits {
            
            // If a new match is found
            if let SpatialRelationship::Overlap(_) = gene.relationship_with(hit) {
                // If gene and blast hit are determined to overlap, either create a new storage vector
                // for blast associations if one doesn't exist, or append new blast association to the
                // existing list
                gene.blast_association = match &mut gene.blast_association {
                    None => {
                        let mut new_association_table: HashSet<genome::BlastAssociationType> = HashSet::new();
                        new_association_table.insert(hit.blast_type.clone());
                        Some(new_association_table)
                    },
                    Some(association_table) => {
                        association_table.insert(hit.blast_type.clone());
                        continue
                    }
                }
            }
        }
    }
}







// UNIT TESTS
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cop_operon_specific::build_cop_permutation_table;
    use crate::permutations::SequencePermutations;
    use crate::genome::{self, GenomeRegion, FeatureType, BlastAssociationType};
    use crate::import::{self, AssemblyMetadata};
    use std::path::PathBuf;
    use std::{collections::HashMap, f64::consts::E};
    use std::time::SystemTime;

    #[test]
    fn test_genome_vector_cross_pdt() {
        let vec_1 = GenomeVector(10.0, -PI);
        let vec_2 = GenomeVector(E, -2.0*PI);
        let vec_3 = GenomeVector(-13.69, -2.1);

        let test_cross_1 = vec_1.cross(&vec_2);
        let test_cross_2 = vec_2.cross(&vec_2);
        let test_cross_3 = vec_2.cross(&vec_3);
        let test_cross_4 = vec_3.cross(&vec_2);

        // Using delta since can't accurately compare floats directly
        let delta_1 = (test_cross_1.abs() - (-54.29211885_f64).abs()).abs();
        assert!(delta_1 < 1e-8);

        assert_eq!(0.0, test_cross_2);

        let delta_3 = (test_cross_3.abs() - (-91.7251987_f64).abs()).abs();
        assert!(delta_3 < 1e-8);

        assert!(test_cross_3 == -test_cross_4);
    }

    #[test]
    fn test_genome_dot_cross_pdt() {
        let vec_1 = GenomeVector(10.0, -PI);
        let vec_2 = GenomeVector(E, -2.0*PI);
        let vec_3 = GenomeVector(-13.69, -2.1);
        let vec_4 = GenomeVector(0.0, 0.0);
        let vec_5 = GenomeVector(1.0, 1.0);

        let test_dot_1 = vec_1.dot(&vec_2);
        let test_dot_2 = vec_2.dot(&vec_2);
        let test_dot_3 = vec_2.dot(&vec_3);
        let test_dot_4 = vec_3.dot(&vec_2);
        let test_dot_5 = vec_4.dot(&vec_2);
        let test_dot_6 = vec_3.dot(&vec_1);
        let test_dot_7 = vec_5.dot(&vec_3);

        let test_answers = [test_dot_1, test_dot_2, test_dot_3, test_dot_4, test_dot_5, test_dot_6, test_dot_7];
        let expected_answers = [46.92202709, 46.8674737, -24.01858909, -24.01858909, 0.0, -130.3026554, -15.79];

        for (left, right) in expected_answers.iter().zip(test_answers.iter()) {
            let delta = (left - right).abs();
            assert!(1e-6 > delta);
        }
    }

    #[test]
    fn test_plus_minus_genome_vector_overload() {
        let vec_1 = GenomeVector(10.0, -8.0);
        let vec_2 = GenomeVector(7.5, -2.33);
        let vec_3 = GenomeVector(-13.69, -2.1777);

        // Test cases
        let test_sum_1 = &vec_1 + &vec_2;
        let test_sum_2 = &vec_2 + &vec_3;
        let test_sum_3 = &vec_3 + &vec_2;
        let test_sum_4 = &vec_1 + &vec_3;

        let test_diff_1 = &vec_1 - &vec_2;
        let test_diff_2 = &vec_2 - &vec_3;
        let test_diff_3 = &vec_3 - &vec_2;
        let test_diff_4 = &vec_1 - &vec_3;
        
        let test_solns = [
            test_sum_1,
            test_sum_2,
            test_sum_3,
            test_sum_4,
            test_diff_1,
            test_diff_2,
            test_diff_3,
            test_diff_4,
        ];

        // Solutions
        let actual_sum_1 = GenomeVector(17.5, -10.33);
        let actual_sum_2 = GenomeVector(-6.19, -4.5077);
        let actual_sum_3 = GenomeVector(-6.19, -4.5077);
        let actual_sum_4 = GenomeVector(-3.69, -10.1777);

        let actual_diff_1 = GenomeVector(2.5, -5.67);
        let actual_diff_2 = GenomeVector(21.19, -0.1523);
        let actual_diff_3 = GenomeVector(-21.19, 0.1523);
        let actual_diff_4 = GenomeVector(23.69, -5.8223);
        
        let actual_solns = [
            actual_sum_1,
            actual_sum_2,
            actual_sum_3,
            actual_sum_4,
            actual_diff_1,
            actual_diff_2,
            actual_diff_3,
            actual_diff_4
        ];


        // Comparisons
        for (left, right) in actual_solns.iter().zip(test_solns.iter()) {
            let left_diff = (left.0 - right.0).abs();
            let right_diff = (left.1 - right.1).abs();

            let test_msg = format!("test_case: {:?}", right);
            println!("{}", test_msg);
            let actual_msg = format!("actual: {:?}\n", left);
            println!("{}", actual_msg);

            assert!(1e-12 > left_diff);
            assert!(1e-12 > right_diff);
        }
    }

    #[test]
    fn test_magnitude_normalize_scale() {

        // Test inputs
        let test_vec_1 = GenomeVector(18.01, -9.76);
        let test_vec_2 = GenomeVector(-1.0, -1.0);
        let test_vec_3 = GenomeVector(299.6, 80.808080);

        // Test cases
        let test_scale_1 = test_vec_1.scalar(-10.0);
        let test_scale_2 = test_vec_2.scalar(-PI);
        let test_scale_3 = test_vec_3.scalar(2.0);

        let test_norm_1 = test_vec_1.normalize();
        let test_norm_2 = test_vec_2.normalize();
        let test_norm_3 = test_vec_3.normalize();

        let test_mag_1 = test_vec_1.magnitude();
        let test_mag_2 = test_vec_2.magnitude();
        let test_mag_3 = test_vec_3.magnitude();

        let test_scaled = [
            test_scale_1,
            test_scale_2,
            test_scale_3,
            test_norm_1,
            test_norm_2,
            test_norm_3,
        ];

        let test_mags = [
            test_mag_1,
            test_mag_2,
            test_mag_3,
        ];

        // Actual values
        let actual_scaled_1 = GenomeVector(-180.1, 97.60);
        let actual_scaled_2 = GenomeVector(PI, PI);
        let actual_scaled_3 = GenomeVector(599.2, 161.61616);

        let actual_norm_1 = GenomeVector(0.8791982467, -0.4764561292);
        let actual_norm_2 = GenomeVector(-2.0_f64.powf(-0.5), -2.0_f64.powf(-0.5));
        let actual_norm_3 = GenomeVector(0.9654971076, 0.2604137767);

        let actual_mag_1 = (test_vec_1.0.powi(2) + test_vec_1.1.powi(2)).sqrt();
        let actual_mag_2 = (test_vec_2.0.powi(2) + test_vec_2.1.powi(2)).sqrt();
        let actual_mag_3 = (test_vec_3.0.powi(2) + test_vec_3.1.powi(2)).sqrt();

        let actual_scaled = [
            actual_scaled_1,
            actual_scaled_2,
            actual_scaled_3,
            actual_norm_1,
            actual_norm_2,
            actual_norm_3,
        ];

        let actual_mags = [
            actual_mag_1,
            actual_mag_2,
            actual_mag_3,
        ];

        // Comparisons
        for (left, right) in actual_scaled.iter().zip(test_scaled.iter()) {
            let diff = left - right;
            let delta = 1e-10;
            
            println!("TEST {:?} :: {:?}", left, right);

            assert!(diff.0.abs() < delta);
            assert!(diff.1.abs()  < delta);
        }

        for (left, right) in actual_mags.iter().zip(test_mags.iter()) {
            println!("TEST {} :: {}", left, right);
            assert_eq!(left, right);
        }
    }

    #[test]
    fn test_genome_vector_angle() {
        let baseline = GenomeVector(1.0, 0.0);

        let x0 = 1_f64; // 0
        let x1 = 3_f64.sqrt() / 2.0; // pi/6
        let x2 = 2_f64.sqrt() / 2.0; // pi/4
        let x3 = 1.0 / 2.0; // pi/3
        let x4 = 0.0; // pi/2
        let base_unit_circle_values = [x0, x1, x2, x3, x4];
        let mut unit_angles = vec![0_f64, (1.0/6.0), 0.25, (1.0/3.0), 0.5, 0.5, (2.0/3.0), 0.75, (5.0/6.0), 1.0];
        let mut unit_angles_half2 = unit_angles.iter().rev().map(|&x| x).collect::<Vec<f64>>();
        unit_angles.append(&mut unit_angles_half2);

        let actual_angles = unit_angles.iter()
                                       .map(|x| x * PI)
                                       .collect::<Vec<f64>>();

        let quadrant_1 = base_unit_circle_values.iter().zip(base_unit_circle_values.iter().rev()).map(|(&x,&y)| (x,y)).collect::<Vec<(f64,f64)>>();
        let quadrant_2 = base_unit_circle_values.iter().zip(base_unit_circle_values.iter().rev()).map(|(&x,&y)| (-x,y)).rev().collect::<Vec<(f64,f64)>>();
        let quadrant_3 = base_unit_circle_values.iter().zip(base_unit_circle_values.iter().rev()).map(|(&x,&y)| (-x,-y)).collect::<Vec<(f64,f64)>>();
        let quadrant_4 = base_unit_circle_values.iter().zip(base_unit_circle_values.iter().rev()).map(|(&x,&y)| (x,-y)).rev().collect::<Vec<(f64,f64)>>();
        let quadrants = [quadrant_1, quadrant_2, quadrant_3, quadrant_4];

        let mut test_angles: Vec<f64> = vec![];
        for quadrant in quadrants {
            for (x,y) in quadrant {
                let test_vec = GenomeVector(x,y);
                let tmp = baseline.angle(&test_vec);
                test_angles.push(tmp);
            }
        }

        for (left, right) in actual_angles.iter().zip(test_angles.iter()) {
            let delta = (left - right).abs();
            assert!(1e-12 > delta);
        }
    }

    #[test]
    fn genome_vector_from_linear_bound_1() {
        const REPLICON_SIZE: usize = 16;

        // Generate all cartesian unit vectors for every bp boundary position in a small circular genome
        let test_cases = (0..=REPLICON_SIZE).map(|n| GenomeVector::new(n, REPLICON_SIZE));

        // Generate cartesian unit vectors via an angle 
        let rotate_vector = |theta: f64, input: GenomeVector| -> GenomeVector {
            let x = input.0 * theta.cos() + input.1 * theta.sin();
            let y = input.0 * theta.sin() + input.1 * theta.cos();
            GenomeVector(x,y)
        };

        let actual_cases = (0..=REPLICON_SIZE).map(|x| {
                                        let theta = 2.0*PI*(x as f64 / REPLICON_SIZE as f64);
                                        let x_hat = GenomeVector(1.0, 0.0);
                                        rotate_vector(theta, x_hat)
                                    });

        for (left, right) in actual_cases.zip(test_cases) {

            println!("Test-Value: {:?}\nActual-Value: {:?}\n", left, right);

            let difference = (&left - &right).magnitude();
            let delta = 1e-12;
            assert!(delta > difference);
        }
    }

    #[test]
    fn genome_vector_from_linear_bound_2() {
        const REPLICON_SIZE: usize = 7_231_001;
        const MAX_INDEX: usize = 7_231_001 * 4;
        const VECTOR_LEN: usize = 23;

        // Generate all cartesian unit vectors for every bp boundary position in a small circular genome
        let test_cases = (0..=MAX_INDEX).map(|n| GenomeVector::new(n, REPLICON_SIZE));

        // Generate cartesian unit vectors via an angle 
        let rotate_vector = |theta: f64, input: GenomeVector| -> GenomeVector {
            let x = input.0 * theta.cos() + input.1 * theta.sin();
            let y = input.0 * theta.sin() + input.1 * theta.cos();

            GenomeVector(x,y)
        };

        let actual_cases = (0..=MAX_INDEX).map(|x| {
                                        let theta = 2.0*PI*(x as f64 / REPLICON_SIZE as f64);
                                        let x_hat = GenomeVector(1.0, 0.0);
                                        rotate_vector(theta, x_hat)
                                    });

        for (left, right) in actual_cases.zip(test_cases) {

            let difference = (&left - &right).magnitude();
            let delta = 1e-12;

            assert!(delta > difference);
        }
    }

    #[test]
    fn genome_vector_from_linear_bound_3() {
        const REPLICON_SIZE: usize = 69_421;
        const MAX_INDEX: usize = 888_888;

        // Generate all cartesian unit vectors for every bp boundary position in a small circular genome
        let test_cases = (0..=MAX_INDEX).map(|n| GenomeVector::new(n, REPLICON_SIZE));

        // Generate cartesian unit vectors via an angle 
        let rotate_vector = |theta: f64, input: GenomeVector| -> GenomeVector {
            let x = input.0 * theta.cos() + input.1 * theta.sin();
            let y = input.0 * theta.sin() + input.1 * theta.cos();

            GenomeVector(x,y)
        };

        let matrix_transformation_calculated_cases = (0..=MAX_INDEX).map(|x| {
                                        let theta = 2.0*PI*(x as f64 / REPLICON_SIZE as f64);
                                        let x_hat = GenomeVector(1.0, 0.0);
                                        rotate_vector(theta, x_hat)
                                    });

        for (left, right) in matrix_transformation_calculated_cases.zip(test_cases) {
            let difference = (&left - &right).magnitude();
            let delta = 1e-4;
            assert!(delta > difference);
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_arc_delta_test_1bp_arc_50Mbp_genome() {

        const REPLICON_SIZE: usize = 50_000_000;

        // Calculate unit arcs (1bp) for every base-pair position in the genome
        let mut list_of_arc_lengths: Vec<f64> = Vec::new();

        for i in 0..=REPLICON_SIZE+1 {
            let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
            let upper_bound = GenomeVector::new(i+1, REPLICON_SIZE);
            let seperating_angle = lower_bound.angle(&upper_bound);
            let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
            let arc_length = seperating_angle * radius;

            list_of_arc_lengths.push(arc_length);
        }

        let delta_list = list_of_arc_lengths.iter().map(|x| (x - 1.0).abs());
        let mut max_deviation = 0_f64;
        for deviation in delta_list {
            if deviation > max_deviation {
                max_deviation = deviation
            }
        }

        println!("MAXIMUM deviation for 1bp arc in {}bp genome = {}", REPLICON_SIZE, max_deviation);
        assert!(max_deviation < 1e-6);
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_arc_delta_test_1bp_arc_25Mbp_genome() {

        const REPLICON_SIZE: usize = 25_000_000;

        // Calculate unit arcs (1bp) for every base-pair position in the genome
        let mut list_of_arc_lengths: Vec<f64> = Vec::new();

        for i in 0..=REPLICON_SIZE+1 {
            let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
            let upper_bound = GenomeVector::new(i+1, REPLICON_SIZE);
            let seperating_angle = lower_bound.angle(&upper_bound);
            let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
            let arc_length = seperating_angle * radius;

            list_of_arc_lengths.push(arc_length);
        }

        let delta_list = list_of_arc_lengths.iter().map(|x| (x - 1.0).abs());
        let mut max_deviation = 0_f64;
        for deviation in delta_list {
            if deviation > max_deviation {
                max_deviation = deviation
            }
        }

        println!("MAXIMUM deviation for 1bp arc in {}bp genome = {}", REPLICON_SIZE, max_deviation);
        assert!(max_deviation < 1e-6);
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_arc_over_continuity() {

        const REPLICON_SIZE: usize = 100_000_000;

        // Standard Variant
        let lower_bound_1 = GenomeVector::new(95_000_000, REPLICON_SIZE);
        let upper_bound_1 = GenomeVector::new(5_000_000, REPLICON_SIZE);
        let seperating_angle_1 = lower_bound_1.angle(&upper_bound_1);
        let radius_1 = (REPLICON_SIZE as f64) / (2.0 * PI);
        let arc_length_1 = seperating_angle_1 * radius_1;

        // Wrap-around Variant
        let lower_bound_2 = GenomeVector::new(95_000_000, REPLICON_SIZE);
        let upper_bound_2 = GenomeVector::new(105_000_000, REPLICON_SIZE);
        let seperating_angle_2 = lower_bound_2.angle(&upper_bound_2);
        let radius_2 = (REPLICON_SIZE as f64) / (2.0 * PI);
        let arc_length_2 = seperating_angle_2 * radius_2;

        // Compare values
        assert_eq!(lower_bound_1, lower_bound_2);
        assert_eq!(upper_bound_1, upper_bound_2);
        assert_eq!(seperating_angle_1, seperating_angle_2);
        assert_eq!(radius_1, radius_2);
        assert_eq!(arc_length_1, arc_length_2);
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_feature_over_continuity_linear_and_circular() {

        let (T4, LACTO) = import_tigr4_lacto_genome();

        // T4
        let test_region_T4 = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 2_160_000,
            end_index_ord: 2_572,
        };

        let t4_linear_test = LinearGenomeLocation::new(&test_region_T4, &T4);

        let actual_t4_linear = LinearGenomeLocation {
            replicon: T4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_bound: 2_159_999,
            end_bound: 2_572,
        };

        let t4_circular = CircularGenomeLocation::new(&t4_linear_test);

        let start_vec = GenomeVector::new(2_159_999, T4.replicons.get(&"NC_003028.3".to_string()).unwrap().len);
        let end_vec = GenomeVector::new(2_572, T4.replicons.get(&"NC_003028.3".to_string()).unwrap().len);
        let actual_t4_circ = CircularGenomeLocation {
            replicon: T4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_unit_vec: start_vec.clone(),
            end_unit_vec: end_vec.clone(),
            center: (&start_vec + &end_vec).normalize(),
            arc_length: 3415.00,
        };

        assert_eq!(actual_t4_linear, t4_linear_test);
        actual_t4_circ.test_compare(&t4_circular);


        // LACTO
        let test_region_LACTO = GenomeRegion {
            replicon_accession: "NC_002662.1".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 2_365_000,
            end_index_ord: 2_572,
        };

        let lacto_linear_test = LinearGenomeLocation::new(&test_region_LACTO, &LACTO);

        let actual_lacto_linear = LinearGenomeLocation {
            replicon: LACTO.replicons.get(&"NC_002662.1".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_bound: 2_364_999,
            end_bound: 2_572,
        };

        let lacto_circular = CircularGenomeLocation::new(&actual_lacto_linear);

        let start_vec = GenomeVector::new(2_364_999, LACTO.replicons.get(&"NC_002662.1".to_string()).unwrap().len);
        let end_vec = GenomeVector::new(2_572, LACTO.replicons.get(&"NC_002662.1".to_string()).unwrap().len);
        let actual_lacto_circ = CircularGenomeLocation {
            replicon: LACTO.replicons.get(&"NC_002662.1".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_unit_vec: start_vec.clone(),
            end_unit_vec: end_vec.clone(),
            center: (&start_vec + &end_vec).normalize(),
            arc_length: 3162.00,
        };

        assert_eq!(actual_lacto_linear, lacto_linear_test);
        actual_lacto_circ.test_compare(&lacto_circular);
    }

    #[test]
    #[allow(non_snake_case)]
    #[allow(unused_parens)]
    fn genome_arc_compare_wrap_arounds() {
        const REPLICON_SIZE: usize = 25_000_000;

        // Calculate unit arcs (1bp) for every base-pair position in the genome
        let mut std_arcs: Vec<f64> = Vec::new();
        for i in 0..REPLICON_SIZE {
            let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
            let upper_bound = GenomeVector::new(i+1, REPLICON_SIZE);
            let seperating_angle = lower_bound.angle(&upper_bound);
            let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
            let arc_length = seperating_angle * radius;

            std_arcs.push(arc_length);
        }

        // Calculate unit arcs (1bp) for every base-pair position in the genome
        let mut wrapped_arcs: Vec<f64> = Vec::new();
        for i in REPLICON_SIZE..2*REPLICON_SIZE {
            let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
            let upper_bound = GenomeVector::new(i+1, REPLICON_SIZE);
            let seperating_angle = lower_bound.angle(&upper_bound);
            let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
            let arc_length = seperating_angle * radius;

            wrapped_arcs.push(arc_length);
        }

        // Compare that arcs for wrap-around indicies are the same for std indicies
        for (index, (std, wrapped)) in std_arcs.iter().zip(wrapped_arcs.iter()).enumerate() {

            // Using exact comparisons b/w floats b/c the output values should be EXACTLY the same
            assert_eq!(*std, *wrapped);

            // Print last 5 comparisons
            if (index >= 24_999_995) {
                println!("std: {:?} | wrapped: {:?}", std, wrapped)
            }
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_arc_delta_test_varying_arc_100Mbp_genome() {

        const REPLICON_SIZE: usize = 100_000_000;
        let arc_length_samples = [1, 2, 3, 10, 25, 100_000, 1_000_000, 10_000_000];

        let mut list_of_max_deviations: Vec<f64> = Vec::new();
        for arc_len in arc_length_samples  {
            
            // Calculate varying length arcs (1bp) for every base-pair position in the genome
            let mut list_of_arc_lengths: Vec<f64> = Vec::new();
            let mut list_of_indicies: Vec<usize> = Vec::new();

            for i in 0..=REPLICON_SIZE+1 {
                let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
                let upper_bound = GenomeVector::new(i+arc_len, REPLICON_SIZE);
                let seperating_angle = lower_bound.angle(&upper_bound);
                let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
                let arc_length = seperating_angle * radius;

                list_of_arc_lengths.push(arc_length);
                list_of_indicies.push(i);
            }

            let delta_list = list_of_arc_lengths.iter().map(|x| (x - arc_len as f64).abs());
            let mut max = 0_f64;
            let mut max_index = 0;
            let mut min = f64::INFINITY;
            let mut min_index = 0;

            for (index, deviation) in delta_list.enumerate() {
                if deviation > max {
                    max = deviation;
                    max_index =  index;
                }
                if deviation < min {
                    min = deviation;
                    min_index =  index;
                }
            }

            println!("MAXIMUM Deviation for {}bp arc length in {}bp genome: {} at pos. {}bp", arc_len, REPLICON_SIZE, max, max_index);
            println!("MINIMUM Deviation for {}bp arc length in {}bp genome: {} at pos. {}bp\n", arc_len, REPLICON_SIZE, min, min_index);
            list_of_max_deviations.push(max);
        }

        let mut max_deviation = 0_f64;
        for deviation in list_of_max_deviations {
            if deviation > max_deviation {
                max_deviation = deviation
            }
        }

        println!("MAXIMUM deviation for various arc lengths  {:?} in {}bp genome across = {}", arc_length_samples, REPLICON_SIZE, max_deviation);
        assert!(max_deviation < 1e-6);
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_arc_delta_test_16bp_arc_20Mbp_genome() {

        const REPLICON_SIZE: usize = 20_000_000;
        let arc_length_samples = [1, 16];

        let mut list_of_max_deviations: Vec<f64> = Vec::new();
        for arc_len in arc_length_samples  {
            // Calculate varying length arcs (1bp) for every base-pair position in the genome
            let mut list_of_arc_lengths: Vec<f64> = Vec::new();

            for i in 0..=REPLICON_SIZE+1 {
                let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
                let upper_bound = GenomeVector::new(i+1, REPLICON_SIZE);
                let seperating_angle = lower_bound.angle(&upper_bound);
                let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
                let arc_length = seperating_angle * radius;

                list_of_arc_lengths.push(arc_length);
            }

            let delta_list = list_of_arc_lengths.iter().map(|x| (x - 1.0).abs());
            let mut max_deviation = 0_f64;
            for deviation in delta_list {
                if deviation > max_deviation {
                    max_deviation = deviation
                }
            }

            println!("Max Deviation for {}bp arc length in {}bp genome: {}", arc_len, REPLICON_SIZE, max_deviation);
            list_of_max_deviations.push(max_deviation);
        }

        let mut max_deviation = 0_f64;
        for deviation in list_of_max_deviations {
            if deviation > max_deviation {
                max_deviation = deviation
            }
        }

        println!("MAXIMUM deviation for 1bp arc in {}bp genome = {}", REPLICON_SIZE, max_deviation);
        assert!(max_deviation < 1e-6);
    }

    // Several of the following tests require a test Genome to work with, so this performs that import operation
    fn import_tigr4_lacto_genome() -> (Genome, Genome) {

        // GLOBAL
        // Import metadata
        let now = SystemTime::now();
        let refseq_metadata_file = PathBuf::from("tests/test_assets/assembly_summary_refseq.txt");
        let genbank_metadata_file = PathBuf::from("tests/test_assets/assembly_summary_genbank.txt");
        let refseq_meta = import::import_database_summary(refseq_metadata_file).into_iter();
        let genbank_meta = import::import_database_summary(genbank_metadata_file).into_iter();
        let full_metadata: HashMap<String, AssemblyMetadata> = refseq_meta.chain(genbank_meta).collect();
        println!("Assembly databases import time: {}ms", now.elapsed().unwrap().as_millis());

        // TIGR4
        // Paths to relevant files to TIGR4
        let now = SystemTime::now();
        let assembly_name = "GCF_000006885.1_ASM688v1".to_string();
        let genome_seq_file = PathBuf::from("tests/test_assets/GCF_000006885.1_ASM688v1/GCF_000006885.1_ASM688v1_genomic.fna");
        let genome_annotation_file = PathBuf::from("tests/test_assets/GCF_000006885.1_ASM688v1/GCF_000006885.1_ASM688v1_genomic.gff");

        // Import data relevant to TIGR4
        let protogenome = import::parse_genome_sequence(&assembly_name, genome_seq_file);
        let protogenes = import::parse_genome_annotation(genome_annotation_file);

        // Build TIGR4 genome
        let asm_pull_error = "ERROR: could not find assembly name in metadata database!";
        let metadata = full_metadata.get(&assembly_name).expect(asm_pull_error);
        let genes = Some(genome::Gene::convert_annotation_entry_list(&protogenes));
        let t4 = genome::Genome::from_proto(protogenome, metadata, genes);
        println!("TIGR4 proto-genome import time: {}ms", now.elapsed().unwrap().as_millis());

        // LACTO
        // Paths to relevant files to lacto
        let assembly_name = "GCF_000006865.1_ASM686v1".to_string();
        let genome_seq_file = PathBuf::from("tests/test_assets/GCF_000006865.1_ASM686v1/GCF_000006865.1_ASM686v1_genomic.fna");
        let genome_annotation_file = PathBuf::from("tests/test_assets/GCF_000006865.1_ASM686v1/GCF_000006865.1_ASM686v1_genomic.gff");

        // Import data relevant to lacto
        let protogenome = import::parse_genome_sequence(&assembly_name, genome_seq_file);
        let protogenes = import::parse_genome_annotation(genome_annotation_file);

        // Build lacto genome
        let asm_pull_error = "ERROR: could not find assembly name in metadata database!";
        let metadata = full_metadata.get(&assembly_name).expect(asm_pull_error);
        let genes = Some(genome::Gene::convert_annotation_entry_list(&protogenes));
        let lacto = genome::Genome::from_proto(protogenome, metadata, genes);

        (t4, lacto)
    }

    #[test]
    #[allow(non_snake_case)]
    fn derive_circular_genome_from_linear_1() {

        let imported_genomes = import_tigr4_lacto_genome();

        let TIGR4 = imported_genomes.0;
        let test_genes = TIGR4.genes.clone().unwrap();
        let test_genes = [test_genes[0].clone(), test_genes[2].clone(), test_genes[1212].clone(), test_genes[4078].clone()];
        let mut test_locs: Vec<CircularGenomeLocation> = Vec::with_capacity(test_genes.len());

        for test_gene in &test_genes {
            let linear_location = LinearGenomeLocation::new(&test_gene.location, &TIGR4);
            let circ_location = CircularGenomeLocation::new(&linear_location);
            test_locs.push(circ_location);
        }

        let circ_loc_1 = CircularGenomeLocation {
            replicon: TIGR4.replicons.get("NC_003028.3").unwrap(),
            strand: StrandSense::Forward,
            start_unit_vec: GenomeVector(1.0, 0.0),
            end_unit_vec: GenomeVector(1.0, 0.0),
            center: GenomeVector(-1.0, 0.0),
            arc_length: 2160842.0,
        };
        let circ_loc_2 = CircularGenomeLocation {
            replicon: TIGR4.replicons.get("NC_003028.3").unwrap(),
            strand: StrandSense::Forward,
            start_unit_vec: GenomeVector::new(196, 2160842),
            end_unit_vec: GenomeVector::new(1558, 2160842),
            center: GenomeVector::new(877, 2160842),
            arc_length: 1362.0,
        };
        let circ_loc_3 = CircularGenomeLocation {
            replicon: TIGR4.replicons.get("NC_003028.3").unwrap(),
            strand: StrandSense::Forward,
            start_unit_vec: GenomeVector::new(603894, 2160842),
            end_unit_vec: GenomeVector::new(610317, 2160842),
            center: GenomeVector::new_arbritary((610317 + 603894) as f64 / 2.0, 2160842),
            arc_length: (610317 - 603894) as f64,
        };
        let circ_loc_4 = CircularGenomeLocation {
            replicon: TIGR4.replicons.get("NC_003028.3").unwrap(),
            strand: StrandSense::Reverse,
            start_unit_vec: GenomeVector::new(1970326, 2160842),
            end_unit_vec: GenomeVector::new(1970400, 2160842),
            center: GenomeVector::new_arbritary((1970400 + 1970326) as f64 / 2.0, 2160842),
            arc_length: (1970400 - 1970326) as f64,
        };
        let actual_values = [circ_loc_1, circ_loc_2, circ_loc_3, circ_loc_4];

        for (left, right) in actual_values.iter().zip(test_locs.iter()) {
            left.test_compare(right);
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_operator_search() {

        // CopY Operators
        let (operator_seq, table) = build_cop_permutation_table();
        let operators = SequencePermutations::new("Cop Operator".to_string(), operator_seq, table);
        let imported_genomes = import_tigr4_lacto_genome();
        
        // TIGR4
        let TIGR4 = imported_genomes.0;
        let TIGR4_search = SearchGenome::new(&TIGR4, &operators, None);

        // T4 PRINTING
        println!("\nTIGR4");
        for (index, op) in TIGR4_search.operators.clone().unwrap().iter().enumerate() {
            println!("[{}]\t{:?}\t{}\t\t5'-{}-3'\t{}:{}", index+1, op.dimension, op.linear_location.strand, op.seq, op.linear_location.start_bound+1, op.linear_location.end_bound);
        }

        // TIGR4 Known Values
        let mut known_operators: Vec<Operator> = Vec::new();
        let dimensional_statuses =   [0,0,1,1,1,1,1,0,0,1,1,0,1,0,1];
        let strand_statuses =        [1,0,1,1,1,1,1,1,1,1,1,1,1,1,1];
        let sequences = [
            "GACTACACTCGTAAGT",
            "ATTGACAACCGTAATC",
            "AGTGACAACTGTCATT",
            "GATTACAGGTGTCAAT",
            "ATTGACAAATGTAGAT",
            "ATTGACAAATGTAGAT",
            "GATTACATCTGTCAGC",
            "GTCTACAGGCGTAGTT",
            "GGCTACAATCGTAATT",
            "AGTGACACCTGTCGTC",
            "ATTGACATTTGTAACT",
            "AATGACATCCGTCATT",
            "AACTACAGGTGTAAAC",
            "AACGACAAGCGTCACC",
            "GGCTACATCTGTAAAT",
        ];
        let starts = [
            70530,
            444016,
            549370,
            568177,
            691554,
            691594,
            775948,
            918522,
            1022871,
            1153103,
            1692084,
            1728688,
            1746226,
            1775201,
            2123135,
        ];
        
        for (((seq, dim), strand), start) in sequences.iter()
                              .zip(dimensional_statuses.iter())
                              .zip(strand_statuses.iter()) 
                              .zip(starts.iter()) {
            
            let new_region = GenomeRegion {
                replicon_accession: "NC_003028.3".to_string(),
                replicon_strand: match strand {
                    0 => StrandSense::Reverse,
                    _ => StrandSense::Forward,
                },
                start_index_ord: *start,
                end_index_ord: start + 15,
            };
            
            let new_op = Operator {
                linear_location: LinearGenomeLocation::new(&new_region, &TIGR4),
                seq: seq.to_string(),
                dimension: match dim {
                    0 => OperatorDimension::Single,
                    _ => OperatorDimension::Double,
                },
            };

            known_operators.push(new_op);
        }

        // TIGR4 Test Values
        let test_operators = match &TIGR4_search.operators {
            Some(operators) => {
                operators
            },
            None => panic!("No operators."),
        };
        assert_eq!(known_operators, *test_operators);


        // LACTO PRINTING
        println!("\nLACTO");
        let lacto = imported_genomes.1;
        let lacto_search = SearchGenome::new(&lacto, &operators, None);
        for (index, op) in lacto_search.operators.clone().unwrap().iter().enumerate() {
            println!("[{}]\t{:?}\t{}\t\t5'-{}-3'\t{}:{}", index+1, op.dimension, op.linear_location.strand, op.seq, op.linear_location.start_bound+1, op.linear_location.end_bound);
        }

        // LACTO Known Values
        let mut known_operators: Vec<Operator> = Vec::new();
        let dimensional_statuses =   [1,1,1,1,1,0,1,0,1,0,1,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0];
        let strand_statuses =        [1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0,1,1,0,0,1,1,1,0,1,0];
        let sequences = [
            "GTTTACAATTGTAAAC",
            "ATTTACACTTGTAAAT",
            "GTTTACATGTGTAAAT",
            "GGCGACAGGTGTAGGC",
            "AATGACAATTGTAGGT",
            "GCTGACAGCCGTCGTC",
            "ATTGACAATTGTCAGT",
            "GCTTACAAGCGTCAAT",
            "GTTTACACGTGTAAAC",
            "ATTTACAACCGTAAAC",
            "GTTTACAATTGTAAAC",
            "ATTTACAATTGTCAGC",
            "ACCTACAAACGTAGTT",
            "GTTTACATCTGTAGAT",
            "GCTGACAATCGTCATC",
            "GTTTACAAACGTAAAC",
            "AGTTACATCTGTCACT",
            "ACCGACAGCCGTCATC",
            "ACTGACAATCGTCGTT",
            "GTTTACATGCGTCAAT",
            "ACCGACACCTGTCAAT",
            "GTTTACATATGTAATC",
            "ATTTACACTTGTAAAC",
            "GCTTACATCCGTAGCC",
            "AATTACAACTGTAATT",
            "ACCTACATGCGTCAAT",
        ];
        let starts = [
            79266,
            205350,
            386028,
            391644,
            470963,
            631419,
            663534,
            773909,
            845055,
            845914,
            873353,
            1168433,
            1262872,
            1280110,
            1384570,
            1511051,
            1580466,
            1700449,
            1787812,
            1926251,
            1986580,
            1991994,
            1992021,
            2161546,
            2182862,
            2293848,
        ];
        
        for (((seq, dim), strand), start) in sequences.iter()
                                .zip(dimensional_statuses.iter())
                                .zip(strand_statuses.iter()) 
                                .zip(starts.iter()) {
            
            let new_region = GenomeRegion {
                replicon_accession: "NC_002662.1".to_string(),
                replicon_strand: match strand {
                    0 => StrandSense::Reverse,
                    _ => StrandSense::Forward,
                },
                start_index_ord: *start as usize,
                end_index_ord: start + 15,
            };
            
            let new_op = Operator {
                linear_location: LinearGenomeLocation::new(&new_region, &lacto),
                seq: seq.to_string(),
                dimension: match dim {
                    0 => OperatorDimension::Single,
                    _ => OperatorDimension::Double,
                },
            };

            known_operators.push(new_op);
        }


        // Lacto Test Values
        let test_operators = match &lacto_search.operators {
            Some(operators) => {
                operators
            },
            None => panic!("No operators."),
        };
        assert_eq!(known_operators, *test_operators);
    }

    #[test]
    #[allow(non_snake_case)]
    fn derive_linear_from_genome_region_then_derive_circular_from_linear() {

        // Import Genomes -- runtime-values are measured in the helper function itself
        let (TIGR4, _) = import_tigr4_lacto_genome(); 

        // Build genome region for BLAST hit
        let proto_test = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 2_012_333,
            end_index_ord: 2_010_867,
        };

        // Linear location for a BLAST hit
        let actual_hit_loc = LinearGenomeLocation {
            replicon: TIGR4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_bound: 2_012_332,
            end_bound: 2_010_867,
        };
        let test_hit_loc = LinearGenomeLocation::new(&proto_test, &TIGR4);
        assert_eq!(actual_hit_loc, test_hit_loc);

        // Build genome region for TEST GENE
        let proto_gene = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 1_922_125,
            end_index_ord: 1_922_997,
        };

        // Linear location for a GENE
        let actual_gene_loc = LinearGenomeLocation {
            replicon: TIGR4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_bound: 1_922_124,
            end_bound: 1_922_997,
        };
        let test_gene_loc = LinearGenomeLocation::new(&proto_gene, &TIGR4);
        assert_eq!(actual_gene_loc, test_gene_loc);


        // Calculate and compare circular locations for HIT + GENE
        let actual_hit_loc = CircularGenomeLocation {
            replicon: TIGR4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_unit_vec: GenomeVector::new(2_012_332, REPLICON_SIZE),
            end_unit_vec: GenomeVector::new(2_010_867, REPLICON_SIZE),
            center: GenomeVector::new_arbritary(931_178.50, REPLICON_SIZE),
            arc_length: 2159377.0,
        };
        let test_hit_circ_loc = CircularGenomeLocation::new(&test_hit_loc);
        println!("{}", actual_hit_loc);
        println!("{}", test_hit_circ_loc);
        actual_hit_loc.test_compare(&test_hit_circ_loc);

        const REPLICON_SIZE: usize = 2_160_842;
        let actual_gene_loc = CircularGenomeLocation {
            replicon: TIGR4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
            strand: StrandSense::Reverse,
            start_unit_vec: GenomeVector::new(1_922_124, REPLICON_SIZE),
            end_unit_vec: GenomeVector::new(1_922_997, REPLICON_SIZE),
            center: GenomeVector::new_arbritary(1_922_560.50, REPLICON_SIZE),
            arc_length: 873.0,
        };
        let test_gene_circ_loc = CircularGenomeLocation::new(&test_gene_loc);
        actual_gene_loc.test_compare(&test_gene_circ_loc);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_spatial_relationship_via_genome_regions() {

        // Import Genomes -- runtime-values are measured in the helper function itself
        let (TIGR4, _) = import_tigr4_lacto_genome(); 

        // Build genome region for BLAST hit
        let proto_hit_test = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 2_012_333,
            end_index_ord: 2_010_867,
        };

        // Linear location for BLAST hit
        let hit_loc = LinearGenomeLocation::new(&proto_hit_test, &TIGR4);
        println!("{}", hit_loc);

        // Build genome region for TEST GENE
        let proto_gene_test = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 1_922_125,
            end_index_ord: 1_922_997,
        };

        // Linear location for GENE
        let gene_loc = LinearGenomeLocation::new(&proto_gene_test, &TIGR4);
        println!("{}", gene_loc);

        // Calculate and compare spatial relationships
        let actual_relationship_1 = SpatialRelationship::Overlap(OverlapType::ContainerOf);
        let actual_relationship_2 = SpatialRelationship::Overlap(OverlapType::EngulfedBy);

        let test_relationship_1 = spatial_relationship(&gene_loc, &hit_loc);
        let test_relationship_2 = spatial_relationship(&hit_loc, &gene_loc);

        assert_eq!(actual_relationship_1, test_relationship_1);
        assert_eq!(actual_relationship_2, test_relationship_2);

        println!("Test Relationship of HIT relative to GENE: {:?}", test_relationship_1);
        println!("Test Relationship of GENE relative to HIT: {:?}", test_relationship_2);
        println!("Actual Relationship of HIT relative to GENE: {:?}", actual_relationship_1);
        println!("Actual Relationship of GENE relative to HIT: {:?}", actual_relationship_2);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_spatial_relationship_1() {

        // Import Genomes
        let (TIGR4, _) = import_tigr4_lacto_genome();

        // 5' Neighbor
        let gene_region_left = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Forward,
            start_index_ord: 1_078_629,
            end_index_ord: 1_082_279,
        };
        let gene_left = LinearGenomeLocation::new(&gene_region_left, &TIGR4);
        println!("\n5' Neighbor\n{}", gene_left);

        // 5' Gap
        let gap_region_1 = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Other,
            start_index_ord: 1_082_280,
            end_index_ord: 1_083_880,
        };
        let gap_1 = LinearGenomeLocation::new(&gap_region_1, &TIGR4);
        let gap_1_len = (gap_1.start_bound as i64 - gap_1.end_bound as i64).abs() as usize;
        println!("5' Gap\n{}", gap_1);

        // Origin
        let gene_region_center = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Forward,
            start_index_ord: 1_083_881,
            end_index_ord: 1_089_895,
        };
        let gene_center = LinearGenomeLocation::new(&gene_region_center, &TIGR4);
        println!("Origin\n{}", gene_center);

        // Left-Center Relationship
        let std_left = spatial_relationship(&gene_center, &gene_left);
        let inv_left = spatial_relationship(&gene_left, &gene_center);
        
        let std_left_len = match &std_left {
            SpatialRelationship::Neighbor(x) => match x {
                NeighborType::FivePrime(x) => *x,
                NeighborType::ThreePrime(x) => *x,
            },
            _ => panic!("ERROR: should be neighbors!"),
        };

        let inv_left_len = match &inv_left {
            SpatialRelationship::Neighbor(x) => match x {
                NeighborType::FivePrime(x) => *x,
                NeighborType::ThreePrime(x) => *x,
            },
            _ => panic!("ERROR: should be neighbors!"),
        };

        println!("std: {:?}\ninv: {:?}\n", std_left, inv_left);

        // 3' Gap
        let gap_region_2 = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Other,
            start_index_ord: 1_089_896,
            end_index_ord: 1_108_257,
        };
        let gap_2 = LinearGenomeLocation::new(&gap_region_2, &TIGR4);
        let gap_2_len = (gap_2.start_bound as i64 - gap_2.end_bound as i64).abs() as usize;
        println!("3' Gap\n{}", gap_2);

        // 3' Neighbor
        let gene_region_right = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 1_108_258,
            end_index_ord: 1_110_717,
        };
        let gene_right = LinearGenomeLocation::new(&gene_region_right, &TIGR4);
        println!("3' Neighbor\n{}", gene_right);

        // Right-Center Relationship
        let std_right = spatial_relationship(&gene_center, &gene_right);
        let inv_right = spatial_relationship(&gene_right, &gene_center);

        let std_right_len = match &std_right {
            SpatialRelationship::Neighbor(x) => match x {
                NeighborType::FivePrime(x) => *x,
                NeighborType::ThreePrime(x) => *x,
            },
            _ => panic!("ERROR: should be neighbors!"),
        };

        let inv_right_len = match &inv_right {
            SpatialRelationship::Neighbor(x) => match x {
                NeighborType::FivePrime(x) => *x,
                NeighborType::ThreePrime(x) => *x,
            },
            _ => panic!("ERROR: should be neighbors!"),
        };
        println!("std: {:?}\ninv: {:?}\n", std_right, inv_right);

        // Test for equivalent gap distance values
        assert_eq!(std_left_len, inv_left_len);
        assert_eq!(std_right_len, inv_right_len);
        assert_eq!(std_left_len, gap_1_len);
        assert_eq!(std_right_len, gap_2_len);
        
        // Test output relationship types
        assert_eq!(SpatialRelationship::Neighbor(NeighborType::FivePrime(1601)), std_left);
        assert_eq!(SpatialRelationship::Neighbor(NeighborType::ThreePrime(1601)), inv_left);
        assert_eq!(SpatialRelationship::Neighbor(NeighborType::ThreePrime(18362)), std_right);
        assert_eq!(SpatialRelationship::Neighbor(NeighborType::ThreePrime(18362)), inv_right);
    }

    use std::fs::File;
    use std::io::{self, BufRead};
    #[test]
    #[allow(non_snake_case)]
    fn test_spatial_relationship_2_genome_locations() {

        // Import manually calculated answer key
        let ans_key = PathBuf::from("tests/test_assets/T4_spatial_relationship_values.txt");
        let ans_file = File::open(ans_key).expect("ERROR: could not open 'T4_spatial_relationship_values.txt'");
        let answers = io::BufReader::new(ans_file).lines();
        let mut known_relationships: Vec<SpatialRelationship> = Vec::new();

        for line in answers {
            if line.as_ref().expect("Could not read line.").chars().rev().last() == Some('#') {
                continue;
            }

            if let Ok(data) = line {
                let tmp = data.split(' ').collect::<Vec<&str>>();

                let relationship = match tmp[0] {
                    "DR" => SpatialRelationship::None,
                    "5N" => SpatialRelationship::Neighbor(NeighborType::FivePrime(tmp[1].parse::<usize>().unwrap())),
                    "3N" => SpatialRelationship::Neighbor(NeighborType::ThreePrime(tmp[1].parse::<usize>().unwrap())),
                    "PE" => SpatialRelationship::Overlap(OverlapType::PerfectEclipse),
                    "5O" => SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary),
                    "3O" => SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary),
                    "EB" => SpatialRelationship::Overlap(OverlapType::EngulfedBy),
                    "CO" => SpatialRelationship::Overlap(OverlapType::ContainerOf),
                    "CR" => SpatialRelationship::Overlap(OverlapType::Crescent),
                    "CG" => SpatialRelationship::Overlap(OverlapType::CrescentGap),
                     _ => panic!("ERROR: ANSWER KEY CONTAINS INVALID SIGNIFIER."),
                };

                known_relationships.push(relationship);
            }
        }

        // CopY Operators
        let (operator_seq, table) = build_cop_permutation_table();
        let operators = SequencePermutations::new("Cop Operator".to_string(), operator_seq, table);

        // Import Genomes
        let (TIGR4, LACTO) = import_tigr4_lacto_genome();

        // Build genome region for TEST GENE + DISCONT. GENE + TEST BUBBLES

        // T4 Region over the discontinuity
        let t4_cross_discontinuity_genome_region = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 2_159_000,
            end_index_ord: 2295,
        };

        // Regular T4 Gene
        let t4_test_gene_genome_region = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Forward,
            start_index_ord: 1_083_881,
            end_index_ord: 1_089_895,
        };

        // Should be engulfing a gene + partially overlapping two others on each side
        let test_bubble_genome_region_1 = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 1_972_971,
            end_index_ord: 1_976_500,
        };

        // Should be completely engulfed by a gene
        let test_bubble_genome_region_2 = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 1_985_150,
            end_index_ord: 1_987_114,
        };

        // Linear location for GENE + Artificial Bubble
        let t4_gene_loc = LinearGenomeLocation::new(&t4_test_gene_genome_region, &TIGR4);
        let t4_cross_loc = LinearGenomeLocation::new(&t4_cross_discontinuity_genome_region, &TIGR4);
        let bubble_loc_1 = LinearGenomeLocation::new(&test_bubble_genome_region_1, &TIGR4);
        let bubble_loc_2 = LinearGenomeLocation::new(&test_bubble_genome_region_2, &TIGR4);

        // Pull a list of elements from genome annotation both near target gene itself and near its antipode
        let mut global_element_list: Vec<SearchGene> = Vec::new();

        let t4_test_genome = SearchGenome::new(&TIGR4, &operators, None);
        let lacto_test_genome = SearchGenome::new(&LACTO, &operators, None);

        let local_start: usize = 2182;
        let local_end: usize= 2237;
        let mut list_of_local_test_elements = match &t4_test_genome.genes {
            None => panic!("ERROR: no gene annotations were detected!"),
            Some(list) => Vec::from_iter(list[local_start..=local_end].iter().cloned())
        };
        
        let antipodal_start: usize = 1;
        let antipodal_end: usize= 24;
        let mut list_of_antipodal_test_elements = match &t4_test_genome.genes {
            None => panic!("ERROR: no gene annotations were detected!"),
            Some(list) => Vec::from_iter(list[antipodal_start..=antipodal_end].iter().cloned())
        };

        let discontinuity_start_1: usize = 4461;
        let discontinuity_end_1: usize = 4466;
        let discontinuity_start_2: usize = 0;
        let discontinuity_end_2: usize = 6;
        let list_of_cross_discontinuity_element = match &t4_test_genome.genes {
            None => panic!("ERROR: no gene annotations were detected!"),
            Some(list) => {
                let mut storage = Vec::new();
                storage.append(&mut Vec::from_iter(list[discontinuity_start_1..=discontinuity_end_1].iter().cloned()));
                storage.append(&mut Vec::from_iter(list[discontinuity_start_2..=discontinuity_end_2].iter().cloned()));

                storage
            }
        };

        let bubble_1_genes_start = 4094;
        let bubble_1_genes_end = 4107;
        let list_of_bubble_1_elements = match &t4_test_genome.genes {
            None => panic!("ERROR: no gene annotations were detected!"),
            Some(list) => Vec::from_iter(list[bubble_1_genes_start..=bubble_1_genes_end].iter().cloned())
        };

        let bubble_2_genes_start = 4122;
        let bubble_2_genes_end = 4123;
        let list_of_bubble_2_elements = match &t4_test_genome.genes {
            None => panic!("ERROR: no gene annotations were detected!"),
            Some(list) => Vec::from_iter(list[bubble_2_genes_start..=bubble_2_genes_end].iter().cloned())
        };

        let lacto_genes_start = 1569;
        let lacto_genes_end = 1573;
        let mut list_of_lacto_genes = match &lacto_test_genome.genes {
            None => panic!("ERROR: no gene annotations were detected!"),
            Some(list) => Vec::from_iter(list[lacto_genes_start..=lacto_genes_end].iter().cloned())
        };

        global_element_list.append(&mut list_of_local_test_elements);
        global_element_list.append(&mut list_of_antipodal_test_elements);
        global_element_list.append(&mut list_of_lacto_genes);

        for (gene, known_relationship) in global_element_list.iter().zip(known_relationships.iter()) {
            let relationship = spatial_relationship(&t4_gene_loc, gene.get_linear_location());
            println!("{:?}", relationship);
            assert_eq!(*known_relationship, relationship);
        }

        for (gene, known_relationship) in list_of_bubble_1_elements.iter().zip(known_relationships[85..=98].iter()) {
            let relationship = spatial_relationship(&bubble_loc_1, gene.get_linear_location());
            println!("{:?}", relationship);
            assert_eq!(*known_relationship, relationship);
        }

        for (gene, known_relationship) in list_of_bubble_2_elements.iter().zip(known_relationships[99..=100].iter()) {
            let relationship = spatial_relationship(&bubble_loc_2, gene.get_linear_location());
            println!("{:?}", relationship);
            assert_eq!(*known_relationship, relationship);
        }

        for (gene, known_relationship) in list_of_cross_discontinuity_element.iter().zip(known_relationships[101..=113].iter()) {
            let relationship = spatial_relationship(&t4_cross_loc, gene.get_linear_location());
            println!("{:?}\n", relationship);
            assert_eq!(*known_relationship, relationship);
        }

        // Test crescent gap
        let relationship =  spatial_relationship(&list_of_cross_discontinuity_element[6].linear_location, &t4_cross_loc);
        assert_eq!(known_relationships[114], relationship);

    }

    #[test]
    #[allow(non_snake_case)]
    fn timing_diagnostic_and_BLAST_linker_test() {

        // Load BLAST validation data
        let validation_path = PathBuf::from("tests/test_assets/blast_hit_validation_values.txt");
        let validation_file = File::open(validation_path).expect("ERROR: could not open 'blast_hit_validation_values.txt'");
        let lines = io::BufReader::new(validation_file).lines();
        let mut prevalidated_assocs: Vec<HashSet<BlastAssociationType>> = Vec::new();

        for line in lines {
            if let Ok(x) = line {
                let mut tmp_set: HashSet<BlastAssociationType> = HashSet::new();

                for association in x.split(' ') {
                    let tmp = match association {
                        "CopA" => BlastAssociationType::CopA,
                        "CupA" => BlastAssociationType::CupA,
                        "CopY" => BlastAssociationType::CopY,
                        "CopZ" => BlastAssociationType::CopZ,
                        _ => panic!("ERROR: Unrecognized BlastType."),
                    };
                    
                    tmp_set.insert(tmp);
                }

                prevalidated_assocs.push(tmp_set);
            }
        }
        
        // Load BLAST hits
        let now = SystemTime::now();
        let blast_files = [
            "tests/test_assets/CopA_blast_result.txt",
            "tests/test_assets/CupA_blast_result.txt",
            "tests/test_assets/CopY_blast_result.txt",
            "tests/test_assets/CopZ-TcrZ_blast_result.txt",
        ];

        let mut blast_files: Vec<PathBuf> = blast_files.into_iter().map(|x| PathBuf::from(x)).collect();
        let CopA_table = genome::BlastHitsTable::build_table_from_file(genome::BlastAssociationType::CopA, blast_files.remove(0));
        let CupA_table = genome::BlastHitsTable::build_table_from_file(genome::BlastAssociationType::CupA, blast_files.remove(0));
        let CopY_table = genome::BlastHitsTable::build_table_from_file(genome::BlastAssociationType::CopY, blast_files.remove(0));
        let CopZ_table = genome::BlastHitsTable::build_table_from_file(genome::BlastAssociationType::CopZ, blast_files.remove(0));
        let BLAST_tables = vec![CopA_table, CupA_table, CopY_table, CopZ_table];
        println!("BLAST results loading time: {}ms", now.elapsed().unwrap().as_millis());

        // CopY Operators
        let now = SystemTime::now();
        let (operator_seq, table) = build_cop_permutation_table();
        let operators = SequencePermutations::new("Cop Operator".to_string(), operator_seq, table);
        println!("Operator calculation time: {}ms", now.elapsed().unwrap().as_millis());

        // Import Genomes -- runtime-values are measured in the helper function itself
        let (TIGR4, LACTO) = import_tigr4_lacto_genome(); 

        // Parse + search genomes
        let t4_now = SystemTime::now();
        let TIGR4_search = SearchGenome::new(&TIGR4, &operators, Some(&BLAST_tables));
        println!("Build T4 search genome struct time: {}ms", t4_now.elapsed().unwrap().as_millis());

        let now = SystemTime::now();
        let LACTO_search = SearchGenome::new(&LACTO, &operators, Some(&BLAST_tables));
        println!("Build LACTO search genome struct time: {}ms", now.elapsed().unwrap().as_millis()); 

        // T4 Genome Blast Hit Validation
        for (gene, validated) in TIGR4_search.genes.unwrap()
                                                   .iter()
                                                   .filter(|&x| x.blast_association != None && x.gene.feature_type != FeatureType::Gene)
                                                   .zip(prevalidated_assocs[0..=5].iter()) {
            
            let test_assoc = gene.blast_association.as_ref().unwrap();

            println!();
            println!("Valid: {:?}", test_assoc);
            println!("Test: {:?}", validated);
            println!();
            
            assert_eq!(validated, test_assoc);
        }

        // LACTO Genome Blast Hit Validation
        for (gene, validated) in LACTO_search.genes.unwrap()
        .iter()
        .filter(|&x| x.blast_association != None && x.gene.feature_type != FeatureType::Gene)
        .zip(prevalidated_assocs[6..].iter()) {
            let test_assoc = gene.blast_association.as_ref().unwrap();

            println!();
            println!("Valid: {:?}", validated);
            println!("Test: {:?}", test_assoc);
            println!();

            assert_eq!(validated, test_assoc);
        }
        println!("\n-----\n");
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_search_bubble_generator() {

        // CopY Operators
        let (operator_seq, table) = build_cop_permutation_table();
        let operators = SequencePermutations::new("Cop Operator".to_string(), operator_seq, table);

        // Import Genomes -- runtime-values are measured in the helper function itself
        let (TIGR4, LACTO) = import_tigr4_lacto_genome(); 

        // Import Search Bubble validation data
        let validation_path = "tests/test_assets/search_bubble_validation.txt";
        let validation_file = File::open(validation_path).expect("ERROR: search bubble validation file is missing.");
        let lines = io::BufReader::new(validation_file).lines();
        let mut validated_vals: Vec<CircularGenomeLocation> = Vec::new();

        for line in lines {
            if let Ok(text_data) = line {
                let mut tmp = text_data.split(' ').map(|x| x.to_string()).collect::<Vec<String>>();
                let (genome, dir, left, right) = (tmp.remove(0), tmp.remove(0), tmp.remove(0).parse::<usize>().unwrap(), tmp.remove(0).parse::<usize>().unwrap());

                let replicon = match genome.as_str() {
                    "L" => LACTO.replicons.get(&"NC_002662.1".to_string()).unwrap(),
                    "S" => TIGR4.replicons.get(&"NC_003028.3".to_string()).unwrap(),
                     _  => panic!("ERROR: INVALID REPLICON IDENTIFIER."),
                };

                let strand = match dir.as_str() {
                    "F" => StrandSense::Forward,
                    "R" => StrandSense::Reverse,
                     _  => panic!("ERROR: INVALID STRAND IDENTIFIER."),
                };

                let loc_tmp = LinearGenomeLocation {
                    replicon,
                    strand,
                    start_bound: left,
                    end_bound: right,
                };

                let loc_tmp = CircularGenomeLocation::new(&loc_tmp);

                validated_vals.push(loc_tmp);

            } else {
                panic!("ERROR: could not parse data from search bubble validation file.")
            };
        }

        // Parse genomes
        let TIGR4_search = SearchGenome::new(&TIGR4, &operators, None);
        let LACTO_search = SearchGenome::new(&LACTO, &operators, None);

        // Pull specific genes from genomes
        let mut target_genes = TIGR4_search.genes.as_ref()
                                                 .unwrap()
                                                 .iter()
                                                 .take(5)
                                                 .filter(|x| x.gene.feature_type != FeatureType::Gene)
                                                 .map(|x| x.clone())
                                                 .collect::<Vec<SearchGene>>();
        let mut last_genes = TIGR4_search.genes.as_ref()
                                               .unwrap()
                                               .iter()
                                               .rev()
                                               .take(4)
                                               .filter(|x| x.gene.feature_type != FeatureType::Gene)
                                               .map(|x| x.clone())
                                               .collect::<Vec<SearchGene>>();
        target_genes.append(&mut last_genes);

        let mut first_lacto_genes = LACTO_search.genes.as_ref()
                                                      .unwrap()
                                                      .iter()
                                                      .map(|x| x.clone())
                                                      .collect::<Vec<SearchGene>>();
        target_genes.append(&mut first_lacto_genes);

        // Compare test search bubbles to actual search bubbles
        const T4_BUBBLE_RADIUS: usize = 2121;
        const LACTO_BUBBLE_RADIUS: usize = 100;

        println!("");
        for (index, (test, validated)) in target_genes.iter().zip(validated_vals[..5].iter()).enumerate() {
            let circ_bubble = CircularGenomeLocation::new(&search_bubble(test, T4_BUBBLE_RADIUS));

            println!("\n[{}]", index+1);
            println!("Test Subject\n{}", test.get_linear_location());
            println!("Test Bubble\n{}", circ_bubble);
            println!("Validated Bubble\n{}", validated);

            validated.test_compare(&circ_bubble);
        }

        println!("");
        for (index, (test, validated)) in target_genes[2692+5..=2698+6].iter().filter(|x| x.gene.feature_type == FeatureType::CDS).zip(validated_vals[5..].iter()).enumerate() {
            let circ_bubble = CircularGenomeLocation::new(&search_bubble(test, LACTO_BUBBLE_RADIUS));

            println!("\n[{}]", index+1);
            println!("Test Subject\n{}", test.get_linear_location());
            println!("Test Bubble\n{}", circ_bubble);
            println!("Validated Bubble\n{}", validated);

            validated.test_compare(&circ_bubble);
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn find_nearby_test() {

        // Import manually calculated answer key
        let ans_key = PathBuf::from("tests/test_assets/nearby_validation.txt");
        let ans_file = File::open(ans_key).expect("ERROR: could not open 'nearby_validation.txt'");
        let answers = io::BufReader::new(ans_file).lines();
        let mut known_relationships: Vec<SpatialRelationship> = Vec::new();

        for line in answers {
            if line.as_ref().expect("Could not read line.").chars().rev().last() == Some('#') {
                continue;
            }

            if let Ok(data) = line {
                let tmp = data.split(' ').collect::<Vec<&str>>();

                let relationship = match tmp[0] {
                    "DR" => SpatialRelationship::None,
                    "5N" => SpatialRelationship::Neighbor(NeighborType::FivePrime(tmp[1].parse::<usize>().unwrap())),
                    "3N" => SpatialRelationship::Neighbor(NeighborType::ThreePrime(tmp[1].parse::<usize>().unwrap())),
                    "PE" => SpatialRelationship::Overlap(OverlapType::PerfectEclipse),
                    "5O" => SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary),
                    "3O" => SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary),
                    "EB" => SpatialRelationship::Overlap(OverlapType::EngulfedBy),
                    "CO" => SpatialRelationship::Overlap(OverlapType::ContainerOf),
                    "CR" => SpatialRelationship::Overlap(OverlapType::Crescent),
                    "CG" => SpatialRelationship::Overlap(OverlapType::CrescentGap),
                        _ => panic!("ERROR: ANSWER KEY CONTAINS INVALID SIGNIFIER."),
                };

                known_relationships.push(relationship);
            }
        }

        // Import Genomes -- runtime-values are measured in the helper function itself
        let (TIGR4, LACTO) = import_tigr4_lacto_genome(); 

        // CopY Operators
        let (operator_seq, table) = build_cop_permutation_table();
        let operators = SequencePermutations::new("Cop Operator".to_string(), operator_seq, table);

        // Parse genomes
        let TIGR4_search = SearchGenome::new(&TIGR4, &operators, None);
        let T4_genome_objects = TIGR4_search.genes.as_ref()
                                                  .unwrap()
                                                  .iter()
                                                  // Adding filter here, but not on LACTO, since this one is closer to start;
                                                  // leave 'CDS' filter off of the LACTO iterator cause need to know item's location
                                                  // which is harder to do after filtering.
                                                  .filter(|&x| x.gene.feature_type == FeatureType::CDS)
                                                  .map(|x| GenomeObject::Gene(x.clone())).collect::<Vec<GenomeObject>>();
        let T4_genes = TIGR4_search.genes.unwrap()
                                         .iter()
                                         .filter(|&x| x.gene.feature_type == FeatureType::CDS)
                                         .map(|x| x.clone())
                                         .collect::<Vec<SearchGene>>();

        let LACTO_search = SearchGenome::new(&LACTO, &operators, None);
        let LACTO_genome_objects = LACTO_search.genes.as_ref()
                                                     .unwrap()
                                                     .iter()
                                                     .map(|x| GenomeObject::Gene(x.clone())).collect::<Vec<GenomeObject>>();
        let LACTO_genes = LACTO_search.genes.unwrap();

        // Verify Nearby T4 Elements
        let nearby_T4 = find_nearby_elements(&T4_genes[0], &T4_genome_objects, 4000).unwrap();
        for (item, known_relationship) in nearby_T4.iter().zip(known_relationships.iter()) {
            if let GenomeObject::Gene(_) = item.0 {
                assert_eq!(*known_relationship, item.1)
            }
        }

        // Verify Nearby LACTO Elements
        let nearby_LACTO = find_nearby_elements(&LACTO_genes[1041], &LACTO_genome_objects, 618).unwrap();
        for (item, known_relationship) in nearby_LACTO.iter().zip(known_relationships[11..].iter()) {
            if let GenomeObject::Gene(_) = item.0 {
                assert_eq!(*known_relationship, item.1);
            }
        }
    }

    // HELPER FUNCTIONS
    impl GenomeVector {
        
        // calculates a genome vector pointing to the center of two bounds 
        fn new_arbritary(bound: f64, replicon_size: usize) -> GenomeVector {

            let start = bound % replicon_size as f64;
            let theta = (2.0 * PI * start) / (replicon_size as f64);

            let x = theta.cos();
            let y = theta.sin();

            GenomeVector(x, y)
        }
    }

    impl<'a> CircularGenomeLocation<'a> {
        fn test_compare(&self, other: &Self) {

            const DELTA: f64 = 1e-9;

            if self.replicon != other.replicon {
                println!("Left: {}\nRight: {}\n", self.replicon.accession_id, other.replicon.accession_id);
                panic!("Genome Locations are not on the same replicon!\n");
            } else {
                println!("Replicons are identical.")
            }

            if self.strand != other.strand {
                println!("Left: {}\nRight: {}\n", self.strand, other.strand);
                panic!("Genome Locations are not on the same strand!")
            } else {
                println!("Strands are identical.")
            }

            if (&self.start_unit_vec - &other.start_unit_vec).magnitude() > DELTA {
                println!("Left: {}\nRight: {}\n", self.start_unit_vec, other.start_unit_vec);
                panic!("Genome Locations do not have the same START vector!")
            } else {
                println!("Start unit vectors are within a distance of {}", DELTA)
            }

            if (&self.end_unit_vec - &other.end_unit_vec).magnitude() > DELTA {
                println!("Left: {}\nRight: {}\n", self.end_unit_vec, other.end_unit_vec);
                panic!("Genome Locations do not have the same END vector!")
            } else {
                println!("Start end vectors are within a distance of {}", DELTA)
            }

            if (&self.center - &other.center).magnitude() > DELTA {
                println!("Left: {}\nRight: {}\n", self.center, other.center);
                panic!("Genome Locations do not have the same CENTER vector!")
            } else {
                println!("Center vectors are within a distance of {}", DELTA)
            }

            if (self.arc_length - other.arc_length).abs() > DELTA {
                println!("Left: {}\nRight: {}\n", self.arc_length, other.arc_length);
                panic!("Genome Locations do not have the same arc lengths!")
            } else {
                println!("Arc lengths differ by less than {}", DELTA)
            }
        }
    }

    use std::fmt::Display;
    impl Display for CircularGenomeLocation<'_> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

            let x_hat = GenomeVector(1.0, 0.0);
            let mut theta_1 = x_hat.angle(&self.start_unit_vec).to_degrees();
            let mut theta_2 = x_hat.angle(&self.end_unit_vec).to_degrees();
            let mut theta_3 = x_hat.angle(&self.center).to_degrees();
            let orientation_1 = x_hat.cross(&self.start_unit_vec);
            let orientation_2 = x_hat.cross(&self.end_unit_vec);
            let orientation_3 = x_hat.cross(&self.center);

            if orientation_1 < 0.0 {
                theta_1 = 360.0 - theta_1;
            }

            if orientation_2 < 0.0 {
                theta_2 = 360.0 - theta_2;
            }

            if orientation_3 < 0.0 {
                theta_3 = 360.0 - theta_3;
            }

            write!(f, "Circular Location\nReplicon: {}\nStrand: {:.15}\nReplicon Size: {}bp\nArc Length: {:.6}bp\nStart:\t{:.15} (θ₁ = {:.5}°)\n\
                       End:\t{:.15} (θ₂ = {:.5}°)\nCenter:\t{:.15} (θ₃ = {:.5}°)\n", 
                       self.replicon.accession_id, self.strand, self.replicon.len, self.arc_length, self.start_unit_vec,
                       theta_1, self.end_unit_vec, theta_2, self.center, theta_3)
        }    
    }

    impl Display for LinearGenomeLocation<'_> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let strand = &self.strand;
            let parent_id = &self.replicon.accession_id;
            let start = self.start_bound;
            let end = self.end_bound;
            let element_len = if end > start {
                (end as i64 - start as i64).abs()
            } else {
                (self.replicon.len as i64) - (end as i64 - start as i64).abs()
            };

            write!(f, "Linear Location\nReplicon: {}\nStrand: {}\nStart Bound: {}\nEnd Bound: {}\nFeature Length: {}\n",
                    parent_id, strand, start, end, element_len)
        }
    }

    impl Display for StrandSense {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            let txt = match self {
                StrandSense::Forward => "Forward",
                StrandSense::Reverse => "Reverse",
                StrandSense::Other => "Other",
            };
            write!(f, "{}", txt)
        }
    }
    
    impl Display for GenomeVector {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "<{}, {}>", self.0, self.1)
        }
    }
}