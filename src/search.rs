#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

use crate::genome::{self, Gene, Genome, GenomeRegion, StrandSense, BlastFragment, Replicon, GetSequence, ReverseComplement, BlastHitsTable};
use crate::permutations::SequencePermutations;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Index};
use std::collections::{HashSet, HashMap};
use std::time::{Duration, SystemTime};
use std::thread::sleep;

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
}

enum GenomeObject<'blast, 'genome> {
    Gene(SearchGene<'genome>),
    Operator(Operator<'genome>),
    Blast(SearchBlastFragment<'blast, 'genome>),
}

#[derive(Debug, Clone, PartialEq)]
enum PalindromeStatus {
    Yes,
    No,
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
        let input_element_len = (input.end_bound as i64 - input.start_bound as i64).abs();

        // Circular positioning logistics
        let start_unit_vec = GenomeVector::new(input.start_bound, replicon_len);
        let end_unit_vec = GenomeVector::new(input.end_bound, replicon_len);
        let arc_orientation = (replicon_len as f64 / 2.0) - (input_element_len as f64);
        let center_unit_vec: GenomeVector;
        let arc_length: f64;

        // The calculated arc-length and the exact center vector depend on the spatial orientation b/w the bounding start-
        // and end-vectors of a given element. To determine this relationship, we use the length of the feature derived from
        // its linear genome location; if the feature is longer than half the replicon, we use the larger angle between
        // its vectors; if the feature is shorter than half the replicon, we use the shorter angle between the features.
        // This approach is slightly better than using the cross_pdt as the determinant of orientation as it is agnostic
        // of which bound is listed as the start bound, and which bound is listed as the end bound. NCBI derived blast data,
        // for elements sitting on the reverse strand, often lists end_ord_indicies that are LESS than start_ord_indices.
        // In the cross product scheme, this meant to take the LARGER arc between the two vectors, even if the gene feature
        // being referenced was really along the smaller one--a result we obviously don't want. This new orientation scheme 
        // where start and end vectors are determined by the length of the feature addresses this shortcoming.
        if arc_orientation > 0.0 {
            center_unit_vec = (&start_unit_vec + &end_unit_vec).normalize();
            arc_length = start_unit_vec.angle(&end_unit_vec) * replicon_radius;
        } else if arc_orientation < 0.0 {
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
        let operators = find_genome_operators(&mut proto_search_genome, permutations);

        proto_search_genome
    }
}

struct SearchGene<'genome> {
    gene: &'genome Gene,
    linear_location: LinearGenomeLocation<'genome>,
    blast_association: Option<HashSet<genome::BlastAssociationType>>,
}

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
    palindromic: PalindromeStatus,
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

/* How object 'self' is related to a given origin; turns spatial relationship into a method call on objects
trait Relationship: LocateOnLinearGenome {
    fn relationship<T: LocateOnLinearGenome>(&self, origin: T) -> SpatialRelationship {
        let origin = origin.get_linear_location();
        let target = &self.get_linear_location();
        spatial_relationship(origin, target)
    }
}
impl<'genome> Relationship for Operator<'genome> {}
impl<'genome> Relationship for SearchGene<'genome> {}
impl<'blast, 'genome> Relationship for SearchBlastFragment<'blast, 'genome> {}
*/

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

            let find = match valid {
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

            let fwd_seq = proto_operator.0.clone();
            let rev_seq = proto_operator.0.chars().rev().collect::<String>();
            
            let new_operator = Operator {
                linear_location: LinearGenomeLocation::new(&tmp_genome_region, parent_genome.genome),
                seq: proto_operator.0,
                palindromic: PalindromeStatus::No,
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

    // Mark any operators in the list as 'palindromic' if their start index appears more than once
    let operator_start_indicies: Vec<usize> = operator_list.iter().map(|x| x.linear_location.start_bound).collect();
    for operator in operator_list.iter_mut() {
        let start_bound = operator.linear_location.start_bound;
        let count = operator_start_indicies.iter().filter(|&n| *n == start_bound).count();
        if count > 1 {
            operator.palindromic = PalindromeStatus::Yes
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
    // Storage logistics
    let output: SpatialRelationship;

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
    let different_replicons = (A.replicon.accession_id != B.replicon.accession_id);
    let no_overlap = Δ_arc_length >= (A_radius + B_radius);
    let engulfing_overlap = f64::max(A_radius, B_radius) > (Δ_arc_length + f64::min(A_radius, B_radius));

    // Apply rules of the algorithm
    let output = if different_replicons {
        SpatialRelationship::None
    } else if no_overlap {
        let arc_distance = Δ_arc_length - A_radius - B_radius;

        if orientation.is_sign_positive() { // AxB > 0 => center of B is to the LEFT of center of A (for fwd strand)
            SpatialRelationship::Neighbor(NeighborType::ThreePrime(arc_distance as usize))
        } else { // AxB < 0 => center of B is to the RIGHT of center of A
            SpatialRelationship::Neighbor(NeighborType::FivePrime(arc_distance as usize))
        }
    } else if engulfing_overlap {
        if A_radius > B_radius {
            SpatialRelationship::Overlap(OverlapType::EngulfedBy)
        } else {
            SpatialRelationship::Overlap(OverlapType::ContainerOf)
        }
    } else {
        if orientation.is_sign_positive() {
            SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary)
        } else {
            SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary)
        }
    };

    output
}

fn search_bubble<'genome, T: LocateOnLinearGenome>(query: &'genome T, search_radius: usize) -> LinearGenomeLocation {
    let mut proto_bubble = query.get_linear_location().clone();
    let start_i64: i64 = (proto_bubble.start_bound as i64) - (search_radius as i64);
   
    proto_bubble.start_bound = start_i64.rem_euclid(proto_bubble.replicon.len as i64) as usize;
    proto_bubble.end_bound += search_radius;

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
        let gene_loc = &gene.linear_location;

        for hit in blast_hits {
            let hit_loc = &hit.linear_location;
            
            // If a new match is found
            if let SpatialRelationship::Overlap(_) = spatial_relationship(gene_loc, hit_loc) {
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
    use crate::genome::{self, GenomeRegion, FeatureType};
    use crate::import::{self, ProtoGenome, AssemblyMetadata};
    use std::path::PathBuf;
    use std::{collections::HashMap, f64::consts::E};

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

        for (index, (left, right)) in matrix_transformation_calculated_cases.zip(test_cases).enumerate() {
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

        // Load BLAST hits
        let blast_file = PathBuf::from("tests/test_assets/CopY_blast_result.txt");
        let blast_table = genome::BlastHitsTable::build_table_from_file(genome::BlastAssociationType::CopY, blast_file);
        
        // TIGR4
        let TIGR4 = imported_genomes.0;
        let TIGR4_search = SearchGenome::new(&TIGR4, &operators, None);

        // T4 PRINTING
        for (index, op) in TIGR4_search.operators.clone().unwrap().iter().enumerate() {
            println!("[{}]\t{:?}\t{}\t\t5'-{}-3'\t{}:{}", index, op.palindromic, op.linear_location.strand, op.seq, op.linear_location.start_bound+1, op.linear_location.end_bound);
        }

        // TIGR4 Known Values
        let mut known_operators: Vec<Operator> = Vec::new();
        let palindrome_statuses =   [0,0,1,1,1,1,1,0,0,1,1,0,1,0,1];
        let strand_statuses =       [1,0,1,1,1,1,1,1,1,1,1,1,1,1,1];
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
        
        for (((seq, palin), strand), start) in sequences.iter()
                              .zip(palindrome_statuses.iter())
                              .zip(strand_statuses.iter()) 
                              .zip(starts.iter()) {
            
            let new_region = GenomeRegion {
                replicon_accession: "NC_003028.3".to_string(),
                replicon_strand: match strand {
                    0 => StrandSense::Reverse,
                    _ => StrandSense::Forward,
                },
                start_index_ord: *start as usize,
                end_index_ord: start + 15,
            };
            
            let new_op = Operator {
                linear_location: LinearGenomeLocation::new(&new_region, &TIGR4),
                seq: seq.to_string(),
                palindromic: match palin {
                    0 => PalindromeStatus::No,
                    _ => PalindromeStatus::Yes,
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
        println!();
        let lacto = imported_genomes.1;
        let lacto_search = SearchGenome::new(&lacto, &operators, None);
        for (index, op) in lacto_search.operators.clone().unwrap().iter().enumerate() {
            println!("[{}]\t{:?}\t{}\t\t5'-{}-3'\t{}:{}", index+1, op.palindromic, op.linear_location.strand, op.seq, op.linear_location.start_bound+1, op.linear_location.end_bound);
        }

        // LACTO Known Values
        let mut known_operators: Vec<Operator> = Vec::new();
        let palindrome_statuses =   [1,1,1,1,1,0,1,0,1,0,1,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0];
        let strand_statuses =       [1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0,1,1,0,0,1,1,1,0,1,0];
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
        
        for (((seq, palin), strand), start) in sequences.iter()
                                .zip(palindrome_statuses.iter())
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
                palindromic: match palin {
                    0 => PalindromeStatus::No,
                    _ => PalindromeStatus::Yes,
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
            center: GenomeVector::new_arbritary(2_011_599.50, REPLICON_SIZE),
            arc_length: 1465.0,
        };
        let test_hit_circ_loc = CircularGenomeLocation::new(&test_hit_loc);
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
    fn test_overlap_1() {

        // Import Genomes -- runtime-values are measured in the helper function itself
        let (TIGR4, _) = import_tigr4_lacto_genome(); 

        // Build genome region for BLAST hit
        let proto_hit_test = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 2012333,
            end_index_ord: 2010867,
        };

        // Linear location for BLAST hit
        let hit_loc = LinearGenomeLocation::new(&proto_hit_test, &TIGR4);

        // Build genome region for TEST GENE
        let proto_gene_test = GenomeRegion {
            replicon_accession: "NC_003028.3".to_string(),
            replicon_strand: StrandSense::Reverse,
            start_index_ord: 1922125,
            end_index_ord: 1922997,
        };

        // Linear location for GENE
        let gene_loc = LinearGenomeLocation::new(&proto_gene_test, &TIGR4);

        // Calculate and compare spatial relationships
        let test_relationship_1 = spatial_relationship(&gene_loc, &hit_loc);
        let test_relationship_2 = spatial_relationship(&hit_loc, &gene_loc);

        let actual_relationship_1 = SpatialRelationship::Neighbor(NeighborType::FivePrime(87_869));
        let actual_relationship_2 = SpatialRelationship::Neighbor(NeighborType::ThreePrime(87_869));

        println!("Test Relationship of HIT relative to GENE: {:?}", test_relationship_1);
        println!("Test Relationship of GENE relative to HIT: {:?}", test_relationship_2);
        println!("Actual Relationship of HIT relative to GENE: {:?}", actual_relationship_1);
        println!("Actual Relationship of GENE relative to HIT: {:?}", actual_relationship_2);
    }

    






    use std::time::Duration; 
    use std::time::SystemTime;
    use std::thread;

    #[test]
    #[allow(non_snake_case)]
    fn timing_diagnostic_and_BLAST_linker_test() {
        
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

        // 
        for gene in TIGR4_search.genes.unwrap() {

            if gene.gene.feature_type == FeatureType::Gene {
                continue;
            }

            if let Some(x) = gene.blast_association {
                println!("{}\t{:?}\t{:?}", gene.linear_location.start_bound+1, gene.gene.feature_type, x);
            }
        }

        println!("\n-----\n");

        // 
        for gene in LACTO_search.genes.unwrap() {

            if gene.gene.feature_type == FeatureType::Gene {
                continue;
            }

            if let Some(x) = gene.blast_association {
                println!("{}\t{:?}\t{:?}", gene.linear_location.start_bound+1, gene.gene.feature_type, x);
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

            if (self.arc_length - other.arc_length) > DELTA {
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
            let element_len = (end as i64 - start as i64).abs();

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