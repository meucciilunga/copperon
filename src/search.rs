#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

use crate::genome::{self, Gene, Genome, Operator, GenomeRegion, StrandSense, BlastFragment, Replicon};
use crate::permutations::SequencePermutations;
use std::f64::consts::PI;
use std::ops::{Add, Sub};

enum SpatialRelationship {
    Neighbor(NeighborType),
    Overlap(OverlapType),
    None,
}

enum NeighborType {
    FivePrime(usize),
    ThreePrime(usize),
}

enum OverlapType {
    FivePrimeBoundary,
    ThreePrimeBoundary,
    EngulfedBy,
    ContainerOf,
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

    // Converts a genome position in linear space to circular space, based on the size of the circular genome / replicon
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

        // Only normalize if vectors aren't already unit vectors; otherwise,
        // normalization would just add unncessary error; this was revealed
        // during unit testing, as the four genome_arc_delta_tests, which have 
        // passing values that were set during my initial round of unit testing
        // with them, failed if an already normalized vector was normalized again.

        let a_diff = (self.magnitude() - 1.0).abs();
        let b_diff = (other.magnitude() - 1.0).abs();
        let delta = 1e-12;

        let a = if a_diff > delta   {
            self.normalize()
        } else {
            self.clone()
        };

        let b = if b_diff > delta   {
            other.normalize()
        } else {
            other.clone()
        };
        
        let theta = a.dot(&b);
        theta.acos()
    }   
}

// Defines a location in a linear genome by its BOUNDARY INDICIES
// NOT BY ITS ORDINAL INDICIES; BOUNDARY INDICIES ∈ ℝ (approximated by f64)
#[derive(Debug, Clone, PartialEq)]
struct LinearGenomeLocation {
    id: String,
    strand: StrandSense,
    start_bound: usize,
    end_bound: usize,
}

impl LinearGenomeLocation {
    fn new(input: &GenomeRegion) -> LinearGenomeLocation {
        let id = input.replicon_accession.clone();
        let strand = input.replicon_strand.clone();
        let left = input.start_index_ord - 1;
        let right = input.end_index_ord;

        LinearGenomeLocation {
            id,
            strand,
            start_bound: left,
            end_bound: right,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct CircularGenomeLocation<'a> {
    replicon: &'a Replicon,
    strand: StrandSense,
    start_unit_vec: GenomeVector,
    end_unit_vec: GenomeVector,
    position: GenomeVector,
}

impl<'a> CircularGenomeLocation<'a> {
    fn new(input: &LinearGenomeLocation, genome: &'a Genome) -> CircularGenomeLocation<'a> {

        let replicon_err = "ERROR: could not find GenomeRegion's replicon in specified Genome!";
        let replicon = genome.replicons.get(&input.id).expect(replicon_err);
        let replicon_len = replicon.fwd_strand.len();
        let replicon_radius = (replicon_len as f64) / (2.0 * PI);

        let strand = input.strand.clone();
        
        let start_unit_vec = GenomeVector::new(input.start_bound, replicon_len);
        let end_unit_vec = GenomeVector::new(input.end_bound, replicon_len);
        let center_unit_vec = (&start_unit_vec + &end_unit_vec).normalize();
        let mut arc_length = start_unit_vec.angle(&end_unit_vec) * replicon_radius;

        if arc_length.abs() < 0.25 {
            arc_length = 1.0;
        }

        let position = center_unit_vec.scalar(arc_length);

        CircularGenomeLocation {
            replicon,
            strand,
            start_unit_vec,
            end_unit_vec,
            position,
        }
    }
}

// Determine what spatial relationship two genomic elements A and B have with each other;
// this function is non-commutative, as the spatial relationship returned is always given
// as the positioning of B relative to A
#[allow(non_snake_case)]
#[allow(unused_parens)]
fn spatial_relationship<A,B>(element_A: &A, element_B: &B) -> SpatialRelationship 
where A: GenomeLocate,
      B: GenomeLocate
{   
    // Storage logistics
    let output: SpatialRelationship;

    // Define circular genome location of each element
    let A = element_A.locate();
    let B = element_B.locate();
    let replicon_radius = (A.replicon.len as f64) / (2.0 * PI);
    // NOTE: replicon.len == circumference of the genome; 2*pi*radius = circumference

    // Determine if B is to the 'left' or 'right' of A
    // For an element A on the forward strand:
    // |AxB| > 0 => B is a 3' neighbor
    // |AxB| < 0 => B is a 5' neighbor
    // For Fwd strand: CCW orientation == 5' -> 3'
    // For Fwd strand: CCW orientation == 3' -> 5'
    let A_hat = A.position.normalize();
    let B_hat = B.position.normalize();
    let orientation = match A.strand {
        StrandSense::Forward | StrandSense::Other => A_hat.cross(&B_hat),
        StrandSense::Reverse => -A_hat.cross(&B_hat),
    }; // Invert orientation result based on which strand element A is on;
    // This inversion of the sign of the orientation X-pdt depending on  
    // which strand a given gene is on allows the spatial relationship
    // algorithm to process relationship results for both Fwd/Rev strand
    // genes without having to duplicate very similar code with only minor
    // changes accounting for strand-specific Left/Right determination

    let A_radius = A.position.magnitude();
    let B_radius = B.position.magnitude();
    
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
            SpatialRelationship::Overlap(OverlapType::ContainerOf)
        } else {
            SpatialRelationship::Overlap(OverlapType::EngulfedBy)
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

// Given a genome and all the sequence permutations of an (operator) consensus sequence,
// locate all places in the genome where some variant of the consensus sequence is found.
fn find_genome_operators(target_genome: &Genome, query: &SequencePermutations) -> Vec<Operator> {
    unimplemented!()
}

// Given a genomic element and a list of all other elements (to be considered) on that genome,
// return a list of all genome elements within a distance 'search_radius' of the query's center
fn find_nearby_elements<Q,S>(query: &Q, possible_elements: &S, search_radius: usize) -> Vec<S>
where Q: GenomeLocate,
      S: GenomeLocate + SearchBubble
{
    unimplemented!()
}

// Given a gene and a BLAST table, update the gene to reflect whether it is associated
// with any of the hits in the BLAST Table
fn blast_link(input_gene: &mut Gene, blast_results: &genome::BlastHitsTable) {
    unimplemented!()
}

trait GetLocation {
    fn get_circ_location(&self) -> LinearGenomeLocation;
}

trait GenomeLocate: GetLocation {
    fn locate(&self) -> CircularGenomeLocation;
    fn locate_rev(&self) -> CircularGenomeLocation;
}
trait SearchBubble: GenomeLocate {}













// UNIT TESTS
#[cfg(test)]
mod tests {
    use super::*;
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

        let actual_cases = (0..=MAX_INDEX).map(|x| {
                                        let theta = 2.0*PI*(x as f64 / REPLICON_SIZE as f64);
                                        let x_hat = GenomeVector(1.0, 0.0);
                                        rotate_vector(theta, x_hat)
                                    });

        for (index, (left, right)) in actual_cases.zip(test_cases).enumerate() {
            let difference = (&left - &right).magnitude();
            let delta = 1e-4;
            assert!(delta > difference);
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn genome_arc_delta_test_1bp_arc_100Mbp_genome() {

        const REPLICON_SIZE: usize = 100_000_000;

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
        assert!(max_deviation < 0.1);
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
        assert!(max_deviation < 5e-3);
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
        let arc_length_samples = [1, 2, 3, 10, 25];

        let mut list_of_max_deviations: Vec<f64> = Vec::new();
        for arc_len in arc_length_samples  {
            
            // Calculate varying length arcs (1bp) for every base-pair position in the genome
            let mut list_of_arc_lengths: Vec<f64> = Vec::new();
            let mut list_of_indicies: Vec<usize> = Vec::new();

            for i in 0..=REPLICON_SIZE+1 {
                let lower_bound = GenomeVector::new(i, REPLICON_SIZE);
                let upper_bound = GenomeVector::new(i+1, REPLICON_SIZE);
                let seperating_angle = lower_bound.angle(&upper_bound);
                let radius = (REPLICON_SIZE as f64) / (2.0 * PI);
                let arc_length = seperating_angle * radius;

                list_of_arc_lengths.push(arc_length);
                list_of_indicies.push(i);
            }

            let delta_list = list_of_arc_lengths.iter().map(|x| (x - 1.0).abs());
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

        println!("MAXIMUM deviation for 1bp arc in {}bp genome across various sequence lengths {:?} = {}", REPLICON_SIZE, arc_length_samples, max_deviation);
        assert!(max_deviation < 0.1);
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
        assert!(max_deviation < 5e-3);
    }

    // Several of the following tests require a test Genome to work with, so this performs that import operation
    fn import_tigr4_genome() -> Genome {

        // Paths to relevant files
        let assembly_name = "GCF_000006885.1_ASM688v1".to_string();
        let genome_seq_file = PathBuf::from("tests/test_assets/GCF_000006885.1_ASM688v1/GCF_000006885.1_ASM688v1_genomic.fna");
        let genome_annotation_file = PathBuf::from("tests/test_assets/GCF_000006885.1_ASM688v1/GCF_000006885.1_ASM688v1_genomic.gff");
        let refseq_metadata_file = PathBuf::from("tests/test_assets/assembly_summary_refseq.txt");
        let genbank_metadata_file = PathBuf::from("tests/test_assets/assembly_summary_genbank.txt");

        // Import data
        let protogenome = import::parse_genome_sequence(&assembly_name, genome_seq_file);
        let protogenes = import::parse_genome_annotation(genome_annotation_file);
        let refseq_meta = import::import_database_summary(refseq_metadata_file).into_iter();
        let genbank_meta = import::import_database_summary(genbank_metadata_file).into_iter();
        let full_metadata: HashMap<String, AssemblyMetadata> = refseq_meta.chain(genbank_meta).collect();

        // Build genome
        let asm_pull_error = "ERROR: could not find assembly name in metadata database!";
        let metadata = full_metadata.get(&assembly_name).expect(asm_pull_error);
        let genes = Some(genome::Gene::convert_annotation_entry_list(&protogenes));
        
        genome::Genome::from_proto(protogenome, metadata, genes)
    }

    #[test]
    #[allow(non_snake_case)]
    fn derive_circular_genome_from_linear_1() {

        let TIGR4 = import_tigr4_genome();
        let mut test_genes = TIGR4.genes.clone().unwrap();
        let test_genes = [test_genes.remove(0), test_genes.remove(2), test_genes.remove(1212), test_genes.remove(4078)];

        for test_gene in &test_genes {
            let linear_location = LinearGenomeLocation::new(&test_gene.location);
            let circ_location = CircularGenomeLocation::new(&linear_location, &TIGR4);
            
            println!("{:?}", linear_location);
            println!("{}\n", circ_location);
        }

    }

    use std::fmt::Display;
    impl Display for CircularGenomeLocation<'_> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

            let theta_1 = GenomeVector(1.0, 0.0).angle(&self.start_unit_vec).to_degrees();
            let theta_2 = GenomeVector(1.0, 0.0).angle(&self.end_unit_vec).to_degrees();
            let theta_3 = GenomeVector(1.0, 0.0).angle(&self.position).to_degrees();

            write!(f, "Strand: {:.15}\nReplicon Size: {}bp\nArc Length: {:.2}bp\nStart:\t{:.15} (θ₁ = {:.5}°)\n\
                       End:\t{:.15} (θ₂ = {:.5}°)\nCenter:\t{:.15} (θ₃ = {:.5}°)", 
                       self.strand, self.replicon.len, self.position.magnitude(), self.start_unit_vec,
                       theta_1, self.end_unit_vec, theta_2, self.position, theta_3)
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