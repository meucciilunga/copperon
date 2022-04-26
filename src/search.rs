#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

use crate::genome::{self, Gene, Genome, Operator, GenomeRegion, StrandSense, BlastFragment};
use crate::permutations::SequencePermutations;

// Special method that pulls a Gene/Operator/BlastFragment interval boundaries
// and center position; implemented via a default implementation of a trait so
// the GenomeLocate method can be managed in a single place
trait GenomeLocate: GetGenomeLocation {

    fn locate(&self) -> GenomeLocation {

        fn bounds(input: &GenomeRegion) -> (f64, f64) {
            let left = (input.start_index_ord - 1) as f64;
            let right = input.end_index_ord as f64;
            (left, right)
        }

        fn radius(input: &GenomeRegion) -> f64 {
            let (x1, x2) = bounds(input);
            (x2 - x1) / 2_f64
        }

        fn center(input: &GenomeRegion) -> f64 {
            let (x1, x2) = bounds(input);
            (x1 + x2) / 2_f64
        }

        let (start, end) = bounds(&self.get_location());
        let radius = radius(&self.get_location());
        let center = center(&self.get_location());

        GenomeLocation {
            id: self.get_location().replicon_accession.clone(),
            strand: self.get_location().replicon_strand.clone(),
            start_bound: start,
            end_bound: end,
            center,
            radius,
        }
    }

    fn locate_rev(&self, replicon_len: usize) -> GenomeLocation {

        fn rev_bounds(input: &GenomeRegion, replicon_len: usize) -> (f64, f64) {
            let offset = replicon_len + 1;
            let ord_start = offset - input.end_index_ord;
            let ord_end =   offset - input.start_index_ord;

            let left = (ord_start - 1) as f64;
            let right = ord_end as f64;
            (left, right)
        }

        fn rev_radius(input: &GenomeRegion, replicon_len: usize) -> f64 {
            let (x1, x2) = rev_bounds(input, replicon_len);
            (x2 - x1) / 2_f64
        }

        fn rev_center(input: &GenomeRegion, replicon_len: usize) -> f64 {
            let (x1, x2) = rev_bounds(input, replicon_len);
            (x1 + x2) / 2_f64
        }

        let (start, end) = rev_bounds(&self.get_location(), replicon_len);
        let radius = rev_radius(&self.get_location(), replicon_len);
        let center = rev_center(&self.get_location(), replicon_len);

        GenomeLocation {
            id: self.get_location().replicon_accession.clone(),
            strand: match &self.get_location().replicon_strand {
                StrandSense::Forward => StrandSense::Reverse,
                StrandSense::Reverse => StrandSense::Forward,
                StrandSense::Other => StrandSense::Other,
            },
            start_bound: start,
            end_bound: end,
            center,
            radius,
        }
    }
}

impl GenomeLocate for Gene {}
impl GenomeLocate for Operator {}
impl GenomeLocate for BlastFragment {}

trait GetGenomeLocation {
    fn get_location(&self) -> GenomeRegion;
}

impl GetGenomeLocation for Gene {
    fn get_location(&self) -> GenomeRegion {
        self.location.clone()
    }
}

impl GetGenomeLocation for Operator {
    fn get_location(&self) -> GenomeRegion {
        self.location.clone()
    }
}

impl GetGenomeLocation for BlastFragment {
    fn get_location(&self) -> GenomeRegion {
        self.location.clone()
    }
}

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
    Inside,
    Around,
}

#[derive(Debug, Clone, PartialEq)]
struct GenomeLocation {
    id: String,
    strand: StrandSense,
    start_bound: f64,
    end_bound: f64,
    center: f64,
    radius: f64,
}

// Given a gene and a BLAST table, update the gene to reflect whether it is associated
// with any of the hits in the BLAST Table
fn blast_link(input_gene: &mut Gene, blast_results: &genome::BlastHitsTable) {
    unimplemented!()
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
      S: GenomeLocate
{
    unimplemented!()
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

    // Define genome locations of elements
    let A = element_A.locate();
    let B = element_B.locate();

    // Define position vector b/w two elements
    let position = B.center - A.center;

    // Define rules of the algorithm
    let different_replicons = (A.id != B.id);
    let no_overlap = position.abs() >= (A.radius + B.radius);
    let engulfing_overlap = f64::max(A.radius, B.radius) > (position.abs() + f64::min(A.radius, B.radius));

    // Apply rules of algorithm
    let output = if different_replicons {
        SpatialRelationship::None
    } else if no_overlap {
        let distance = position.abs() - A.radius - B.radius;
        if position.is_sign_positive() {
            SpatialRelationship::Neighbor(NeighborType::ThreePrime(distance as usize))
        } else {
            SpatialRelationship::Neighbor(NeighborType::FivePrime(distance as usize))
        }
    } else if engulfing_overlap {
        if A.radius > B.radius {
            SpatialRelationship::Overlap(OverlapType::Inside)
        } else {
            SpatialRelationship::Overlap(OverlapType::Around)
        }
    } else {
        if position.is_sign_positive() {
            SpatialRelationship::Overlap(OverlapType::ThreePrimeBoundary)
        } else {
            SpatialRelationship::Overlap(OverlapType::FivePrimeBoundary)
        }
    };

    output
}

#[cfg(test)]
mod tests {
    use crate::genome::{self, GenomeRegion};
    use std::collections::HashMap;

    use super::{GenomeLocate, GenomeLocation};

    #[test]
    fn test_locate_boundaries_gene() {

        let mut attributes: HashMap<String, String> = HashMap::new();
        attributes.insert("test_key".to_string(), "test_val".to_string());

        let new_gene = genome::Gene {
            location: GenomeRegion {
                replicon_accession: "Alphabet".to_string(),
                replicon_strand: genome::StrandSense::Forward,
                start_index_ord: 7,
                end_index_ord: 52,
            },
            feature_type: genome::FeatureType::Other,
            attributes,
            blast_association: None,
        };

        let actual_fwd_location = GenomeLocation {
            id: "Alphabet".to_string(),
            strand: genome::StrandSense::Forward,
            start_bound: 6_f64,
            end_bound: 52_f64,
            center: 29_f64,
            radius: 23_f64,
        };

        let actual_rev_location = GenomeLocation {
            id: "Alphabet".to_string(),
            strand: genome::StrandSense::Reverse,
            start_bound: 48_f64,
            end_bound: 94_f64,
            center: 71_f64,
            radius: 23_f64,
        };

        let test_fwd_location = new_gene.locate();
        let test_rev_location = new_gene.locate_rev(100);

        assert_eq!(actual_fwd_location, test_fwd_location);
        assert_eq!(actual_rev_location, test_rev_location);
    }

    #[test]
    fn test_locate_boundaries_operator() {
        let test_operator = genome::Operator {
            location: genome::GenomeRegion {
                replicon_accession: "TEST_OPERATOR".to_string(),
                replicon_strand: genome::StrandSense::Forward,
                start_index_ord: 2_135_413,
                end_index_ord: 2_135_433,
            },
            seq: "AAAAATTTTTCCCCCGGGGG".to_string(),
            potential_operon: None,
        };

        let actual_fwd_location = GenomeLocation {
            id: "TEST_OPERATOR".to_string(),
            strand: genome::StrandSense::Forward,
            start_bound: 2_135_412_f64,
            end_bound: 2_135_433_f64,
            center: 2_135_422.50_f64,
            radius: 10.5_f64,
        };

        let actual_rev_location = GenomeLocation {
            id: "TEST_OPERATOR".to_string(),
            strand: genome::StrandSense::Reverse,
            start_bound: 1_424_881_f64,
            end_bound: 1_424_902_f64,
            center: 1_424_891.50_f64,
            radius: 10.5_f64,
        };

        let test_fwd_location = test_operator.locate();
        let test_rev_location = test_operator.locate_rev(3_560_314);

        assert_eq!(actual_fwd_location, test_fwd_location);
        assert_eq!(actual_rev_location, test_rev_location);
    }

    #[test]
    fn test_locate_boundaries_blast_fragment() {
        let test_operator = genome::BlastFragment {
            location: genome::GenomeRegion {
                replicon_accession: "TEST_FRAGMENT".to_string(),
                replicon_strand: genome::StrandSense::Forward,
                start_index_ord: 7_000_001,
                end_index_ord: 8_000_000,
            },
            length: 1_000_000,
            pident: 100.00,
            evalue: 19e-10,
            qseqid: "TEST_QUERY".to_string(),
            qstart: 1,
            qend: 1_000_000
        };

        let actual_fwd_location = GenomeLocation {
            id: "TEST_FRAGMENT".to_string(),
            strand: genome::StrandSense::Forward,
            start_bound: 7_000_000_f64,
            end_bound: 8_000_000_f64,
            center: 7_500_000_f64,
            radius: 500_000_f64,
        };

        let actual_rev_location = GenomeLocation {
            id: "TEST_FRAGMENT".to_string(),
            strand: genome::StrandSense::Reverse,
            start_bound: 11_763_508_f64,
            end_bound: 12_763_508_f64,
            center: 12_263_508_f64,
            radius: 500_000_f64,
        };

        let test_fwd_location = test_operator.locate();
        let test_rev_location = test_operator.locate_rev(19_763_508);

        assert_eq!(actual_fwd_location, test_fwd_location);
        assert_eq!(actual_rev_location, test_rev_location);
    }
}