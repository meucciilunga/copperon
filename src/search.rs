#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]

use crate::genome::{Genome, Operator, GenomeLocation};
use crate::permutations::SequencePermutations;

fn find_genome_operators(target_genome: Genome, query: SequencePermutations) -> Vec<Operator> {
    unimplemented!()
}

// Given a genomic element and a list of all other elements (to be considered) on that genome,
// return a list of all genome elements within a distance 'search_radius' of the query's center
fn find_nearby_elements<Q,S>(query: Q, possible_elements: S, search_radius: usize) -> Vec<S>
where Q: GenomeLocate,
      S: GenomeLocate
{
    unimplemented!()
}

// Determine what spatial relationship two genomic elements A and B have with each other;
// this function is non-commutative, as the spatial relationship returned is always given
// as the positioning of B relative to A
#[allow(non_snake_case)]
fn spatial_relationship<A,B>(element_A: A, element_B: B) -> SpatialRelationship 
where A: GenomeLocate,
      B: GenomeLocate
{
    unimplemented!()
}

enum SpatialRelationship {
    FivePrimeNeighbor(usize),
    ThreePrimeNeighbor(usize),
    FivePrimeBoundaryOverlap,
    ThreePrimeBoundaryOverlap,
    Engulfing,
    DifferentReplicons,
}

trait GenomeLocate {
    fn locate(&self) -> GenomeLocation;
}

#[cfg(test)]
mod tests {

}

/*
    fn link_blast_data(&mut self, blast_map: BlastHitsTable) {

        // Pull out list of annotated genes from genome
        let genes = match &mut self.genes {
            None => panic!("ERROR: genome has no annotation data to cross-reference BLAST data against."),
            Some(x) => x,
        };

        // Pull out genome accession IDs for each replicon in genome
        let replicon_accessions = self.replicons.keys();

        //
        for replicon in replicon_accessions {

            // Pull list of BLAST hits for the given genome
            let replicon_specific_results = blast_map.results_table.get(replicon);

            match replicon_specific_results {
                None => continue, // If no results, test next replicon
                Some(hits) => {
                    for hit in hits {
                        for gene in genes.iter_mut() {
                            gene.blast_link(hit);
                        }
                    }
                }
            }


        }

    }
*/