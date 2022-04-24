#![allow(dead_code)]

use std::iter::Product;

use crate::import::{self, ProtoGenome, AssemblyMetadata, ProtoReplicon, ProtoRepliconType};

#[derive(Debug, PartialEq)]
enum Strand {
    Forward,
    ReverseComp,
}

#[derive(Debug, PartialEq, Clone)]
enum RepliconType {
    Chromosome,
    Plasmid,
    Other,
}

struct Genome {
    metadata: AssemblyMetadata,
    replicons: Vec<Replicon>,
    features: Option<Vec<Box<dyn GenomeLocation>>>
}

impl Genome {
    fn from_proto(input_data: ProtoGenome, metadata: AssemblyMetadata) -> Genome {

        let replicons = {
            let mut tmp: Vec<Replicon> = Vec::with_capacity(input_data.genomic_elements.len());

            for proto_replicon in input_data.genomic_elements {
                let new_replicon = Replicon::from_proto(proto_replicon);
                tmp.push(new_replicon);
            }
            
            tmp
        };

        Genome {
            metadata,
            replicons,
            features: None,
        };

        unimplemented!()
    }
}

#[derive(Debug, PartialEq)]
struct Replicon {
    accession_id:   String,
    sequence:       String,
    strand:         Strand,
    replicon_type:  RepliconType,
}

impl Replicon {
    fn from_proto(input_data: ProtoReplicon) -> Replicon {
        let accession_id = input_data.replicon_accession;
        let sequence = input_data.replicon_sequence;
        let strand = Strand::Forward;
        let replicon_type = match input_data.proto_replicon_type {
            ProtoRepliconType::Chromosome => RepliconType::Chromosome,
            ProtoRepliconType::Plasmid => RepliconType::Plasmid,
        };

        Replicon { 
            accession_id,
            sequence, 
            strand,
            replicon_type
        }
    }

    fn reverse_complement(&self) -> Replicon {
        let id = self.accession_id.clone();
        let replicon_type = self.replicon_type.clone();

        let strand = match self.strand {
            Strand::Forward =>      Strand::ReverseComp,
            Strand::ReverseComp =>  Strand::Forward,
        };

        let sequence = {
            let na_complement_remapping = |x| {
                match x {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    'U' => 'A',
                    'R' => 'Y',
                    'Y' => 'R',
                    'N' => 'N',
                     _  =>  x,
                    // Any other characters are turned back into themselves
                }
            };

            self.sequence.chars()
                         .rev()
                         .map(na_complement_remapping)
                         .collect::<String>()
        };

        Replicon {
            accession_id: id,
            sequence,
            strand,
            replicon_type,
        }
    }
}

trait GenomeLocation {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_1() {
        let test_replicon = Replicon {
            accession_id: "Test Replicon #1".to_string(),
            sequence: "AAATTTGGGCCC".to_string(),
            strand: Strand::Forward,
            replicon_type: RepliconType::Chromosome,
        };

        let correct_answer = Replicon {
            accession_id: "Test Replicon #1".to_string(),
            sequence: "GGGCCCAAATTT".to_string(),
            strand: Strand::ReverseComp,
            replicon_type: RepliconType::Chromosome,
        };

        let test_answer = test_replicon.reverse_complement();

        assert_eq!(correct_answer, test_answer);
    }

    #[test]
    fn test_reverse_complement_2() {
        let test_replicon = Replicon {
            accession_id: "Test Replicon #2".to_string(),
            sequence: "AAAAGGGGATATCGCG".to_string(),
            strand: Strand::ReverseComp,
            replicon_type: RepliconType::Plasmid,
        };

        let correct_answer = Replicon {
            accession_id: "Test Replicon #2".to_string(),
            sequence: "CGCGATATCCCCTTTT".to_string(),
            strand: Strand::Forward,
            replicon_type: RepliconType::Plasmid,
        };

        let test_answer = test_replicon.reverse_complement();

        assert_eq!(correct_answer, test_answer);
    }

    #[test]
    fn test_reverse_complement_3() {
        let test_replicon = Replicon {
            accession_id: "Test Replicon #3".to_string(),
            sequence: "AAAAGGGGATATCGCG".to_string(),
            strand: Strand::ReverseComp,
            replicon_type: RepliconType::Plasmid,
        };

        let intermediate = Replicon {
            accession_id: "Test Replicon #3".to_string(),
            sequence: "CGCGATATCCCCTTTT".to_string(),
            strand: Strand::Forward,
            replicon_type: RepliconType::Plasmid,
        };

        let test_answer = intermediate.reverse_complement();

        assert_eq!(test_replicon, test_answer);
    }

    #[test]
    #[should_panic]
    fn test_reverse_complement_4() {
        let test_replicon = Replicon {
            accession_id: "Test Replicon #4".to_string(),
            sequence: "AAAAGGGGATATCGCG".to_string(),
            strand: Strand::ReverseComp,
            replicon_type: RepliconType::Chromosome,
        };

        let intermediate = Replicon {
            accession_id: "Test Replicon #4".to_string(),
            sequence: "CGCGATATCCCCTTTT".to_string(),
            strand: Strand::Forward,
            replicon_type: RepliconType::Other,
        };

        let test_answer = intermediate.reverse_complement();

        assert_eq!(test_replicon, test_answer);
    }
}