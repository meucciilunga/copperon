#![allow(dead_code)]
#![allow(unused_imports)]

use crate::import::{self, ProtoGenome, AssemblyMetadata, ProtoReplicon, ProtoRepliconType, AnnotationEntry, BlastDerivedAnnotation};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq)]
pub enum StrandSense {
    Forward,
    Reverse,
    Other,
}

#[derive(Debug, PartialEq, Clone)] 
pub enum RepliconType {
    Chromosome,
    Plasmid,
    Other,
}

#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum FeatureType {
    Region,
    Gene,
    CDS,
    Pseudogene,
    tRNA,
    BlastResult,
    Other,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GenomeLocation {
    pub replicon_accession: String,
    pub replicon_strand:    StrandSense,
    pub start_index_ord:    usize,
    pub end_index_ord:     usize,
}

pub struct Operator {
    pub location: GenomeLocation,
    pub seq: String,
    pub potential_operon: Vec<Gene>
}

// Wrapper for BLAST annotation data
#[derive(Clone, Debug, PartialEq)]
pub struct BlastFragment {
    pub location: GenomeLocation,
    pub length: usize,
    pub pident: f64,
    pub evalue: f64,
    pub qseqid: String,
    pub qstart: usize,
    pub qend:   usize,
}

impl BlastFragment {
    fn new(input: &BlastDerivedAnnotation) -> BlastFragment {
        let location = GenomeLocation {
            replicon_accession: input.sseqid.clone(),
            replicon_strand: if input.sframe.is_positive() {
                                 StrandSense::Forward
                             } else if input.sframe.is_negative() {
                                 StrandSense::Reverse
                             } else {
                                 StrandSense::Other
                             },
            start_index_ord: input.sstart,
            end_index_ord: input.send,
        };

        let length = input.match_length;
        let pident = input.pident;
        let evalue = input.evalue;
        let qseqid = input.qseqid.clone();
        let qstart = input.qstart;
        let qend = input.qend;
        
        BlastFragment {
            location,
            length,
            pident,
            evalue,
            qseqid,
            qstart,
            qend,
        }
    }
}

// Wrapper for accessing BLAST annotation data via genome accession id
pub struct BlastHitsTable {
    pub table_name: String,
    pub results_table: HashMap<String, Vec<BlastFragment>>,
    // The HashMap links every genome (via its accession id) to all known
    // BLAST hits found within that genome for a given BLAST search; the
    // table is derived from a vector containing all known BLAST results
    // from a given search
}

impl BlastHitsTable {

    // Turns a vector of parsed BLAST Fragments into a table containing
    // the BLAST Fragments linked to (and accessible by) their corresponding 
    // genome accession id
    pub fn new(table_name: String, input_blast_results: Vec<BlastFragment>) -> BlastHitsTable {
        
        let results_table = {

            // List of unique genome accessions
            let mut list_of_genome_accessions: HashSet<String> = HashSet::with_capacity(input_blast_results.len());
            for result in &input_blast_results {
                list_of_genome_accessions.insert(result.location.replicon_accession.clone());
            }

            // Primary storage structure
            let mut genome_to_hits: HashMap<String, Vec<BlastFragment>> = HashMap::with_capacity(list_of_genome_accessions.len());

            // Map a list of BlastFragements derived from the
            // same Genome to that Genome's accession id 
            for genome_id in list_of_genome_accessions {
                
                // For storing BlastFragments corresponding to the same genome
                let mut tmp_blast_results: Vec<BlastFragment> = Vec::new();

                for result in &input_blast_results {
                    if result.location.replicon_accession == genome_id {
                        tmp_blast_results.push(result.clone())
                    }
                }

                genome_to_hits.insert(genome_id, tmp_blast_results);
            }

            genome_to_hits
        };

        BlastHitsTable {
            table_name,
            results_table,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Gene {
    pub location: GenomeLocation,
    pub feature_type: FeatureType,
    pub attributes: HashMap<String, String>,
    pub blast_association: Option<BlastFragment>,
}

impl Gene {
    pub fn from_annotation_entry(input: &AnnotationEntry) -> Gene {

        let location = GenomeLocation {
            replicon_accession: input.genomic_accession.clone(),
            replicon_strand: match input.replicon_strand {
                '+' => StrandSense::Forward,
                '-' => StrandSense::Reverse,
                 _  => StrandSense::Other,
            },
            start_index_ord: input.start_index_ord,
            end_index_ord: input.end_index_ord,
        };

        let feature_type = match input.feature_type.to_lowercase().as_str() {
            "region" =>     FeatureType::Region,
            "gene" =>       FeatureType::Gene,
            "cds" =>        FeatureType::CDS,
            "pseudogene" => FeatureType::Pseudogene,
            "trna" =>       FeatureType::tRNA,
            _ =>            FeatureType::Other,
        };

        Gene {
            location,
            feature_type,
            attributes: input.attributes.clone(),
            blast_association: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Replicon {
    pub accession_id: String,
    pub replicon_type: RepliconType,
    pub fwd_strand: String,
    pub rev_strand: String,
}

// Special method specific to this module that allows for the generation
// of the reverse complement of a given ProtoReplicon's sequence
impl ProtoReplicon {

    fn reverse_complement(&self) -> String {
        let seq = {
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
                     _  =>  panic!("ERROR: replicon sequence contains invalid character."),
                    // Any other characters are turned back into themselves
                }
            };

            self.replicon_sequence.chars()
                                  .rev()
                                  .map(na_complement_remapping)
                                  .collect::<String>()
        };

        seq
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Genome {
    pub metadata:   AssemblyMetadata,
    pub replicons:  HashMap<String, Replicon>,
    pub genes:      Option<Vec<Gene>>, 
}

impl Genome {
    pub fn from_proto(input_data: ProtoGenome, 
                  metadata: AssemblyMetadata,
                  genes: Option<Vec<Gene>>) -> Genome {
        
        // Storage Logistics
        let mut replicons: HashMap<String, Replicon> = HashMap::with_capacity(2 * input_data.proto_replicons.len());

        // Convert ProtoReplicons into Replicons by combining initial replicon strands with their RevComps into a struct
        for proto_replicon_strand in input_data.proto_replicons {

            // Build replicon elements
            let replicon_type = match proto_replicon_strand.proto_replicon_type {
                ProtoRepliconType::Chromosome => RepliconType::Chromosome,
                ProtoRepliconType::Plasmid => RepliconType::Plasmid,
            };
            let fwd_strand = proto_replicon_strand.replicon_sequence.clone();
            let rev_strand = proto_replicon_strand.reverse_complement();

            let new_replicon = Replicon {
                accession_id: proto_replicon_strand.replicon_accession.clone(),
                replicon_type,
                fwd_strand,
                rev_strand,
            };

            // Map new replicon to its accession id
            replicons.insert(new_replicon.accession_id.clone(), new_replicon);
        }

        // Link genome metadata, genome sequence data, and genome annotation data into unified data structure
        Genome {
            metadata,
            replicons,
            genes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::import;
    use std::path::PathBuf;
    
    #[test]
    fn build_blast_fragments_1() {
        let test_blast_results_file = PathBuf::from("tests/test_assets/CopA_blast_result.txt");
        let test_annotation = import::parse_annotations_from_blast_results(test_blast_results_file);

        let test_fragment_1 = BlastFragment::new(&test_annotation[0]);
        let test_fragment_2 = BlastFragment::new(&test_annotation[314]);
        let test_fragment_3 = BlastFragment::new(&test_annotation[628]);
        let test_fragments = [test_fragment_1, test_fragment_2, test_fragment_3];

        let actual_fragment_1 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NC_003028.3".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 692412,
                end_index_ord: 694652,
            },
            length: 747,
            pident: 100.000,
            evalue: 0.0,
            qseqid: "WP_000136310.1".to_string(),
            qstart: 1,
            qend: 747,
        };
        let actual_fragment_2 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_CP038808.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 1446972,
                end_index_ord: 1445413,
            },
            length: 553,
            pident: 25.136,
            evalue: 1.30e-19,
            qseqid: "WP_000136310.1".to_string(),
            qstart: 160,
            qend: 655,
        };
        let actual_fragment_3 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_AP018338.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 1762342,
                end_index_ord: 1760876,
            },
            length: 523,
            pident: 27.533,
            evalue: 4.20e-53,
            qseqid: "WP_000136310.1".to_string(),
            qstart: 157,
            qend: 664,
        };
        let actual_fragments = [actual_fragment_1, actual_fragment_2, actual_fragment_3];

        for (actual, test) in actual_fragments.iter().zip(test_fragments.iter()) {
            assert_eq!(actual, test);
        }
    }

    #[test]
    fn build_blast_fragments_2() {
        let test_blast_results_file = PathBuf::from("tests/test_assets/CopY_blast_result.txt");
        let test_annotation = import::parse_annotations_from_blast_results(test_blast_results_file);

        let test_fragment_1 = BlastFragment::new(&test_annotation[48175]);
        let test_fragment_2 = BlastFragment::new(&test_annotation[77060]);
        let test_fragment_3 = BlastFragment::new(&test_annotation[117439]);
        let test_fragments = [test_fragment_1, test_fragment_2, test_fragment_3];

        let actual_fragment_1 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NC_021181.2".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 1973909,
                end_index_ord: 1974298,
            },
            length: 130,
            pident: 40.000,
            evalue: 1.40e-22,
            qseqid: "WP_003097211.1".to_string(),
            qstart: 3,
            qend: 132,
        };
        let actual_fragment_2 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "CP031018.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 2155972,
                end_index_ord: 2156385,
            },
            length: 138,
            pident: 40.580,
            evalue: 2.59e-26,
            qseqid: "WP_108033676.1".to_string(),
            qstart: 3,
            qend: 140,
        };
        let actual_fragment_3 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NC_020450.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 818375,
                end_index_ord: 818758,
            },
            length: 130,
            pident: 43.846,
            evalue: 6.13e-24,
            qseqid: "WP_042492010.1".to_string(),
            qstart: 5,
            qend: 134,
        };
        let actual_fragments = [actual_fragment_1, actual_fragment_2, actual_fragment_3];

        for (actual, test) in actual_fragments.iter().zip(test_fragments.iter()) {
            assert_eq!(actual, test);
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_blast_hits_table_CopA() {
        
        // Parse BLAST results file
        let test_blast_results_file = PathBuf::from("tests/test_assets/CopA_blast_result.txt");
        let test_annotations = import::parse_annotations_from_blast_results(test_blast_results_file).iter()
                                                                                                    .map(|x| BlastFragment::new(x))
                                                                                                    .collect::<Vec<BlastFragment>>();
        let test_blast_table = BlastHitsTable::new("CopA".to_string(), test_annotations);

        // Test that we have the expected number of keys for CopA; value for 
        // correct answer was derived from Python Script whose source code
        // is located here: tests/test_assets/CopA_blast_result.txt
        assert_eq!(25000, test_blast_table.results_table.len());

        // Test that we have the expect number of BlastFragment entries after building the hashtable
        let mut tot_num_entries = 0;
        for blast_results_list in test_blast_table.results_table.values() {
            tot_num_entries += blast_results_list.len();
        }
        assert_eq!(101_314, tot_num_entries);
        // 101,314 == total number of lines in CopA_blast_result.txt file;
        // each line should correspond to a unique BLAST result

        // Test-case 1
        let mut test_binding_1: HashMap<String, Vec<BlastFragment>> = HashMap::with_capacity(1);
        let binding_name = "NZ_AKVY01000001.1".to_string();
        let actual_fragment_1 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_AKVY01000001.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 694341,
                end_index_ord: 696581,
            },
            length: 747,
            pident: 100.000,
            evalue: 0.0,
            qseqid: "WP_000136310.1".to_string(),
            qstart: 1,
            qend:   747,
        };
        let actual_fragment_2 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_AKVY01000001.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 2014268,
                end_index_ord: 2012802,
            },
            length: 523,
            pident: 27.916,
            evalue: 3.53e-53,
            qseqid: "WP_000136310.1".to_string(),
            qstart: 157,
            qend:   664,
        };
        let actual_fragment_3 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_AKVY01000001.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 1527501,
                end_index_ord: 1525942,
            },
            length: 553,
            pident: 25.136,
            evalue: 1.89e-19,
            qseqid: "WP_000136310.1".to_string(),
            qstart: 160,
            qend:   655,
        };
        let binding_vec = vec![actual_fragment_1, actual_fragment_2, actual_fragment_3];
        test_binding_1.insert(binding_name, binding_vec);

        let actual = test_blast_table.results_table.get("NZ_AKVY01000001.1").unwrap();
        let test = test_binding_1.get("NZ_AKVY01000001.1").unwrap();
        assert_eq!(actual, test);
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_blast_hits_table_CopY() {

        // Parse BLAST results file
        let test_blast_results_file = PathBuf::from("tests/test_assets/CopY_blast_result.txt");
        let test_annotations = import::parse_annotations_from_blast_results(test_blast_results_file).iter()
                                                                                                    .map(|x| BlastFragment::new(x))
                                                                                                    .collect::<Vec<BlastFragment>>();
        let test_blast_table = BlastHitsTable::new("CopY".to_string(), test_annotations);
        
        // Test that we have the expected number of keys for CopY; value for 
        // correct answer was derived from Python Script whose source code
        // is located here: tests/test_assets/CopA_blast_result.txt
        assert_eq!(4966, test_blast_table.results_table.len());

        // Test that we have the expect number of BlastFragment entries after building the hashtable
        let mut tot_num_entries = 0;
        for blast_results_list in test_blast_table.results_table.values() {
            tot_num_entries += blast_results_list.len();
        }
        assert_eq!(123_173, tot_num_entries);
        // 123,173 == total number of lines in CopY_blast_result.txt file;
        // each line should correspond to a unique BLAST result


        // Test-case 1
        let actual_fragment = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_CP061021.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 1527375,
                end_index_ord: 1526989,
            },
            length: 131,
            pident: 38.931,
            evalue: 6.49e-26,
            qseqid: "WP_001167211.1".to_string(),
            qstart: 1,
            qend:   131,
        };

        // Test to see if first element of results table matches
        let test_fragment = test_blast_table.results_table.get("NZ_CP061021.1").unwrap()[0].clone();
        assert_eq!(actual_fragment, test_fragment);



        // Test-case 2
        let actual_fragment = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_AP024837.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 1823308,
                end_index_ord: 1822901,
            },
            length: 138,
            pident: 34.783,
            evalue: 4.95e-25,
            qseqid: "WP_018030946.1".to_string(),
            qstart: 3,
            qend:   140,
        };

        // Test to see if first element of results table matches
        let test_fragment = test_blast_table.results_table.get("NZ_AP024837.1").unwrap()[15].clone();
        assert_eq!(actual_fragment, test_fragment);



        // Test-case 3 -- interesting opportunity to test two fragments at once since they are consecutive
        // in the actual BLAST results table
        let actual_fragment_1 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_CP014144.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 107183,
                end_index_ord: 107587,
            },
            length: 135,
            pident: 57.037,
            evalue: 1.46e-43,
            qseqid: "WP_108033676.1".to_string(),
            qstart: 9,
            qend:   143,
        };

        let actual_fragment_2 = BlastFragment {
            location: GenomeLocation {
                replicon_accession: "NZ_CP014144.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 1641099,
                end_index_ord: 1641497,
            },
            length: 135,
            pident: 46.667,
            evalue: 3.13e-32,
            qseqid: "WP_108033676.1".to_string(),
            qstart: 9,
            qend:   143,
        };
        let actual_fragments = vec![actual_fragment_1, actual_fragment_2];
        //let actual = &actual_fragments[..];

        // Test to see if first element of results table matches
        let test = test_blast_table.results_table.get("NZ_CP014144.1").unwrap();

        assert_eq!(&actual_fragments[..], &test[40..=41]);
    }
    
    #[test]
    fn gene_from_annotation_entry() {

        // Parse genome annotation file
        let test_blast_results_file = PathBuf::from("tests/test_assets/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.gff");
        let test_annotations = import::parse_genome_annotation(test_blast_results_file);

        let test_gene_1 = Gene::from_annotation_entry(&test_annotations[0]);
        let test_gene_2 = Gene::from_annotation_entry(&test_annotations[6738]);
        let test_gene_3 = Gene::from_annotation_entry(&test_annotations[331]);
        let test_gene_4 = Gene::from_annotation_entry(&test_annotations[6998]);

        // Test-case 1
        let actual_gene_1 = Gene {
            location: GenomeLocation {
                replicon_accession: "NC_000911.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 1,
                end_index_ord: 3573470,
            },
            feature_type: FeatureType::Region,
            attributes: {
                "ID=NC_000911.1:1..3573470;\
                Dbxref=taxon:1148;\
                Is_circular=true;\
                Name=ANONYMOUS;\
                gbkey=Src;\
                genome=chromosome;\
                mol_type=genomic DNA;\
                note=synonym:Synechocystis PCC6803;\
                old-name=Synechocystis sp. PCC 6803;\
                strain=PCC 6803".split(';')
                                .map(|x| x.split('=')
                                          .map(|s| s.to_string())
                                          .collect::<Vec<String>>())
                                .collect::<Vec<Vec<String>>>()
                                .into_iter()
                                .map(|x| (x[0].clone(), x[1].clone()))
                                .collect::<HashMap<String, String>>()
            },
            blast_association: None,
        };
        assert_eq!(actual_gene_1, test_gene_1);
        

        // Test-case 2
        let actual_gene_2 = Gene {
            location: GenomeLocation {
                replicon_accession: "NC_005230.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 1,
                end_index_ord: 103307,
            },
            feature_type: FeatureType::Region,
            attributes: {
                "ID=NC_005230.1:1..103307;\
                Dbxref=taxon:1148;Is_circular=true;\
                Name=pSYSA;gbkey=Src;\
                genome=plasmid;\
                mol_type=genomic DNA;\
                note=synonym:Synechocystis PCC6803;\
                plasmid-name=pSYSA;strain=PCC 6803".split(';')
                                                   .map(|x| x.split('=')
                                                             .map(|s| s.to_string())
                                                             .collect::<Vec<String>>())
                                                    .collect::<Vec<Vec<String>>>()
                                                    .into_iter()
                                                    .map(|x| (x[0].clone(), x[1].clone()))
                                                    .collect::<HashMap<String, String>>()
            },
            blast_association: None,
        };
        assert_eq!(actual_gene_2, test_gene_2);


        // Test-case 3
        let actual_gene_3 = Gene {
            location: GenomeLocation {
                replicon_accession: "NC_000911.1".to_string(),
                replicon_strand: StrandSense::Reverse,
                start_index_ord: 175143,
                end_index_ord: 175682,
            },
            feature_type: FeatureType::CDS,
            attributes: {
                "ID=cds-WP_010871367.1;\
                Parent=gene-SGL_RS02685;\
                Dbxref=Genbank:WP_010871367.1;\
                Name=WP_010871367.1;\
                gbkey=CDS;\
                inference=COORDINATES: similar to AA sequence:RefSeq:WP_010871367.1;\
                locus_tag=SGL_RS02685;\
                product=F0F1 ATP synthase subunit B;protein_id=WP_010871367.1;\
                transl_table=11".split(';')
                                .map(|x| x.split('=')
                                            .map(|s| s.to_string())
                                            .collect::<Vec<String>>())
                                .collect::<Vec<Vec<String>>>()
                                .into_iter()
                                .map(|x| (x[0].clone(), x[1].clone()))
                                .collect::<HashMap<String, String>>()
            },
            blast_association: None,
        };
        assert_eq!(actual_gene_3, test_gene_3);


        // Test-case 4
        let actual_gene_4 = Gene {
            location: GenomeLocation {
                replicon_accession: "NC_005231.1".to_string(),
                replicon_strand: StrandSense::Forward,
                start_index_ord: 25147,
                end_index_ord: 27345,
            },
            feature_type: FeatureType::CDS,
            attributes: {
                "ID=cds-WP_011153801.1;\
                Parent=gene-SGL_RS01275;\
                Dbxref=Genbank:WP_011153801.1;\
                Name=WP_011153801.1;\
                gbkey=CDS;\
                inference=COORDINATES: similar to AA sequence:RefSeq:WP_011153801.1;\
                locus_tag=SGL_RS01275;\
                product=DUF839 domain-containing protein;\
                protein_id=WP_011153801.1;\
                transl_table=11".split(';')
                                .map(|x| x.split('=')
                                            .map(|s| s.to_string())
                                            .collect::<Vec<String>>())
                                .collect::<Vec<Vec<String>>>()
                                .into_iter()
                                .map(|x| (x[0].clone(), x[1].clone()))
                                .collect::<HashMap<String, String>>()
            },
            blast_association: None,
        };
        assert_eq!(actual_gene_4, test_gene_4);
    }

    #[test]
    fn test_proto_replicon_rev_comp_1() {
        let test_proto = import::ProtoReplicon {
            replicon_accession: "TEST_1".to_string(),
            replicon_sequence:  "ATGCACNNNNNATCG".to_string(),
            proto_replicon_type: ProtoRepliconType::Plasmid
        };

        let actual_rc = "CGATNNNNNGTGCAT".to_string();
        let test_rc = test_proto.reverse_complement();

        assert_eq!(actual_rc, test_rc);
    }

    #[test]
    fn test_proto_replicon_rev_comp_2() {
        let test_proto = import::ProtoReplicon {
            replicon_accession: "TEST_2".to_string(),
            replicon_sequence:  "NNNAAAYYYNNNGGGCCCAAATTTAAARRRUUUU".to_string(),
            proto_replicon_type: ProtoRepliconType::Plasmid
        };

        let actual_rc = "AAAAYYYTTTAAATTTGGGCCCNNNRRRTTTNNN".to_string();
        let test_rc = test_proto.reverse_complement();

        assert_eq!(actual_rc, test_rc);
    }

    #[test]
    #[should_panic]
    fn test_proto_replicon_rev_comp_3() {
        let test_proto = import::ProtoReplicon {
            replicon_accession: "TEST_3".to_string(),
            replicon_sequence:  "aacAdadcc129".to_string(),
            proto_replicon_type: ProtoRepliconType::Plasmid
        };

        test_proto.reverse_complement();
    }

    #[test]
    fn test_genome_from_proto_genome() {

        // Defining input
        let test_proto_1 = import::ProtoReplicon {
            replicon_accession: "TEST_1".to_string(),
            replicon_sequence:  "ATGCACNNNNNATCG".to_string(),
            proto_replicon_type: ProtoRepliconType::Chromosome
        };

        let test_proto_2 = import::ProtoReplicon {
            replicon_accession: "TEST_2".to_string(),
            replicon_sequence:  "NNNAAAYYYNNNGGGCCCAAATTTAAARRRUUUU".to_string(),
            proto_replicon_type: ProtoRepliconType::Plasmid
        };
        
        let protogenome = ProtoGenome {
            assembly_name: "TEST_PROTOGENOME_1".to_string(),
            proto_replicons: vec![test_proto_1, test_proto_2]
        };

        let metadata = AssemblyMetadata {
            assembly_accession: "TEST_ASSEMBLY".to_string(),
            biosample_accession: "TEST_BIOSAMPLE".to_string(),
            ref_seq_cat: "TEST_REF_SEQ_CAT".to_string(),
            tax_id: "1000".to_string(),
            species_tax_id: "2000".to_string(),
            organism_name: "TEST-ORGANISM-1".to_string(),
            infraspecifc_name: "TEST-INFRASPECIFIC-NAME".to_string(),
            assembly_level: "COMPLETE".to_string(),
            assembly_name: "TEST-ASSEMBLY".to_string(),
            ftp_location: "TEST-FTP-LOC".to_string(),
            subdir_name: "TEST-SUBDIR-NAME".to_string(),
            ref_seq_status: "TEST-STATUS".to_string(),
        };

        let actual_replicon_1 = Replicon {
            accession_id: "TEST_1".to_string(),
            replicon_type: RepliconType::Chromosome,
            fwd_strand: "ATGCACNNNNNATCG".to_string(),
            rev_strand: "CGATNNNNNGTGCAT".to_string(),
        };

        let actual_replicon_2 = Replicon {
            accession_id: "TEST_2".to_string(),
            replicon_type: RepliconType::Plasmid,
            fwd_strand: "NNNAAAYYYNNNGGGCCCAAATTTAAARRRUUUU".to_string(),
            rev_strand: "AAAAYYYTTTAAATTTGGGCCCNNNRRRTTTNNN".to_string(),
        };

        let mut replicon_map: HashMap<String, Replicon> = HashMap::with_capacity(2);
        replicon_map.insert("TEST_1".to_string(), actual_replicon_1);
        replicon_map.insert("TEST_2".to_string(), actual_replicon_2);

        let actual_genome = Genome {
            metadata: metadata.clone(),
            replicons: replicon_map,
            genes: None,
        };

        let test_genome = Genome::from_proto(protogenome, metadata, None);

        assert_eq!(actual_genome, test_genome);
    }
}