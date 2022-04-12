#![allow(non_snake_case)]
use std::collections::HashMap;
use std::io::Write;
use std::fs::File;

// Given a consensus sequence, performs a single substitution at
// the first substitutable position found along length of the string
fn single_site_permutator(consensus_seq: &str) -> Vec<String> {

    // sequence scaffold; only one site on this string will be mutated
    let mut permutable_seq = String::from(consensus_seq);

    // Arrays containing corresponding values to substitute
    let R_arr = ['A', 'G'];
    let N_arr = ['A', 'C', 'G', 'T'];
    let Y_arr = ['C', 'T'];
    let K_arr = ['G', 'T'];
    let M_arr = ['A', 'C'];
    let ignore_arr = ['I']; // 'I' for ignore; needed b/c match statements
                            //  require all possible outputs be idenitcal

    // Storage vector
    let mut permutation_list: Vec<String> = Vec::new();

    // Consider (and number) every position in consensus sequence string
    for (index, c) in consensus_seq.chars().enumerate() {

        // For a given position, determine list of required substitutions
        let substitution_arr = match c {
            'R' => &R_arr[..],
            'N' => &N_arr[..],
            'Y' => &Y_arr[..],
            'K' => &K_arr[..],
            'M' => &M_arr[..],
             _  => &ignore_arr[..],
        };
        
        // Move on to next character if the ignore flag is found
        if substitution_arr[0] == 'I' {
            continue;
        }

        // If valid substitution flag is found, perform substitutions
        for char in substitution_arr {
            let sub = char.to_string();
            permutable_seq.replace_range(index..index+1, &sub[..]);
            permutation_list.push(permutable_seq.clone());
        }
        
        break; // forces loop to stop after first substitutable position is located
    }

    permutation_list
}


// For a given consensus sequence, calculate all possible 
// permutations and store the resulting sequences in a vector. 
pub fn build_full_set(init_consensus_seq: &str) -> Vec<String> {
    
    // Calculate number of substitutable positions
    let sub_flags = "RNYKM";
    let mut sub_count = 0;

    for x in init_consensus_seq.chars() {
        if sub_flags.contains(x) {
            sub_count = sub_count + 1;
        }
    }

    // Perform substitutions
    let mut scope_transfer_list: Vec<String> = vec![]; // enables moving sequences generated during loop-scope into root scope of function
    let mut substituent_list = vec![String::from(init_consensus_seq)]; // this vector stores all strings undergoing the NEXT round of substitutions
    let mut tmp: Vec<String>;

    // Run the single-site substitution code the same number of times
    // as there are substitutable positions in the consensus sequence
    for _count in 0..sub_count {

        // Generate substitution permutations--considering only a single site 
        // at a time--and make a vector to store the resulting set of sequences.
        for seq in &substituent_list {
            tmp = single_site_permutator(&seq[..]);
            scope_transfer_list.append(&mut tmp);
        }

        // For next iteration, make newly generated consensus sequences
        // available for permutation processing; then, remove everything
        // from the loop-storage vector after the transfer
        substituent_list = scope_transfer_list.clone();         
        scope_transfer_list.retain(|_x| false);
    }

    // Last scope transfer from loop contains list of 
    // all possible permutations; return this list
    substituent_list
}


// Given a list of permutation sequences as a vector, export it to a text file.
pub fn _pregenerate_permutations(path: &str, permutation_list: Vec<String>) {
    let err_msg = "ERROR: Could not create file to output permuatation list data to.";
    let mut file = File::create(path).expect(err_msg);
    let mut line: String;

    for (i, seq) in permutation_list.iter().enumerate() {
        
        if i+1 != permutation_list.len() {
            line = format!("{}{}", seq, "\n");
        } else {
            line = seq.to_string();
        }
        
        file.write_all(line.as_bytes()).expect("ERROR: Failed to write line to file.");
    }
}

// Create a symmetric hash-table from a list of permutation_sequences
pub fn HashZip(seq_permutations: &Vec<String>) -> HashMap<&str, &str> {

    // Create HashMap where we can see if a seq is in the HashTable
    // by directly querying that seq
    let mut HashZip: HashMap<&str, &str> = HashMap::new();
    for seq in seq_permutations {
        HashZip.insert(&seq[..], &seq[..]);
    }

    HashZip
}


// Unit Tests
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn simple_R_substitution() {
        let test_val = single_site_permutator("RRR");
        let actual_val = vec!["ARR".to_string(), "GRR".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn simple_N_substitution() {
        let test_val = single_site_permutator("NNN");
        let actual_val = vec!["ANN".to_string(), "CNN".to_string(), "GNN".to_string(), "TNN".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn simple_Y_substitution() {
        let test_val = single_site_permutator("YYY");
        let actual_val = vec!["CYY".to_string(), "TYY".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn simple_K_substitution() {
        let test_val = single_site_permutator("KKK");
        let actual_val = vec!["GKK".to_string(), "TKK".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn simple_M_substitution() {
        let test_val = single_site_permutator("MMM");
        let actual_val = vec!["AMM".to_string(), "CMM".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn arbitrary_pos_R_substitution() {
        let test_val = single_site_permutator("AATLATLLARSMMMM");
        let actual_val = vec![String::from("AATLATLLAASMMMM"), String::from("AATLATLLAGSMMMM")];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn arbitrary_pos_N_substitution() {
        let test_val = single_site_permutator("OINUEWRQWE");
        let actual_val = vec!["OIAUEWRQWE".to_string(), "OICUEWRQWE".to_string(), "OIGUEWRQWE".to_string(), "OITUEWRQWE".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn arbitrary_pos_Y_substitution() {
        let test_val = single_site_permutator("AGTCSYYYSDFSDFSDFSDFSDFS");
        let actual_val = vec!["AGTCSCYYSDFSDFSDFSDFSDFS".to_string(), "AGTCSTYYSDFSDFSDFSDFSDFS".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn arbitrary_pos_K_substitution() {
        let test_val = single_site_permutator("ZXCZXCZXASADAKZCZCXZCZCZCASDA");
        let actual_val = vec!["ZXCZXCZXASADAGZCZCXZCZCZCASDA".to_string(), "ZXCZXCZXASADATZCZCXZCZCZCASDA".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn arbitrary_pos_M_substitution() {
        let test_val = single_site_permutator("QWEEMQ");
        let actual_val = vec!["QWEEAQ".to_string(), "QWEECQ".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn full_substitution_1() {
        let test_val = build_full_set("KR");
        let actual_val = vec!["GA".to_string(), "GG".to_string(), "TA".to_string(), "TG".to_string()];
        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn full_substitution_2() {
        let output_vals = ["AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                           "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT"];
        let test_val = build_full_set("MKN");
        let mut actual_val: Vec<String> = vec![];

        for seq in output_vals {
            actual_val.push(seq.to_string());
        }

        assert_eq!(test_val, actual_val);
    }

    #[test]
    fn expected_num_of_CopY_Operator_permutations() {

        let substituent_list = build_full_set("RNYKACANNYGTMRNY");
        const EXPECTED_PERM_COUNT: usize = 32_768;

        if substituent_list.len() != EXPECTED_PERM_COUNT {
            // Length of raw output list
            println!("Length of Calculated Permutations Length: {}", substituent_list.len());
    
            // Length of hash-set
            println!("Expected Number of Permutations for 'RNYKACANNYGTMRNY': {}", EXPECTED_PERM_COUNT);
    
            let msg = "ERROR: permutation generator does not produce expected number of permutations for CopY consensus sequence!";
            panic!("{}", msg);
        }
    }

    #[test]
    fn list_of_CopY_Operator_permutations_is_unique() {

        let substituent_list = build_full_set("RNYKACANNYGTMRNY");
        let mut permutation_verify = HashSet::new();
        for seq in &substituent_list {
            permutation_verify.insert(seq);
        }

        if substituent_list.len() != permutation_verify.len() {
            
            println!("\nResults of Uniqueness Test:");

            // Length of raw output list
            println!("Substitution List Length: {}", substituent_list.len());
    
            // Length of hash-set
            println!("Substitution Hash-Set Length: {}", permutation_verify.len());
    
            let msg = "ERROR: list of permutations contains non-unique element(s)!";
            panic!("{}", msg);
        }
    }

}