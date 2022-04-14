mod multithreader;
mod permutations;
mod genome;
mod import;

pub mod cop_operon_specific {
    use std::collections::HashMap;

    pub fn build_cop_permutation_table() -> (String, HashMap<char, Vec<char>>) {
        let operator_seq = "RNYKACANNYGTMRNY".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
    
            let character_flags = ['R', 'N', 'Y', 'K', 'M'];
            let substituents = [
                vec!['A', 'G'],
                vec!['A', 'C', 'G', 'T'],
                vec!['C', 'T'],
                vec!['G', 'T'],
                vec!['A', 'C'],
            ];
    
            for (key, val) in character_flags.into_iter().zip(substituents) {
                tmp_table.insert(key, val);
            }
            tmp_table
        };
    
        return(operator_seq, permutation_table)
    }
}