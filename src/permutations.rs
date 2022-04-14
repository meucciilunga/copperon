#![allow(dead_code)]
pub use std::collections::HashMap;

// Given a scaffold sequence and a remapping table defined for each target substitutable
// character in the scaffold, calculate all permutations of the string related to the FIRST 
// SUBSTITUTABLE POSITION found along the length of the string
fn next_site_permutator(sequence: &String, remapping_table: &HashMap<char, Vec<char>>) -> Option<Vec<String>> {

    // Storage logistics
    let mut permutation_list: Vec<String> = Vec::new();
    let mut tmp_seq = String::from(sequence).chars().collect::<Vec<char>>();

    // Substitution logic
    for (index, c) in sequence.chars().enumerate() {

        // See if remapping table has a value for the current character
        if !remapping_table.keys().collect::<Vec<&char>>().contains(&&c) {
            continue;
        }

        // Pull substitution table from hashmap for current char
        let err_msg = format!("ERROR: could not pull substitution table for character: '{}'!", &c);
        let substitution_table = remapping_table.get(&c).expect(&err_msg);
        
        // Perform replacement for each item in the substitution table
        for new_char in substitution_table  {
            tmp_seq[index] = *new_char;
            let new_permutation = tmp_seq.iter().collect::<String>();
            permutation_list.push(new_permutation);
        }

        // exit loop after first succesful substitution
        break;
    }

    if permutation_list.len() == 0 {
        None
    } else {
        Some(permutation_list)
    }
}


// Given a sequence scaffold and a remapping table defined for each target substitutable position,
// return a vector containing every possible permutation of the sequence
pub fn full_seq_permutator(sequence: String, remapping_table: HashMap<char, Vec<char>>) -> Vec<String> {
    
    // Storage logistics
    let mut curr_permutations: Vec<String> = vec![sequence];
    let mut resulting_permutations: Vec<String> = vec![];

    loop {
        for item in curr_permutations.iter().map(|x| next_site_permutator(x, &remapping_table)) {
            
            match item {
                Some(mut results_list) => {
                    resulting_permutations.append(&mut results_list)
                },
                None => break // when can confirm loop has reached the end of substitutable positions
            }
        }

        if resulting_permutations.len() == 0 {
            // If 'none' was reached in previous loop, means that curr_permutations list is complete
            // and thus can be returned to the caller
            break curr_permutations
        } else {
            curr_permutations = resulting_permutations;
            resulting_permutations = vec![];
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io::Read;
    use std::path::PathBuf;
    use std::fs::File;
    use std::collections::HashSet;
    use crate::cop_operon_specific;

    #[test]
    fn next_site_test_1() {
        let test_seq = "CXC".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
            let x_subs = vec!['1', '2', '3', '4'];
            tmp_table.insert('X', x_subs);
            tmp_table
        };

        let expected_result = {
            let permutations: [String; 4] = [
                "C1C".to_string(),
                "C2C".to_string(),
                "C3C".to_string(),
                "C4C".to_string(),
            ];

            Vec::from(permutations)
        };

        let actual_result = next_site_permutator(&test_seq, &permutation_table);

        assert_eq!(actual_result.unwrap(), expected_result);
    }

    #[test]
    fn next_site_test_2() {
        let test_seq = "CXCX".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
            let x_subs = vec!['1', '2', '3', '4'];
            tmp_table.insert('X', x_subs);
            tmp_table
        };

        let expected_result = {
            let permutations: [String; 4] = [
                "C1CX".to_string(),
                "C2CX".to_string(),
                "C3CX".to_string(),
                "C4CX".to_string(),
            ];

            Vec::from(permutations)
        };

        let actual_result = next_site_permutator(&test_seq, &permutation_table);

        assert_eq!(actual_result.unwrap(), expected_result);
    }

    #[test]
    fn next_site_test_3() {
        let test_seq = "1XCX".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
            let subs = vec!['L', 'M', 'N'];
            tmp_table.insert('1', subs);
            tmp_table
        };

        let expected_result = {
            let permutations: [String; 3] = [
                "LXCX".to_string(),
                "MXCX".to_string(),
                "NXCX".to_string(),
            ];

            Vec::from(permutations)
        };

        let actual_result = next_site_permutator(&test_seq, &permutation_table);

        assert_eq!(actual_result.unwrap(), expected_result);
    }

    #[test]
    #[should_panic]
    fn next_site_test_4() {
        let test_seq = "2X2Xasdas3YY".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
            let subs = vec!['L', 'M', 'N'];
            tmp_table.insert('1', subs);
            tmp_table
        };

        let actual_result = next_site_permutator(&test_seq, &permutation_table);

        // Since there is no avaialble position to be substituted in the test_seq, the next_site function
        // should return an Option::None--which will panic when unwrapped!
        actual_result.unwrap();
    }

    #[test]
    fn full_permutator_test_1() {
        let test_seq = "1X1XABCD".to_string();
        let permutation_table = {
            let mut tmp_table: HashMap<char, Vec<char>> = HashMap::new();
            let subs = vec!['A', 'B', 'C'];
            tmp_table.insert('1', subs);
            tmp_table
        };

        let expected_result = {
            let permutations: [String; 9] = [
                "AXAXABCD".to_string(),
                "AXBXABCD".to_string(),
                "AXCXABCD".to_string(),
                "BXAXABCD".to_string(),
                "BXBXABCD".to_string(),
                "BXCXABCD".to_string(),
                "CXAXABCD".to_string(),
                "CXBXABCD".to_string(),
                "CXCXABCD".to_string(),
            ];
            
            let set_of_permutations: HashSet<String> = permutations.into_iter().collect();
            set_of_permutations
        };

        let actual_result: HashSet<String> = full_seq_permutator(test_seq, permutation_table).into_iter().collect();

        assert_eq!(actual_result, expected_result);
    }

    #[test]
    fn full_permutator_test_2() {

        // Expected
        let expected_result = {
            let permutation_list_path = PathBuf::from("test_assets/RNYKACANNYGTMRNY_permutations.txt");
            let err_msg = "ERROR: could not open file containing pre-generated list of CopY permutations!";
            let mut file = File::open(permutation_list_path).expect(err_msg);
            let mut contents = String::new();

            file.read_to_string(&mut contents).expect("ERROR: could not pull data from file!");
            let permutations = contents.split('\n');
            
            let set_of_permutations: HashSet<String> = permutations.map(|x| x.to_string()).collect();
            set_of_permutations
        };

        // Actual
        let (test_seq, permutation_table) = cop_operon_specific::build_cop_permutation_table();
        let actual_result: HashSet<String> = full_seq_permutator(test_seq, permutation_table).into_iter().collect();

        assert_eq!(actual_result, expected_result);
    }

}