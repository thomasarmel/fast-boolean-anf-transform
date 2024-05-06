use num_traits::Unsigned;

/// Fast ANF transformation for cellular automata truth table rules expressed as unsigned integers
/// # Arguments
/// * `rule_number` - The rule number truth table to transform
/// * `num_variables_function` - The number of variables in the cellular automata function
/// # Returns
/// The ANF transformed rule number, as unsigned integer
/// # Panics
/// Panics if the rule number is greater than or equal to 2^(2^n), n being the number of variables in the function, as the rule number cannot exist
/// # Example
/// ```
/// use fast_boolean_anf_transform::fast_bool_anf_transform_unsigned;
/// assert_eq!(fast_bool_anf_transform_unsigned(3u32, 2), 5); // rule 3: 1 ^ x1
/// ```
pub fn fast_bool_anf_transform_unsigned<
    U: Unsigned
        + std::ops::Shr<u32>
        + std::ops::BitOr<u32, Output = U>
        + std::ops::BitAnd<u32, Output = U>
        + PartialOrd<u32>
        + Copy,
>(
    rule_number: U,
    num_variables_function: usize,
) -> U
where
    <U as std::ops::Shr<u32>>::Output: std::ops::BitAnd<u32>,
    <<U as std::ops::Shr<u32>>::Output as std::ops::BitAnd<u32>>::Output: PartialEq<u32>,
{
    #[cfg(debug_assertions)]
    if rule_number >= ((1 << (1 << num_variables_function)) as u32) {
        panic!("The rule number must be less than 2^(2^n), n being the number of variables in the function");
    }

    let mut blocksize = 1;
    let mut final_f = rule_number;
    for _ in 0..num_variables_function {
        let mut source: u32 = 0;
        while source < (1 << num_variables_function) {
            let target = source + blocksize;
            for i in 0..blocksize {
                let f_target_i: bool = ((final_f >> (target + i)) & 1) != 0;
                let f_source_i: bool = ((final_f >> (source + i)) & 1) != 0;
                let f_target_i_xor_f_source_i = f_target_i ^ f_source_i;
                if f_target_i_xor_f_source_i {
                    final_f = final_f | (1 << (target + i));
                } else {
                    final_f = final_f & !(1 << (target + i));
                }
            }
            source = source + (blocksize << 1);
        }
        blocksize = blocksize << 1;
    }
    final_f
}

/// Fast ANF transformation for cellular automata truth table rules expressed as boolean arrays
/// # Arguments
/// * `rule_truth_table` - The rule truth table to transform, first element is the output for input 0, second element is the output for input 1, etc.
/// # Returns
/// The ANF transformed rule truth table, as boolean array
/// # Panics
/// Panics if the rule truth table length is not equal to 2^n, n being the number of variables in the function
/// # Example
/// ```
/// use fast_boolean_anf_transform::fast_bool_anf_transform_bool_array;
/// let mut rule_truth_table = [false, true, false, false];
/// fast_bool_anf_transform_bool_array(&mut rule_truth_table);
/// assert_eq!(rule_truth_table, [false, true, false, true]); // rule 2: x0 ^ (x0 . x1)
pub fn fast_bool_anf_transform_bool_array(rule_truth_table: &mut [bool]) {
    let mut blocksize = 1;
    let input_size = rule_truth_table.len().trailing_zeros() as usize;

    #[cfg(debug_assertions)]
    if rule_truth_table.len() != 1 << input_size {
        panic!("The input truth table must have a size of 2^n, n being the number of variables in the function");
    }

    for _ in 0..input_size {
        let mut source = 0;
        while source < (1<<input_size) {
            let target = source + blocksize;
            for i in 0..blocksize {
                rule_truth_table[target + i] ^= rule_truth_table[source + i];
            }
            source = source + (blocksize << 1);
        }
        blocksize = blocksize << 1;
    }
}

#[cfg(test)]
mod tests {
    use super::{fast_bool_anf_transform_bool_array, fast_bool_anf_transform_unsigned};
    #[test]
    fn test_fast_bool_anf_transform_unsigned() {
        assert_eq!(fast_bool_anf_transform_unsigned(0u32, 2), 0);
        assert_eq!(fast_bool_anf_transform_unsigned(1u32, 2), 15);
        assert_eq!(fast_bool_anf_transform_unsigned(2u32, 2), 10);
        assert_eq!(fast_bool_anf_transform_unsigned(3u32, 2), 5);
        assert_eq!(fast_bool_anf_transform_unsigned(4u32, 2), 12);
        assert_eq!(fast_bool_anf_transform_unsigned(5u32, 2), 3);
        assert_eq!(fast_bool_anf_transform_unsigned(6u32, 2), 6);
        assert_eq!(fast_bool_anf_transform_unsigned(7u32, 2), 9);
        assert_eq!(fast_bool_anf_transform_unsigned(8u32, 2), 8);
        assert_eq!(fast_bool_anf_transform_unsigned(9u32, 2), 7);
        assert_eq!(fast_bool_anf_transform_unsigned(10u32, 2), 2);
        assert_eq!(fast_bool_anf_transform_unsigned(11u32, 2), 13);
        assert_eq!(fast_bool_anf_transform_unsigned(12u32, 2), 4);
        assert_eq!(fast_bool_anf_transform_unsigned(13u32, 2), 11);
        assert_eq!(fast_bool_anf_transform_unsigned(14u32, 2), 14);
        assert_eq!(fast_bool_anf_transform_unsigned(15u32, 2), 1);

        assert_eq!(fast_bool_anf_transform_unsigned(240u32, 3), 16);
        assert_eq!(fast_bool_anf_transform_unsigned(30u32, 3), 30);
    }

    #[test]
    #[should_panic]
    fn test_fast_bool_anf_transform_unsigned_rule_number_too_large() {
        let _ = fast_bool_anf_transform_unsigned(16 as u32, 2);
    }

    #[test]
    fn test_fast_bool_anf_transform_bool_array() {
        let mut rule_truth_table = [false, false, false, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, false, false, false]);

        rule_truth_table = [true, false, false, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, true, true, true]);

        rule_truth_table = [false, true, false, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, true, false, true]);

        rule_truth_table = [true, true, false, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, false, true, false]);

        let mut rule_truth_table = [false, false, true, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, false, true, true]);

        rule_truth_table = [true, false, true, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, true, false, false]);

        rule_truth_table = [false, true, true, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, true, true, false]);

        rule_truth_table = [true, true, true, false];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, false, false, true]);

        let mut rule_truth_table = [false, false, false, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, false, false, true]);

        rule_truth_table = [true, false, false, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, true, true, false]);

        rule_truth_table = [false, true, false, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, true, false, false]);

        rule_truth_table = [true, true, false, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, false, true, true]);

        let mut rule_truth_table = [false, false, true, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, false, true, false]);

        rule_truth_table = [true, false, true, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, true, false, true]);

        rule_truth_table = [false, true, true, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, true, true, true]);

        rule_truth_table = [true, true, true, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [true, false, false, false]);

        let mut rule_truth_table = [false, false, false, false, true, true, true, true]; // 240
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, false, false, false, true, false, false, false]);

        let mut rule_truth_table = [false, true, true, true, true, false, false, false]; // 30
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
        assert_eq!(rule_truth_table, [false, true, true, true, true, false, false, false]);
    }

    #[test]
    #[should_panic]
    fn test_fast_bool_anf_transform_bool_array_wrong_input_size() {
        let mut rule_truth_table = [false, false, false, false, true, true, true];
        fast_bool_anf_transform_bool_array(&mut rule_truth_table);
    }
}
