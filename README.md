# Fast Algebraic Normal Form Transform for boolean functions

*We used the function presented in this paper: [[1]](#1), of which we would like to thank the authors.*

---

## ANF transform

ANF transform is a method to convert a boolean function **from its truth table** representation **to its Algebraic Normal Form** (ANF) representation.
ANF is a representation of a boolean function as **a XOR sum of AND monomials**.

For example, let's consider the following 2-variables boolean function truth table:

| x1 | x0 | f(x1, x2) |
|----|----|-----------|
| 0  | 0  | 1         |
| 0  | 1  | 1         |
| 1  | 0  | 0         |
| 1  | 1  | 0         |

The ANF representation of this function is: `f(x1, x0) = 1 XOR x1`.

### Boolean function degree

The degree of a boolean function is the maximum degree of its monomials in its ANF representation.

For example, the degree of the function `f(x1, x0) = x0 XOR x0.x1` is 2.

### Unsigned integer representation

Truth tables are often represented as unsigned integers.
For example, the truth table of the following function:

| x1 | x0 | f(x1, x2) |
|----|----|-----------|
| 0  | 0  | 1         |
| 0  | 1  | 1         |
| 1  | 0  | 0         |
| 1  | 1  | 0         |

has the unsigned integer representation: `0b0011 = 3`.

In the same way, the ANF representation of boolean functions can be represented as unsigned integers.
For example, the ANF representation of the previous function `f(x1, x0) = 1 XOR x1`.
This can be written as `0b0101 = 5`.
To find the function from the ANF representation in binary, separate each bit and indicate its position starting **from the least significant bit**.

Position **0 corresponds to constant 1**, for the other positions convert them to binary.
This binary number gives the monomial, starting with x0 for the least significant bit.
All you have to do is multiply this monomial by the corresponding bit of the ANF representation.

|                               |               |            |            |           |
|-------------------------------|---------------|------------|------------|-----------|
| **ANF binary representation** | 0             | 1          | 0          | 1         |
| **bit position**              | 3             | 2          | 1          | 0         |
| **bit position in binary**    | 0b11          | 0b10       | 0b01       | /         |
| **monomial**                  | x1.x0         | x1         | x0         | 1         |
| **function**                  | **0**.(x1.x0) | **1**.(x1) | **0**.(x0) | **1**.(1) |

The final function is: `f(x1, x0) = x1 XOR 1`. Its degree is 1.


## Library usage

This library provides 2 functions to convert a boolean function from its truth table representation to its ANF representation.

### Boolean array representation

The first function `fast_bool_anf_transform_bool_array` takes a boolean array as input.

```rust
use fast_boolean_anf_transform::fast_bool_anf_transform_bool_array;

// f(0, 0) = true, f(0, 1) = true, f(1, 0) = false, f(1, 1) = false
let mut rule_truth_table = [true, true, false, false];
fast_bool_anf_transform_bool_array(&mut rule_truth_table);
// f(x1, x0) = 1 XOR x1
assert_eq!(rule_truth_table, [true, false, true, false]);
```

### Unsigned integer representation

The second function `fast_bool_anf_transform_unsigned` takes an unsigned integer as input.
You have to specify the number of variables of the boolean function.

#### Warning: if the size of the specified unsigned type is too small, the function will panic (debug only).

```rust
use fast_boolean_anf_transform::fast_bool_anf_transform_unsigned;
assert_eq!(fast_bool_anf_transform_unsigned(3u32, 2), 5);
```

## References
<a id="1">[1]</a>
Bakoev, Valentin. Fast implementation of the ANF transform.
