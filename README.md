Check bimodality of a histogram using `van_der_eijk` function.

## Limiations

-   The first half of the histogram should not not be too heavy (>2x) or way too
    light (<0.5x) of the second half of the histogram i.e. if you sum first half
    of the histogram and second half of the histogram, either of them should not
    be more than twice the value than the other.
-   The peaks of the bimodal histogram must be at the start and end.

# Usage

This library provides two functions for both Rust and Python.

-   `van_der_eijk` returns A-score as described in [1]. If it is less than 0.0,
    the distribution is very likely bimodal.
-   `is_histogram_bimodal` is a wrapper on `van_der_eijk` function and it
    returns `True` if histogram is bimodal.

## Python

Here are some run on small histogram. This library typically performs much
better on larger histogram.

```bash
>>> is_bimodal.van_der_eijk([3, 1, 3])  # A-score (negative means bimodal)
-0.14285714285714285
>>> is_bimodal.is_histogram_bimodal([3, 1, 3])
True
>>>>> is_bimodal.is_histogram_bimodal([3, 2, 3])  # the peaks are not very high
False
>>> is_bimodal.is_histogram_bimodal([3, 2, 4])
False
>>> is_bimodal.is_histogram_bimodal([4, 2, 2, 4])
False
>>> is_bimodal.is_histogram_bimodal([4, 1, 2, 4])
True
```

See <python/tests/test_is_bimodal.py>. It has some test cases where the method works
poorly.


## Rust

Here are some examples from unit tests.

```
assert!(van_der_eijk(vec![30, 40, 210, 130, 530, 50, 10]) > 0.0);
assert!(van_der_eijk(vec![30, 40, 210, 10, 530, 50, 10]) > 0.0);
assert!(van_der_eijk(vec![30, 40, 10, 10, 30, 50, 100]) > 0.0);
assert!(van_der_eijk(vec![3, 4, 1, 1, 3, 5, 10]) > 0.0);
assert!(van_der_eijk(vec![3, 4, 1, 1, 3, 5, 1]) > 0.0);
assert!(van_der_eijk(vec![1, 1, 1, 1, 1, 1, 1]) > 0.0);
assert!(van_der_eijk(vec![1, 1, 1, 1, 1, 1, 1000]) > 0.0);

// bimodal and detected as bimodal.
assert!(van_der_eijk(vec![10000, 1, 1, 1, 1, 1, 10]) < 0.0);
assert!(van_der_eijk(vec![10, 10, 0, 0, 0, 10, 10]) < 0.0);
assert!(van_der_eijk(vec![10, 10, 0, 0, 0, 0, 10]) < 0.0);
assert!(van_der_eijk(vec![1, 1, 1, 0, 0, 1, 1]) < 0.0);
assert!(van_der_eijk(vec![1, 1, 1, 0, 1, 1, 1]) < 0.0);

// Test cases that bring the limitations of the algorithm.
// This should be bi-modal. Algo fails because weights are not balanced here.
// One side of the see-saw is 2x heavier.
assert!(van_der_eijk(vec![10, 11, 0, 0, 0, 0, 3, 3]) > 0.0);
assert!(van_der_eijk(vec![10, 11, 0, 0, 0, 0, 30, 31]) > 0.0);
assert!(van_der_eijk(vec![10, 11, 0, 0, 0, 0, 20, 11]) < 0.0);
```

# References

[1] Eijk, Cees. (2001). Measuring Agreement in Ordered Rating Scales. Quality
and Quantity. 35. 325-341. 10.1023/A:1010374114305.
