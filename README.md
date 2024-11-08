Check bimodality of a distribution using `van_der_eijk` function.

Return the A score ref <https://www.researchgate.net/publication/225958476_Measuring_Agreement_in_Ordered_Rating_Scales>

## Limiations

The first half of the histogram should not not be too heavy (>2x) or way too
light (<0.5x) of the second half of the histogram.

# Usage

This crate provides two functions

- `van_der_eijk` return A score. If it is less than 0.0, the distribution is
  very likely bimodal.
- `is_histogram_bimodal` is a wrapper on `van_der_eijk` function and it returns
  `True` if histogram is bimodal.


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

## TODO

- [ ] python bindings
