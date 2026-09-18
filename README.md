Check bimodality of a histogram using the `van_der_eijk` A-score [1] or the
(Sarle/SAS) bimodality coefficient.

The A-score ranges from `-1` to `1`:

-   `1`: perfect unimodality (agreement),
-   `0`: perfect uniformity,
-   `-1`: perfect bimodality (maximum disagreement).

The bimodality coefficient is the moment-based quantity

```text
         skewness^2 + 1
BC = ----------------------
          kurtosis
```

Values above `5/9 ≈ 0.555` suggest bimodality, values below suggest
unimodality. It needs no layer decomposition, so it can complement the
A-score when the latter's balance assumption is violated.

# Limitations

-   The two halves of the histogram should be roughly balanced. If one half
    holds more than twice the mass of the other, the score is biased towards
    unimodality and a genuinely bimodal histogram can be missed.
-   Peaks that are not at the extremes can still be missed, e.g.
    `[1, 4, 1, 1, 4, 4, 1]` scores `0.10` and is classified as unimodal.
-   The score does not report *where* the modes are.
-   At least three categories are required.

# Usage

This library provides three functions for both Rust and Python.

-   `van_der_eijk` returns the A-score as described in [1]. If it is less than
    `0.0`, the distribution is very likely bimodal.
-   `is_histogram_bimodal` is a wrapper on `van_der_eijk` and returns `True`
    when the A-score is `<= 0.0`.
-   `bimodality_coefficient` returns the (Sarle/SAS) bimodality coefficient;
    values above `5/9` suggest bimodality.

Invalid input (fewer than three bins, or an all-zero histogram) is rejected:
the Python functions raise `ValueError` and the Rust functions return `None`.

## Python

Here are some runs on small histograms. This library typically performs much
better on larger histograms.

```bash
>>> is_bimodal.van_der_eijk([3, 1, 3])  # A-score (negative means bimodal)
-0.5714285714285714
>>> is_bimodal.is_histogram_bimodal([3, 1, 3])
True
>>> is_bimodal.is_histogram_bimodal([3, 2, 3])
True
>>> is_bimodal.is_histogram_bimodal([3, 2, 4])
True
>>> is_bimodal.is_histogram_bimodal([4, 2, 2, 4])
True
>>> is_bimodal.is_histogram_bimodal([4, 1, 2, 4])
True
>>> is_bimodal.is_histogram_bimodal([4, 1, 2, 4, 1])
True
```

Peaks in the middle of the histogram are detected as well:

```
>>> is_bimodal.is_histogram_bimodal([4, 1, 1, 4, 1])
True
>>> is_bimodal.is_histogram_bimodal([1, 4, 2, 1, 1, 4, 1])
True
>>> is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4, 1])
-0.03809523809523809
>>> is_bimodal.van_der_eijk([1, 4, 2, 1, 1, 4, 4])
-0.020915032679738564
>>> is_bimodal.van_der_eijk([1, 4, 1, 1, 1, 4])
-0.275
>>> is_bimodal.van_der_eijk([4, 4, 1, 1, 4])
-0.17857142857142858
>>> is_bimodal.van_der_eijk([4, 4, 1, 1, 4, 4])
-0.13333333333333333
```

This one is bimodal but is **missed** because its two halves are not balanced
(the first half is much heavier than the second):

```
>>> is_bimodal.is_histogram_bimodal([1, 4, 1, 1, 4, 4, 1])
False
>>> is_bimodal.van_der_eijk([1, 4, 1, 1, 4, 4, 1])
0.1
```

Note that a negative A-score means bimodal.

The bimodality coefficient can be used on its own:

```
>>> is_bimodal.bimodality_coefficient([1, 1, 1, 1, 1])  # uniform, ~5/9
0.5882352941176471
>>> is_bimodal.bimodality_coefficient([10, 0, 0, 0, 10])  # two spikes
1.0
>>> is_bimodal.bimodality_coefficient([0, 0, 10, 0, 0])  # single spike
0.0
```

Invalid input is rejected instead of panicking or returning `NaN`:

```
>>> is_bimodal.is_histogram_bimodal([0, 0, 0])
ValueError: histogram must have at least 3 bins and a positive total
```

## Rust

Here are some examples from the unit tests. The first one is the example from
van der Eijk's paper (p. 331), which should score `0.6113333`.

```rust
assert_float_eq!(0.6113333, van_der_eijk(&[30, 40, 210, 130, 530, 50, 10]).unwrap(), abs <= 1e-6);

// unimodal and detected as unimodal.
assert!(van_der_eijk(&[30, 40, 210, 10, 530, 50, 10]).unwrap() > 0.0);
assert!(van_der_eijk(&[30, 40, 10, 10, 30, 50, 100]).unwrap() > 0.0);
assert!(van_der_eijk(&[3, 4, 1, 1, 3, 5, 10]).unwrap() > 0.0);
assert!(van_der_eijk(&[3, 4, 1, 1, 3, 5, 1]).unwrap() > 0.0);
assert!(van_der_eijk(&[1, 1, 1, 1, 1, 1, 1000]).unwrap() > 0.0);

// perfectly uniform -> A == 0.
assert_float_eq!(0.0, van_der_eijk(&[1, 1, 1, 1, 1, 1, 1]).unwrap(), abs <= 0.001);

// bimodal and detected as bimodal.
assert!(van_der_eijk(&[10, 10, 0, 0, 0, 10, 10]).unwrap() < 0.0);
assert!(van_der_eijk(&[10, 10, 0, 0, 0, 0, 10]).unwrap() < 0.0);
assert!(van_der_eijk(&[1, 1, 1, 0, 0, 1, 1]).unwrap() < 0.0);
assert!(van_der_eijk(&[1, 1, 1, 0, 1, 1, 1]).unwrap() < 0.0);

// Test cases that bring out the limitations of the algorithm.
// These should be bimodal, but the algorithm fails because the weights are
// not balanced: one side of the see-saw is more than 2x heavier.
assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 3, 3]).unwrap() > 0.0);
assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 30, 31]).unwrap() > 0.0);
assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 20, 11]).unwrap() < 0.0);

// Invalid input is rejected.
assert!(van_der_eijk(&[0, 0, 0]).is_none());
```

# References

[1] Eijk, Cees. (2001). Measuring Agreement in Ordered Rating Scales. Quality
and Quantity. 35. 325-341. 10.1023/A:1010374114305.
