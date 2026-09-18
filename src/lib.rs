//! Check bimodality of a distribution using the `van_der_eijk` A-score and the
//! bimodality coefficient.

// PyO3's `#[pyfunction]`/`#[pymodule]` macro expansion triggers a
// false-positive `clippy::useless_conversion` warning on the generated
// wrappers, which cannot be silenced on the individual functions.
#![allow(clippy::useless_conversion)]

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
fn is_bimodal(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_is_histogram_bimodal, m)?)?;
    m.add_function(wrap_pyfunction!(py_van_der_eijk, m)?)?;
    m.add_function(wrap_pyfunction!(py_bimodality_coefficient, m)?)?;
    Ok(())
}

/// Return the A score ref <https://www.researchgate.net/publication/225958476_Measuring_Agreement_in_Ordered_Rating_Scales>
///
/// The A-score lies in `[-1, 1]`: `1` is perfect unimodality, `0` is perfect
/// uniformity and `-1` is perfect bimodality.
///
/// Returns `None` when the input is not a valid histogram: it must have at
/// least three bins and a positive total count.
///
/// # Limiations
///
/// 1. The first half of the histogram should not not be way too heavy (>2x) or
///    way too light (<0.5x) of the second half of the histogram.
pub fn van_der_eijk(histogram: &[u32]) -> Option<f64> {
    if histogram.len() < 3 {
        return None;
    }
    // Use u64 for the total to avoid overflow on large histograms.
    let total: u64 = histogram.iter().map(|&x| u64::from(x)).sum();
    if total == 0 {
        return None;
    }
    let total = total as f64;

    // Ensure that minimum value is 0.
    let min_value = *histogram.iter().min()?;
    let mut layer: Vec<_> = histogram.iter().map(|x| x - min_value).collect();
    let mut a_score = 0.0;

    while let Some(n_min) = non_zero_min(&layer) {
        let mut layer_bin = vec![false; layer.len()];
        let mut weight = 0.0;
        layer.iter_mut().enumerate().for_each(|(ei, e)| {
            if *e > 0 {
                layer_bin[ei] = true;
                weight += (n_min as f64) / total;
                *e -= n_min;
            }
        });
        let a = compute_a(&layer_bin);
        a_score += weight * a;
    }
    Some(a_score)
}

/// Check if given histogram is bimodal.
///
/// If the `A` score is negative then the histogram is very likely to be
/// bimodal. Returns `None` for invalid input (see [`van_der_eijk`]).
pub fn is_histogram_bimodal(histogram: &[u32]) -> Option<bool> {
    van_der_eijk(histogram).map(|a| a <= 0.0)
}

/// Return the (Sarle/SAS) bimodality coefficient.
///
/// ```text
///         skewness^2 + 1
/// BC = -------------------
///          kurtosis
/// ```
///
/// where `kurtosis` is the (non-excess) fourth standardized moment. Values
/// above `5/9 ≈ 0.555` suggest a bimodal distribution; values below suggest a
/// unimodal one. A constant histogram (zero variance) returns `0.0`.
///
/// Returns `None` when the input is not a valid histogram: it must have at
/// least three bins and a positive total count.
pub fn bimodality_coefficient(histogram: &[u32]) -> Option<f64> {
    if histogram.len() < 3 {
        return None;
    }
    let total: u64 = histogram.iter().map(|&x| u64::from(x)).sum();
    if total == 0 {
        return None;
    }
    let n = total as f64;

    // Treat the bin index as the value on the ordered scale.
    let mean = histogram
        .iter()
        .enumerate()
        .map(|(i, &f)| i as f64 * f as f64)
        .sum::<f64>()
        / n;

    let central = |power: i32| {
        histogram
            .iter()
            .enumerate()
            .map(|(i, &f)| f as f64 * (i as f64 - mean).powi(power))
            .sum::<f64>()
            / n
    };

    let m2 = central(2);
    if m2 == 0.0 {
        // All the mass sits on a single bin: perfectly unimodal.
        return Some(0.0);
    }
    let m3 = central(3);
    let m4 = central(4);

    let skewness = m3 / m2.powf(1.5);
    let kurtosis = m4 / m2.powi(2);
    Some((skewness * skewness + 1.0) / kurtosis)
}

/// Compute A score. The score is computed using following equations. See
/// reference [1].
///
/// ```ignore
///               S  - 1
/// A = U . (1 - -------)
///               K - 1
///
///      (K - 2) . TU - (K - 1) . TDU
/// U = -----------------------------
///          (K - 2) . (TU + TDU)
/// ```
///
/// Where,
///
/// K: layer length
/// S: No of '1'
/// TU: No of times 101 was seen.
/// TDU: no of times 110 and 011 was seen.
///
/// [1]. https://www.researchgate.net/publication/225958476_Measuring_Agreement_in_Ordered_Rating_Scales
fn compute_a(layer: &[bool]) -> f64 {
    let k = layer.len();
    let s = layer.iter().filter(|e| **e).count();
    // Every category is non-empty: the layer is perfectly uniform, A = 0.
    if s == k {
        return 0.0;
    }
    // A single non-empty category is perfect agreement, A = 1.
    // Ref: R `agrmt::patternAgreement`: `if (sum(P) == 1) (A <- 1)`.
    if s == 1 {
        return 1.0;
    }
    if s == 0 {
        return 0.0;
    }

    let (tu, tdu) = count_tu_tdu(layer);
    let k = k as f64;
    let num = (k - 2.0) * tu - (k - 1.0) * tdu;
    let den = (k - 2.0) * (tu + tdu);
    if den == 0.0 {
        return 0.0;
    }
    let u = num / den;

    u * (1.0 - (s as f64 - 1.0) / (k - 1.0))
}

/// Count `TU` (`110` and `011`) and `TDU` (`101`) subsequence triples in a
/// pattern in `O(n)` time using prefix counts.
///
/// For every position `j`:
///
/// - if `j` is a zero it can be the middle of a `101` triple, contributing
///   `ones_before[j] * ones_after[j]` to `TDU`,
/// - if `j` is a one it can be the middle of a `110` triple, contributing
///   `ones_before[j] * zeros_after[j]`, or of a `011` triple, contributing
///   `zeros_before[j] * ones_after[j]`, to `TU`.
fn count_tu_tdu(layer: &[bool]) -> (f64, f64) {
    let n = layer.len();
    let mut ones_before = vec![0usize; n];
    let mut zeros_before = vec![0usize; n];
    let (mut ones, mut zeros) = (0usize, 0usize);
    for (i, &b) in layer.iter().enumerate() {
        ones_before[i] = ones;
        zeros_before[i] = zeros;
        if b {
            ones += 1;
        } else {
            zeros += 1;
        }
    }

    let (mut ones_after, mut zeros_after) = (0usize, 0usize);
    let mut tu = 0usize;
    let mut tdu = 0usize;
    for i in (0..n).rev() {
        if layer[i] {
            tu += ones_before[i] * zeros_after + zeros_before[i] * ones_after;
            ones_after += 1;
        } else {
            tdu += ones_before[i] * ones_after;
            zeros_after += 1;
        }
    }

    (tu as f64, tdu as f64)
}

/// Find non-zero minimum element. Returns `Some(min)` if there is at least one
/// non-zero value, `None` otherwise.
#[inline]
fn non_zero_min(elements: &[u32]) -> Option<u32> {
    let mut min = u32::MAX;
    for e in elements {
        if *e == 0 {
            continue;
        }
        if *e < min {
            min = *e;
        }
    }
    if min == u32::MAX {
        None
    } else {
        Some(min)
    }
}

#[pyfunction]
#[pyo3(name = "van_der_eijk")]
fn py_van_der_eijk(histogram: Vec<u32>) -> PyResult<f64> {
    van_der_eijk(&histogram).ok_or_else(|| {
        PyValueError::new_err("histogram must have at least 3 bins and a positive total")
    })
}

#[pyfunction]
#[pyo3(name = "is_histogram_bimodal")]
fn py_is_histogram_bimodal(histogram: Vec<u32>) -> PyResult<bool> {
    is_histogram_bimodal(&histogram).ok_or_else(|| {
        PyValueError::new_err("histogram must have at least 3 bins and a positive total")
    })
}

#[pyfunction]
#[pyo3(name = "bimodality_coefficient")]
fn py_bimodality_coefficient(histogram: Vec<u32>) -> PyResult<f64> {
    bimodality_coefficient(&histogram).ok_or_else(|| {
        PyValueError::new_err("histogram must have at least 3 bins and a positive total")
    })
}

#[cfg(test)]
mod tests {
    use float_eq::assert_float_eq;
    use tracing_test::traced_test;

    use super::*;

    #[traced_test]
    #[test]
    fn test_van_der_eijk() {
        // unimodal and detected as unimodal.
        assert!(van_der_eijk(&[30, 40, 210, 130, 530, 50, 10]).unwrap() > 0.0);
        assert!(van_der_eijk(&[30, 40, 210, 10, 530, 50, 10]).unwrap() > 0.0);
        assert!(van_der_eijk(&[30, 40, 10, 10, 30, 50, 100]).unwrap() > 0.0);
        assert!(van_der_eijk(&[3, 4, 1, 1, 3, 5, 10]).unwrap() > 0.0);
        assert!(van_der_eijk(&[3, 4, 1, 1, 3, 5, 1]).unwrap() > 0.0);
        // perfectly uniform -> A == 0.
        assert_float_eq!(
            0.0,
            van_der_eijk(&[1, 1, 1, 1, 1, 1, 1]).unwrap(),
            abs <= 0.001
        );
        assert!(van_der_eijk(&[1, 1, 1, 1, 1, 1, 1000]).unwrap() > 0.0);

        // A single dominant peak at one end is unimodal, not bimodal.
        assert!(van_der_eijk(&[10000, 1, 1, 1, 1, 1, 10]).unwrap() > 0.0);

        // bimodal and detected as bimodal.
        assert!(van_der_eijk(&[10, 10, 0, 0, 0, 10, 10]).unwrap() < 0.0);
        assert!(van_der_eijk(&[10, 10, 0, 0, 0, 0, 10]).unwrap() < 0.0);
        assert!(van_der_eijk(&[1, 1, 1, 0, 0, 1, 1]).unwrap() < 0.0);
        assert!(van_der_eijk(&[1, 1, 1, 0, 1, 1, 1]).unwrap() < 0.0);

        // Test cases that bring the limitations of the algorithm.
        // This should be bi-modal. Algo fails because weights are not balanced here.
        // One side of the see-saw is 2x heavier.
        assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 3, 3]).unwrap() > 0.0);
        assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 30, 31]).unwrap() > 0.0);

        // fixed versions of above tests.
        assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 10, 2]).unwrap() < 0.0);
        assert!(van_der_eijk(&[10, 11, 0, 0, 0, 0, 20, 11]).unwrap() < 0.0);
    }

    /// The worked example from van der Eijk (2001, p. 331). The R package
    /// `agrmt` reports an A-score of 0.6113333 for it.
    #[test]
    fn test_reference_example() {
        assert_float_eq!(
            0.6113333,
            van_der_eijk(&[30, 40, 210, 130, 530, 50, 10]).unwrap(),
            abs <= 1e-6
        );
    }

    /// A perfectly uniform histogram has A == 0.
    #[test]
    fn test_uniform_histogram_has_zero_score() {
        assert_float_eq!(0.0, van_der_eijk(&[5, 5, 5, 5, 5]).unwrap(), abs <= 1e-9);
        assert_float_eq!(0.0, van_der_eijk(&[1, 1, 1]).unwrap(), abs <= 1e-9);
    }

    /// A layer with a single non-empty category has perfect agreement (A == 1),
    /// while a layer where every category is non-empty is uniform (A == 0).
    #[test]
    fn test_layer_edge_cases() {
        assert_float_eq!(1.0, compute_a_str("1000000"), abs <= 1e-9);
        assert_float_eq!(1.0, compute_a_str("0001000"), abs <= 1e-9);
        assert_float_eq!(0.0, compute_a_str("1111111"), abs <= 1e-9);
    }

    /// Peaks in the middle of the histogram are detected.
    #[test]
    fn test_middle_peaks_are_bimodal() {
        assert!(is_histogram_bimodal(&[4, 1, 1, 4, 1]).unwrap());
        assert!(is_histogram_bimodal(&[1, 4, 1, 1, 4, 1]).unwrap());
        assert!(is_histogram_bimodal(&[4, 1, 2, 4, 1]).unwrap());
        assert!(is_histogram_bimodal(&[4, 4, 1, 1, 4]).unwrap());
        assert!(is_histogram_bimodal(&[1, 4, 2, 1, 1, 4, 1]).unwrap());
        assert!(is_histogram_bimodal(&[1, 4, 1, 1, 1, 4]).unwrap());
    }

    /// The A-score is always within [-1, 1].
    #[test]
    fn test_a_score_is_bounded() {
        let histograms: Vec<Vec<u32>> = vec![
            vec![1, 1, 1, 1, 1],
            vec![0, 0, 10, 0, 0],
            vec![10, 0, 0, 0, 10],
            vec![1, 2, 3, 4, 5, 6, 7],
            vec![7, 6, 5, 4, 3, 2, 1],
            vec![100, 1, 1, 1, 100],
        ];
        for h in histograms {
            let a = van_der_eijk(&h).unwrap();
            assert!(
                (-1.0..=1.0).contains(&a),
                "A-score {a} out of range for {h:?}"
            );
        }
    }

    /// Invalid input is rejected rather than panicking or returning NaN.
    #[test]
    fn test_invalid_input_returns_none() {
        assert!(van_der_eijk(&[]).is_none());
        assert!(van_der_eijk(&[1, 2]).is_none());
        assert!(van_der_eijk(&[0, 0, 0]).is_none());
        assert!(is_histogram_bimodal(&[]).is_none());
        assert!(bimodality_coefficient(&[1, 2]).is_none());
        assert!(bimodality_coefficient(&[0, 0, 0, 0]).is_none());
    }

    /// A uniform histogram has the theoretical bimodality coefficient 5/9.
    #[test]
    fn test_bimodality_coefficient_uniform() {
        // For a discrete uniform histogram the coefficient is close to, but not
        // exactly, 5/9; it converges to 5/9 as the number of bins grows.
        assert_float_eq!(
            0.5882352941176471,
            bimodality_coefficient(&[1, 1, 1, 1, 1]).unwrap(),
            abs <= 1e-12
        );
        assert_float_eq!(
            5.0 / 9.0,
            bimodality_coefficient(&[1; 1000]).unwrap(),
            abs <= 1e-4
        );
        // A single spike is perfectly unimodal.
        assert_float_eq!(
            0.0,
            bimodality_coefficient(&[0, 0, 10, 0, 0]).unwrap(),
            abs <= 1e-9
        );
    }

    /// A two-spike histogram has a bimodality coefficient above 5/9.
    #[test]
    fn test_bimodality_coefficient_two_spikes() {
        assert!(bimodality_coefficient(&[10, 0, 0, 0, 10]).unwrap() > 5.0 / 9.0);
        assert!(bimodality_coefficient(&[10, 0, 0, 0, 10]).unwrap() > 0.9);
    }

    fn compute_a_str(layer: &str) -> f64 {
        compute_a(&layer.chars().map(|x| x == '1').collect::<Vec<_>>())
    }

    /// Brute-force reference for [`count_tu_tdu`].
    fn count_triple(layer: &[bool], triple: [bool; 3]) -> usize {
        let mut count = 0usize;
        for i in 0..layer.len() {
            if layer[i] != triple[0] {
                continue;
            }
            for j in i + 1..layer.len() {
                if layer[j] != triple[1] {
                    continue;
                }
                for &c in layer.iter().skip(j + 1) {
                    if c == triple[2] {
                        count += 1;
                    }
                }
            }
        }
        count
    }

    fn count_triple_str(layer: &str, pattern: [bool; 3]) -> usize {
        count_triple(
            &layer.chars().map(|x| x == '1').collect::<Vec<_>>(),
            pattern,
        )
    }

    fn count_tu_tdu_str(layer: &str) -> (f64, f64) {
        count_tu_tdu(&layer.chars().map(|x| x == '1').collect::<Vec<_>>())
    }

    #[test]
    fn test_count_triples() {
        assert_eq!(count_triple_str("1001000", [true, false, true]), 2);
        assert_eq!(count_triple_str("1001000", [true, true, false]), 3);
        assert_eq!(count_triple_str("1001000", [false, true, true]), 0);
        // 011 patterns.
        assert_eq!(count_triple_str("0011000", [false, true, true]), 2);
        assert_eq!(count_triple_str("0110000", [false, true, true]), 1);
        // 101 patterns.
        assert_eq!(count_triple_str("1010100", [true, false, true]), 4);
    }

    /// The O(n) counter must agree with the brute-force triple counter on every
    /// pattern up to length 10.
    #[test]
    fn test_count_tu_tdu_matches_brute_force() {
        for len in 3..=10 {
            for mask in 0u32..(1 << len) {
                let layer: Vec<bool> = (0..len).map(|i| (mask >> i) & 1 == 1).collect();
                let (tu, tdu) = count_tu_tdu(&layer);
                let expected_tu = (count_triple(&layer, [true, true, false])
                    + count_triple(&layer, [false, true, true]))
                    as f64;
                let expected_tdu = count_triple(&layer, [true, false, true]) as f64;
                assert_eq!((tu, tdu), (expected_tu, expected_tdu), "layer={layer:?}");
            }
        }
    }

    #[test]
    fn test_count_tu_tdu_str() {
        // 110 = 3, 011 = 0, 101 = 2.
        assert_eq!(count_tu_tdu_str("1001000"), (3.0, 2.0));
        // 110 = 3, 011 = 2, 101 = 0.
        assert_eq!(count_tu_tdu_str("0011000"), (5.0, 0.0));
    }

    #[traced_test]
    #[test]
    fn test_compute_a() {
        assert_float_eq!(0.000, compute_a_str("1111111"), abs <= 0.001);
        assert_float_eq!(0.833, compute_a_str("1100000"), abs <= 0.001);
        assert_float_eq!(0.467, compute_a_str("1010000"), abs <= 0.001);
        assert_float_eq!(0.100, compute_a_str("1001000"), abs <= 0.001);
        assert_float_eq!(-0.267, compute_a_str("1000100"), abs <= 0.001);
        assert_float_eq!(-0.633, compute_a_str("1000010"), abs <= 0.001);
        assert_float_eq!(-1.000, compute_a_str("1000001"), abs <= 0.001);
    }
}
