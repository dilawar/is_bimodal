//! Check bimodality of a distribution using `van_der_eijk` function.

use pyo3::prelude::*;

/// A Python module implemented in Rust.
#[pymodule]
fn is_bimodal(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(is_histogram_bimodal, m)?)?;
    m.add_function(wrap_pyfunction!(van_der_eijk, m)?)?;
    Ok(())
}

/// Return the A score ref <https://www.researchgate.net/publication/225958476_Measuring_Agreement_in_Ordered_Rating_Scales>
///
/// # Limiations
///
/// 1. The first half of the histogram should not not be way too heavy (>2x) or
///    way too light (<0.5x) of the second half of the histogram.
#[pyfunction]
pub fn van_der_eijk(histogram: Vec<u32>) -> f64 {
    let mut a_score = 0.0;
    // Ensure that minimum value is 0.
    let min_value = histogram.iter().min().expect("No min value?");
    let mut layer: Vec<_> = histogram.iter().map(|x| x - min_value).collect();
    let total = histogram.iter().sum::<u32>() as f64;

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
        // now
        let a = compute_a(&layer_bin);
        a_score += weight * a;
        // tracing::trace!("> {n_min}, weight={weight:.3} a={a}:
        // layer={layer_bin:?}");
    }
    tracing::debug!("a_score={a_score} for sequence {histogram:?}");
    a_score
}

/// Check if given histogram is bimodal.
///
/// If `A` score is negative then the histogram is very likely to be bimodal.
#[pyfunction]
pub fn is_histogram_bimodal(histogram: Vec<u32>) -> bool {
    van_der_eijk(histogram) <= 0.0
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
    let s = layer.iter().filter(|e| **e).count() as f64;
    if s == layer.len() as f64 {
        return 1.0;
    }

    // tracing::trace!("> layer={layer:?} s={s}");
    // Paper doesn't say anything about this case. I discard this by returning 0.0
    if s <= 1.0 {
        return 0.0;
    }

    let tdu = count_triple(layer, [true, false, true]) as f64;
    let tu = count_triple(layer, [true, true, false]) as f64
        + count_triple(layer, [false, true, true]) as f64;
    let k = layer.len() as f64;
    // tracing::trace!("{layer}: k={k} s={s} tu={tu} tdu={tdu}");

    let num = (k - 2.0) * tu - (k - 1.0) * tdu;
    let den = (k - 2.0) * (tu + tdu);
    let u = num / den;
    // tracing::trace!("> u={u}={num}/{den}, k={k}, s={s}, tu={tu}, tdu={tdu}.");

    u * (1.0 - (s - 1.0) / (k - 1.0))
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

/// Count number of triples found in layer.
fn count_triple(layer: &[bool], triple: [bool; 3]) -> usize {
    let mut count = 0usize;
    let mut anyfind = false;
    for (ai, a) in layer.iter().enumerate() {
        if *a == triple[0] {
            for (bi, b) in layer.iter().skip(ai + 1).enumerate() {
                if *b == triple[1] {
                    for c in layer.iter().skip(ai + bi + 2) {
                        if *c == triple[2] {
                            count += 1;
                            anyfind = true;
                        }
                    }
                    // if no 'c' was found for a 'b', there is no 'c' for any other b.
                    if !anyfind {
                        // eprintln!("No c was found for this b, quitting...");
                        break;
                    }
                    anyfind = true;
                }
            }
            // if no 'b' was found for a 'a', there is no 'b' for any other a.
            if !anyfind {
                // eprintln!("No b was found for this a, quitting...");
                break;
            }
        }
    }
    count
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

        // fixed versions of above tests.
        assert!(van_der_eijk(vec![10, 11, 0, 0, 0, 0, 10, 2]) < 0.0);
        assert!(van_der_eijk(vec![10, 11, 0, 0, 0, 0, 20, 11]) < 0.0);
    }

    // These should pass but currently failing.
    #[test]
    fn test_these_should_be_bimodal() {
        assert!(is_histogram_bimodal(vec![4, 1, 1, 4, 1]));
        assert!(is_histogram_bimodal(vec![1, 4, 1, 1, 4, 1]));
    }

    fn compute_a_str(layer: &str) -> f64 {
        compute_a(&layer.chars().map(|x| x == '1').collect::<Vec<_>>())
    }

    fn count_triple_str(layer: &str, pattern: [bool; 3]) -> usize {
        count_triple(
            &layer.chars().map(|x| x == '1').collect::<Vec<_>>(),
            pattern,
        )
    }

    #[test]
    fn test_count_triples() {
        assert_eq!(count_triple_str("1001000", [true, false, true]), 2);
        assert_eq!(count_triple_str("1001000", [true, true, false]), 3);
        assert_eq!(count_triple_str("1001000", [false, true, true]), 0);
    }

    #[traced_test]
    #[test]
    fn test_compute_a() {
        assert_float_eq!(1.000, compute_a_str("1111111"), abs <= 0.001);
        assert_float_eq!(0.833, compute_a_str("1100000"), abs <= 0.001);
        assert_float_eq!(0.467, compute_a_str("1010000"), abs <= 0.001);
        assert_float_eq!(0.100, compute_a_str("1001000"), abs <= 0.001);
        assert_float_eq!(-0.267, compute_a_str("1000100"), abs <= 0.001);
        assert_float_eq!(-0.633, compute_a_str("1000010"), abs <= 0.001);
        assert_float_eq!(-1.000, compute_a_str("1000001"), abs <= 0.001);
    }
}
