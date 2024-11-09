use pyo3::prelude::*;

#[pyfunction]
fn van_der_eijk(histogram: Vec<usize>) -> f64 {
    ::is_bimodal::van_der_eijk(histogram)
}

#[pyfunction]
fn is_histogram_bimodal(histogram: Vec<usize>) -> bool {
    ::is_bimodal::is_histogram_bimodal(histogram)
}

#[pymodule]
fn is_bimodal(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(van_der_eijk, m)?)?;
    m.add_function(wrap_pyfunction!(is_histogram_bimodal, m)?)?;
    Ok(())
}
