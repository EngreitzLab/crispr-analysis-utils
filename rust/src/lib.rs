//! The compiled core of crispr-analysis-utils.
//!
//! Imported as `crispr_analysis_utils._core`. It holds nothing but the build
//! identity so far; the guide scan lands here next.

use pyo3::prelude::*;

/// The crate version, so Python can check which core it loaded.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", version())?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    Ok(())
}
