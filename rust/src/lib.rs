//! The compiled core of crispr-analysis-utils.
//!
//! Imported as `crispr_analysis_utils._core`. It holds the exact guide scan: a
//! guide is designed, not sequenced, so its own site is an exact match of the
//! spacer, and the only other questions are what sits 3' of it (the PAM) and
//! what sits 5' of it (the added G).

mod scan;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;

use scan::{Spacers, matches_pam, pam_length};

/// The crate version, so Python can check which core it loaded.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// One exact match, as Python receives it.
type PyHit = (u32, String, usize, &'static str, String, String, bool);

/// Every exact match of `spacers` in a reference FASTA.
///
/// Returns one tuple per match: the spacer's index in `spacers`, the contig,
/// the 0-based start of the protospacer on the forward strand, the strand the
/// guide reads on, the genomic PAM in guide orientation, the `upstream` bases
/// 5' of the spacer in guide orientation, and whether the PAM matches.
///
/// A flank that would run off the end of its contig comes back empty, and a
/// PAM that does is not a match.
#[pyfunction]
#[pyo3(signature = (fasta, spacers, pam = "NGG", upstream = 1, contigs = None))]
fn scan_fasta(
    py: Python<'_>,
    fasta: &str,
    spacers: Vec<String>,
    pam: &str,
    upstream: usize,
    contigs: Option<Vec<String>>,
) -> PyResult<Vec<PyHit>> {
    let width = pam_length(pam).map_err(PyValueError::new_err)?;
    let table = Spacers::new(&spacers).map_err(PyValueError::new_err)?;
    let wanted: Option<HashSet<String>> = contigs.map(|names| names.into_iter().collect());
    let mut seen: HashSet<String> = HashSet::new();

    let file = File::open(fasta)
        .map_err(|error| PyOSError::new_err(format!("Cannot read {fasta}: {error}")))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);

    let mut hits = Vec::new();
    let mut name = String::new();
    let mut sequence: Vec<u8> = Vec::new();
    let mut line = Vec::new();

    loop {
        line.clear();
        let read = reader
            .read_until(b'\n', &mut line)
            .map_err(|error| PyOSError::new_err(format!("Cannot read {fasta}: {error}")))?;
        let last = read == 0;
        if last || line.first() == Some(&b'>') {
            if !name.is_empty() {
                seen.insert(name.clone());
            }
            if !name.is_empty() && wanted.as_ref().is_none_or(|set| set.contains(&name)) {
                // Releasing the GIL here is what makes a whole genome bearable
                // from Python: the scan touches nothing the interpreter owns.
                let found = py.detach(|| scan::scan(&sequence, &table, width, upstream));
                for hit in found {
                    let matched = matches_pam(&hit.pam, pam);
                    hits.push((
                        hit.guide,
                        name.clone(),
                        hit.start,
                        if hit.reverse { "-" } else { "+" },
                        hit.pam,
                        hit.upstream,
                        matched,
                    ));
                }
            }
            if last {
                // A contig that was asked for and never appeared is a typo, not
                // an empty result: without this the whole run would quietly
                // report that no guide maps anywhere.
                if let Some(set) = wanted.as_ref() {
                    let mut missing: Vec<&str> =
                        set.iter().filter(|n| !seen.contains(*n)).map(|n| n.as_str()).collect();
                    if !missing.is_empty() {
                        missing.sort_unstable();
                        return Err(PyValueError::new_err(format!(
                            "{fasta} has no contig named {}.",
                            missing
                                .iter()
                                .map(|n| format!("{n:?}"))
                                .collect::<Vec<_>>()
                                .join(", ")
                        )));
                    }
                }
                break;
            }
            name = String::from_utf8_lossy(&line[1..])
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            sequence.clear();
            continue;
        }
        sequence.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
    }
    Ok(hits)
}

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", version())?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_function(wrap_pyfunction!(scan_fasta, module)?)?;
    Ok(())
}
