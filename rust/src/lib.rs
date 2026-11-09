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

use scan::{matches_pam, pam_length, Spacers};

/// The crate version, so Python can check which core it loaded.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// One exact match, as Python receives it.
type PyHit = (u32, String, usize, &'static str, String, String, bool);

/// One spacer's tally: matches found, how many ran off a contig end, and how
/// many of the rest had no PAM. Kept whether or not the matches are returned.
type PyCounts = (u64, u64, u64);

/// Every exact match of `spacers` in a reference FASTA.
///
/// Returns the matches and, beside them, one `(found, off_end, no_pam)` per
/// spacer. Each match is a tuple: the spacer's index in `spacers`, the
/// contig, the 0-based start of the protospacer on the forward strand, the
/// strand the guide reads on, the genomic PAM in guide orientation, the
/// `upstream` bases 5' of the spacer in guide orientation, and whether the
/// PAM matches.
///
/// A flank that would run off the end of its contig comes back empty, and a
/// PAM that does is not a match.
///
/// `pam_only` returns only the matches whose PAM matches. A short spacer
/// matches far more places than a nuclease could act at, and those are what
/// fills memory before Python ever sees them; the counts keep them, so
/// nothing the caller reports has to change.
#[pyfunction]
#[pyo3(signature = (fasta, spacers, pam = "NGG", upstream = 1, contigs = None, pam_only = false))]
fn scan_fasta(
    py: Python<'_>,
    fasta: &str,
    spacers: Vec<String>,
    pam: &str,
    upstream: usize,
    contigs: Option<Vec<String>>,
    pam_only: bool,
) -> PyResult<(Vec<PyHit>, Vec<PyCounts>)> {
    let width = pam_length(pam).map_err(PyValueError::new_err)?;
    let table = Spacers::new(&spacers).map_err(PyValueError::new_err)?;
    let wanted: Option<HashSet<String>> = contigs.map(|names| names.into_iter().collect());
    let mut seen: HashSet<String> = HashSet::new();

    let file = File::open(fasta)
        .map_err(|error| PyOSError::new_err(format!("Cannot read {fasta}: {error}")))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);

    let mut hits = Vec::new();
    let mut counts: Vec<PyCounts> = vec![(0, 0, 0); spacers.len()];
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
                    let tally = &mut counts[hit.guide as usize];
                    tally.0 += 1;
                    if hit.pam.len() != width {
                        tally.1 += 1;
                    } else if !matched {
                        tally.2 += 1;
                    }
                    if pam_only && !matched {
                        continue;
                    }
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
                    let mut missing: Vec<&str> = set
                        .iter()
                        .filter(|n| !seen.contains(*n))
                        .map(|n| n.as_str())
                        .collect();
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
    Ok((hits, counts))
}

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", version())?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_function(wrap_pyfunction!(scan_fasta, module)?)?;
    Ok(())
}
