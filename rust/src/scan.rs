//! The exact guide scan: one pass over a reference, no alignment.
//!
//! A guide is designed, not sequenced, so its own site is an exact match. The
//! only questions are where the spacer matches, what sits 3' of it (the PAM)
//! and what sits 5' of it (the added G). Both orientations of every spacer go
//! into one table, so a single forward pass answers all three.

use std::collections::HashMap;

/// Base codes: 0-3 for ACGT, 4 for anything else.
const INVALID: u8 = 4;

fn code(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => INVALID,
    }
}

fn base_of(c: u8) -> u8 {
    b"ACGT"[c as usize]
}

/// Which bases an IUPAC code stands for, as a 4-bit mask over ACGT.
fn iupac_mask(symbol: u8) -> Option<u8> {
    Some(match symbol.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => return None,
    })
}

/// One guide's spacer, in the orientation it was put into the table.
#[derive(Clone, Copy)]
struct Entry {
    guide: u32,
    /// True when the key is the reverse complement, so the guide reads on `-`.
    reverse: bool,
}

/// Every spacer and its reverse complement, keyed by the 2-bit packed sequence.
///
/// A bitmap over a hash of the key rejects almost every genomic position
/// without touching the map, which is what keeps the inner loop cheap.
pub struct Spacers {
    length: usize,
    table: HashMap<u64, Vec<Entry>>,
    filter: Vec<u64>,
    mask: u64,
}

const FILTER_BITS: u32 = 24;

impl Spacers {
    pub fn new(spacers: &[String]) -> Result<Self, String> {
        let length = match spacers.first() {
            Some(first) => first.len(),
            None => return Err("No spacers were given.".to_string()),
        };
        if length == 0 || length > 31 {
            return Err(format!("A spacer must be 1 to 31 bases, not {length}."));
        }
        let slots = 1usize << FILTER_BITS;
        let mut found = Spacers {
            length,
            table: HashMap::with_capacity(spacers.len() * 2),
            filter: vec![0u64; slots / 64],
            mask: (1u64 << (2 * length)) - 1,
        };
        for (index, spacer) in spacers.iter().enumerate() {
            if spacer.len() != length {
                return Err(format!(
                    "Every spacer must be the same length: {:?} is {} bases, not {}.",
                    spacer,
                    spacer.len(),
                    length
                ));
            }
            let mut forward = 0u64;
            let mut reverse = 0u64;
            for (offset, byte) in spacer.bytes().enumerate() {
                let c = code(byte);
                if c == INVALID {
                    return Err(format!("The spacer {spacer:?} holds a base outside ACGT."));
                }
                forward = (forward << 2) | c as u64;
                // The complement of a 2-bit code is 3 - c, and the reverse
                // complement is built from the other end.
                reverse |= ((3 - c) as u64) << (2 * offset);
            }
            let guide = index as u32;
            found.insert(forward, Entry { guide, reverse: false });
            // A palindromic spacer keeps both entries on purpose: it really does
            // read 5' to 3' on either strand, and the two readings have
            // different flanks, so one can have a PAM where the other does not.
            found.insert(reverse, Entry { guide, reverse: true });
        }
        Ok(found)
    }

    fn insert(&mut self, key: u64, entry: Entry) {
        let slot = Self::slot(key);
        self.filter[slot / 64] |= 1u64 << (slot % 64);
        self.table.entry(key).or_default().push(entry);
    }

    #[inline(always)]
    fn slot(key: u64) -> usize {
        // Multiply-shift: the low bits of a packed spacer are far from random,
        // so mix before taking the index.
        let mixed = key.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        (mixed >> (64 - FILTER_BITS)) as usize
    }

    #[inline(always)]
    fn maybe(&self, key: u64) -> bool {
        let slot = Self::slot(key);
        self.filter[slot / 64] >> (slot % 64) & 1 == 1
    }

    #[inline(always)]
    fn get(&self, key: u64) -> Option<&Vec<Entry>> {
        if !self.maybe(key) {
            return None;
        }
        self.table.get(&key)
    }
}

/// One exact match of a spacer, with the bases around it in guide orientation.
pub struct Hit {
    pub guide: u32,
    pub start: usize,
    pub reverse: bool,
    pub pam: String,
    pub upstream: String,
}

/// Every exact match of `spacers` in `sequence`, with its flanks.
///
/// `start` is the 0-based start of the protospacer on the forward strand, so
/// it is the left end whichever strand the guide reads on. `pam` and
/// `upstream` are in guide orientation, and are empty when the flank would
/// run off the end of the contig.
pub fn scan(sequence: &[u8], spacers: &Spacers, pam_length: usize, upstream: usize) -> Vec<Hit> {
    let length = spacers.length;
    let mut hits = Vec::new();
    if sequence.len() < length {
        return hits;
    }
    let mut key = 0u64;
    let mut run = 0usize;
    for (index, &byte) in sequence.iter().enumerate() {
        let c = code(byte);
        if c == INVALID {
            run = 0;
            key = 0;
            continue;
        }
        key = ((key << 2) | c as u64) & spacers.mask;
        run += 1;
        if run < length {
            continue;
        }
        let Some(entries) = spacers.get(key) else {
            continue;
        };
        let start = index + 1 - length;
        for entry in entries {
            let (pam, up) = if entry.reverse {
                // The guide reads right to left here: its PAM is 5' of the
                // match on the forward strand, reverse-complemented.
                (
                    revcomp_slice(sequence, start.wrapping_sub(pam_length), pam_length, start),
                    revcomp_slice(sequence, start + length, upstream, sequence.len()),
                )
            } else {
                (
                    plain_slice(sequence, start + length, pam_length, sequence.len()),
                    plain_slice(sequence, start.wrapping_sub(upstream), upstream, start),
                )
            };
            hits.push(Hit {
                guide: entry.guide,
                start,
                reverse: entry.reverse,
                pam,
                upstream: up,
            });
        }
    }
    hits
}

/// `length` bases from `from`, uppercased, or empty if they do not all exist.
fn plain_slice(sequence: &[u8], from: usize, length: usize, limit: usize) -> String {
    if length == 0 {
        return String::new();
    }
    if from > limit || from.saturating_add(length) > limit || from.saturating_add(length) > sequence.len() {
        return String::new();
    }
    sequence[from..from + length]
        .iter()
        .map(|&b| {
            let c = code(b);
            if c == INVALID {
                b'N'
            } else {
                base_of(c)
            }
        })
        .map(char::from)
        .collect()
}

/// The same bases, reverse-complemented, for a guide that reads on `-`.
fn revcomp_slice(sequence: &[u8], from: usize, length: usize, limit: usize) -> String {
    let plain = plain_slice(sequence, from, length, limit);
    if plain.is_empty() {
        return plain;
    }
    plain
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

/// Check a PAM pattern and return its length.
pub fn pam_length(pattern: &str) -> Result<usize, String> {
    if pattern.is_empty() {
        return Err("The PAM must not be empty.".to_string());
    }
    for symbol in pattern.bytes() {
        if iupac_mask(symbol).is_none() {
            return Err(format!(
                "The PAM {pattern:?} holds {:?}, which is not an IUPAC code.",
                symbol as char
            ));
        }
    }
    Ok(pattern.len())
}

/// Whether a stretch of genomic bases matches an IUPAC pattern.
pub fn matches_pam(bases: &str, pattern: &str) -> bool {
    if bases.len() != pattern.len() {
        return false;
    }
    bases.bytes().zip(pattern.bytes()).all(|(base, symbol)| {
        let c = code(base);
        match (c, iupac_mask(symbol)) {
            (INVALID, _) | (_, None) => false,
            (c, Some(mask)) => mask >> c & 1 == 1,
        }
    })
}
