# Guide alignment outputs

<!-- --8<-- [start:tables] -->
Everything goes into the output folder (`--outdir`, the `outdir` argument of
`run`), which is created when missing; nothing is written anywhere else. The
tables are tab-separated with a header line, except the BED files. A value
that does not apply is an empty field, and booleans are `TRUE` or `FALSE`.

Coordinates are 0-based and half-open on the reference's forward strand, as
in BED, and the protospacer excludes the PAM and the added G. Sequences
(spacers, protospacers, PAMs, flanks) read 5' to 3' in guide orientation:
the reverse complement of the forward strand for a `-` site.

## `guides.tsv`

One row per guide, in the order of the guide table, duplicates included. The
"Name" column says which names follow the lab's IGVF guide metadata and which
are this tool's own; the IGVF names keep this tool's definitions, given here.

| Column | Name | Content |
| --- | --- | --- |
| `guide_id` | IGVF | The id from the guide table, as text |
| `spacer` | IGVF | The spacer: the 3'-most `--spacer-length` bases of the input sequence |
| `pam` | IGVF | The genomic PAM of the unique site, such as `TGG` |
| `targeting` | IGVF | `TRUE` when the guide has at least one perfect primary-PAM site (class `unique` or `multi`) |
| `guide_chr` | IGVF | The unique site's contig |
| `guide_start` | IGVF | The unique site's protospacer start |
| `guide_end` | IGVF | The unique site's protospacer end, exclusive |
| `strand` | IGVF | The unique site's strand: `+` when the protospacer reads 5' to 3' on the forward strand, else `-` |
| `cut_site` | tool | The unique site's cut boundary `b`: the cut falls between bases `b - 1` and `b`, 3 bases 5' of the PAM, so `b` is `guide_end - 3` on `+` and `guide_start + 3` on `-` |
| `alignment_class` | tool | `unique`, `multi`, `off_target_only` or `no_site` (below) |
| `orientation` | tool | The orientation verdict (see `orientation.tsv`), or `not_checked` when the check was off |
| `leading_g` | tool | `TRUE` when the guide's site reads carry the added 5' G |
| `min_ngg_mismatches` | tool | The fewest spacer mismatches among the guide's primary-PAM sites; empty without any |
| `n_ngg_0mm` ... `n_ngg_3mm` | tool | The guide's primary-PAM sites with 0, 1, 2 and 3 spacer mismatches |
| `n_alt_pam_0mm` ... `n_alt_pam_3mm` | tool | Its alternative-PAM sites, counted the same way |
| `duplicate_of` | tool | The id of the first guide with the same spacer; empty for that first guide |
| `input_sequence` | tool | The sequence as given, stripped of surrounding whitespace and uppercased |

`pam`, `guide_chr`, `guide_start`, `guide_end`, `strand` and `cut_site` are
filled only for `unique` guides. The count columns run from `0mm` to
`--max-mismatches` (3 by default). `n_ngg_*` and `min_ngg_mismatches` count
the primary PAM's sites whatever `--pam` is.

The class is decided on primary-PAM sites only:

| `alignment_class` | Primary-PAM sites |
| --- | --- |
| `unique` | Exactly one with 0 spacer mismatches |
| `multi` | More than one with 0 spacer mismatches |
| `off_target_only` | None with 0 mismatches, at least one with 1 to `--max-mismatches` |
| `no_site` | None within `--max-mismatches` |

Alternative-PAM sites are counted but never decide the class, and a site
whose genomic PAM matches no PAM pattern is not counted at all. A guide whose
spacer repeats an earlier guide's gets that guide's sites, class and verdict.

## `sites.tsv`

Every site kept for every guide: one row per site with at most
`--max-mismatches` spacer mismatches, whatever its PAM. Guides come in table
order, duplicates included under their own id; a guide's sites are sorted by
contig in reference order, then start, then strand. A site found by the reads
of several PAMs is listed once.

| Column | Content |
| --- | --- |
| `guide_id` | The guide's id |
| `chr` | The contig |
| `start` | The protospacer's first base |
| `end` | One past the protospacer's last base |
| `strand` | `+` or `-`: the strand the protospacer reads 5' to 3' on |
| `pam` | The genomic PAM, uppercased |
| `pam_class` | The first PAM pattern the genomic PAM matches, trying `--pam` and then `--alt-pams` in their order, or `none`. A genomic `N` matches only an `N` of the pattern |
| `n_mismatches` | Mismatches between the spacer and the genomic protospacer. A genomic `N` counts as one; the added G never counts |
| `mismatch_positions` | The 1-based spacer positions of those mismatches, counted from the guide's 5' end, comma-separated; empty for a perfect site |
| `genomic_protospacer` | The reference bases of the protospacer, uppercased |
| `leading_g_base` | For a guide with the added G, the genomic base under that G; empty otherwise |

Sites whose `pam_class` is `none` are the ones GEM reported within its error
budget for the PAM reads, so they are not a complete list of such sites, and
they count nowhere in `guides.tsv`.

## `orientation.tsv`

The orientation check (pass 1), written unless `--orientation-check off`.
Each guide, duplicates included, has one row per perfect hit of its bare
spacer, or a single row with empty hit fields when it has none; `verdict`
repeats on every row of a guide.

| Column | Content |
| --- | --- |
| `guide_id` | The guide's id |
| `chr`, `start`, `end`, `strand` | The perfect hit, as in `sites.tsv` |
| `flank_5p` | The PAM-length stretch of reference 5' of the hit, uppercased; shorter where the contig ends |
| `flank_3p` | The same on the hit's 3' side |
| `verdict` | The guide's verdict, below |

| `verdict` | Meaning |
| --- | --- |
| `ok` | A perfect hit has a PAM (`--pam` or one of `--alt-pams`) on its 3' side |
| `reversed` | No hit is `ok`, but one has the reverse complement of a PAM on its 5' side: the sequence is probably given reverse-complemented |
| `no_pam` | Perfect hits, none of them `ok` or `reversed` |
| `no_perfect_hit` | No perfect hit |

A flank cut short by the end of a contig matches no PAM.

## `guides.bed` and `cut_sites.bed`

BED6 without a header, one row per `unique` guide (a duplicate under its own
id): contig, start, end, guide id, score `0`, strand. Rows are sorted by
contig name as plain text, then by position and guide id.

- `guides.bed`: the protospacer, `guide_start` to `guide_end`.
- `cut_sites.bed`: the one base right of the cut on the forward strand,
  `cut_site` to `cut_site + 1`.

## `summary.tsv`

Three columns, `section`, `name` and `value`, with the sections in this order:

| `section` | Rows |
| --- | --- |
| `guides` | `n_guides`; `n_distinct_spacers`; `n_duplicates`, the guides with `duplicate_of` set |
| `class` | One row per `alignment_class`; they sum to `n_guides` |
| `orientation` | One row per verdict, `not_checked` included; they sum to `n_guides` |
| `orientation_check` | The alignment counts of pass 1, below; absent when the check is off |
| `site_search` | The alignment counts of pass 2, below |
| `parameter` | The run's parameters, then `guides`, `reference`, `index`, `gem_version` and `cau_version` |

The alignment counts:

| Name | Alignments |
| --- | --- |
| `n_out_of_range` | gem-mapper records placed outside their contig (a read that overhangs a contig end), dropped before anything else |
| `n_records` | SAM records read |
| `n_unmapped` | Records of reads with no alignment |
| `n_imperfect` | Pass 1 only: alignments whose reference bases differ from the spacer, set aside |
| `n_outside_contigs` | Alignments on contigs that `--contigs` leaves out |
| `n_off_contig_end` | Alignments whose G, protospacer, PAM or bare spacer runs past the end of the contig |
| `n_over_max_mismatches` | Pass 2 only: alignments with more spacer mismatches than `--max-mismatches` once recounted against the reference, dropped |
| `n_repeated` | Alignments at a site already kept for the same guide, merged into it |
| `n_nm_disagreements` | Pass 2 only: alignments whose `NM` tag differs from the recount (spacer mismatches, plus each PAM position whose code is not A, C, G or T or differs from the genome, plus a mismatched added G) |
| `n_hits` | Pass 1 only: perfect hits kept |
| `n_sites` | Pass 2 only: sites kept, every PAM class included, before duplicates get theirs |

## `run.json`

The record of the run, rewritten after each GEM pass and at the end.

| Key | Content |
| --- | --- |
| `complete` | `true` once every output is written |
| `error` | The exception that stopped the run, as `Type: message`, or `null` |
| `started`, `finished` | UTC timestamps |
| `cau_version`, `gem_version` | The package and the gem-mapper versions |
| `inputs` | The guide table (path and MD5 checksum; `null` for a DataFrame), the reference (path, size, MD5 checksum), the index (path, size, modification time, and the reference checksum its `PREFIX.gem.json` records) and the FASTA index used |
| `parameters` | The parameters, as given |
| `passes` | Each GEM pass (`orientation`, `sites`, `sites_leading_g`) that ran or was reused: its FASTQ, SAM and log, the read count, the gem-mapper command, the wall time, the record counts, the `key` a rerun compares, and `reused` |
| `counts` | The alignment counts of `summary.tsv`, by pass |
| `classes`, `orientation` | The class and verdict counts of `summary.tsv` |
| `seconds` | The wall time of each step, and `total` |
| `alignments_per_guide` | For each guide that is not a duplicate, the alignment records of its reads kept after the `n_out_of_range` drop: `orientation` (when the check ran) and `sites`, both site passes together |

A rerun into the same folder reuses a GEM pass when the reads, the GEM
options, the index and the GEM version recorded for it match, and its SAM is
unchanged. Until a run completes, the records of the passes it has not reached
are kept from the earlier run, marked `carried_over`.

## The other files

| Path | Content |
| --- | --- |
| `alignments/orientation.fastq`, `.sam.gz` | Pass 1: one read per distinct spacer, named by its guide's row (0-based, in table order), and its alignments |
| `alignments/sites.fastq`, `.sam.gz` | Pass 2, the reads without an added G: `spacer + PAM`, one per distinct spacer and PAM pattern (the PAM written as the pattern, such as `NGG`), named `<row>:<PAM>` |
| `alignments/sites_leading_g.fastq`, `.sam.gz` | Pass 2, the reads with the added G: `G + spacer + PAM` |
| `logs/gem-mapper.<pass>.log` | The gem-mapper command, then its messages |
| `logs/guide_alignment.log` | The run's log, appended to by every run into the folder |
| `<reference name>.fai` | The FASTA index, built here when the reference has none next to it |

A pass with no reads (such as `sites_leading_g` without `--add-leading-g`)
writes no files.

## The index

`cau guide-alignment index REFERENCE PREFIX` (Python: `build_index`) writes
`PREFIX.gem`, gem-indexer's `PREFIX.info`, its output in `PREFIX.log`, and
`PREFIX.gem.json`: the reference's path, size and MD5 checksum, the GEM
version, the command and its wall time. `run` compares the reference's MD5
checksum with this file and stops on a mismatch.

## The Python return value

`run` returns a `RunResult`:

| Attribute | Content |
| --- | --- |
| `outdir` | The output folder |
| `guides` | The library: one `Guide` per row of the table |
| `summaries` | One `GuideSummary` per guide, in table order: the fields of `guides.tsv`, with the counts as the tuples `n_ngg` and `n_alt_pam` |
| `sites` | A `SiteResult`: the kept sites by guide row (`.sites`), and the counts of `site_search` but `n_out_of_range` (`.counts()`) |
| `orientation` | An `OrientationResult` (perfect hits and verdicts by guide row), or `None` when the check was off |
| `passes` | The `passes` record of `run.json` |
<!-- --8<-- [end:tables] -->
