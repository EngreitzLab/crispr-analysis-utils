# Guide alignment

`cau guide-alignment` finds where a CRISPR guide library maps in a genome. A
guide is designed, not sequenced, so its site is an **exact** match of the
protospacer: the tool scans the reference once and answers three questions per
guide — where does the protospacer match, is there a PAM 3' of it, and what
sits 5' of it. There is no aligner, no index and no error budget.

A 44,000-guide library against hg38 takes about half a minute.

In Claude Code the `cau-guide-alignment` skill teaches Claude to run it and
read the outputs ([Agents and skills](agents-and-skills.md)).

## Quick start

```bash
cau guide-alignment run \
    --guides guides.tsv \
    --reference GRCh38.fa \
    --outdir results/
```

The guide table is two columns by default, a name and a sequence, with no
header:

```text
240911_EC_TAPseq_Enhancers_1    GAAAAGCCAACATGAATGCAG
240911_EC_TAPseq_Enhancers_2    GTTTGAGCACATAGCTCTAAG
```

(the separator above is a tab)

The reference must be an uncompressed FASTA. Nothing is written outside
`--outdir`, and nothing is written beside your inputs.

## What the sequences must be

The guide as ordered, 5' to 3', **without** the PAM. The tool appends `--pam`
(default `NGG`) itself.

`--spacer-length` (default 20) says how long the protospacer is. Anything 5'
of it — typically the G added for the U6 promoter — is aligned too, but a
mismatch there is tolerated and reported rather than counted. That is what
makes the two sequence columns of the output worth comparing.

If your table holds bare protospacers and the construct adds a G that is not
in the file, pass `--add-leading-g`.

!!! warning "One spacer in sixteen ends in `GG`"

    So a table that silently includes the PAM will run without complaint and
    be wrong by three bases. Check the column before the first run.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--guides`, `--reference`, `--outdir` | required | |
| `--pam` | `NGG` | IUPAC codes, 3' of the protospacer |
| `--spacer-length` | `20` | The protospacer is the 3'-most this many bases |
| `--add-leading-g` | off | Put a G in front of every sequence |
| `--orientation-check` | `error` | `error`, `warn` or `off` |
| `--fix-orientation` | off | Align the reverse complement of flagged entries |
| `--contigs` | every contig | Comma-separated names |
| `--id-col`, `--spacer-col` | `0`, `1` | A 0-based position, or a column name |
| `--sep` | tab | `tab`, `comma`, or one character |
| `--header` | off | The first line holds column names |

For hg38, give `--contigs` the primary chromosomes: without it the alt and
random scaffolds are read too, and a guide will match its own alt copy.

## Coordinates

0-based and half-open, on the reference's forward strand, and the protospacer
alone — the PAM and any base 5' of it are outside `[start, end)`.

`cut_site` is the boundary the blunt SpCas9 cut falls on, three bases from the
PAM: the cut lies between `cut_site - 1` and `cut_site`. On `+` that is
`end - 3`, on `-` it is `start + 3`.

## Orientation

Both readings of every guide go into the same scan. An entry is
`reverse_complement` only when its reverse complement has a PAM-bearing site
and the sequence as given does not. Reading a `CCN` 5' of a match is **not**
the test — that falls there by chance one time in sixteen.

A table reversed throughout stops the run. A handful of entries only warns,
because a single guide can look reversed when its PAM simply differs between
assemblies. `orientation_evidence` records whether the test was `distinct` or
`ambiguous`; ambiguous means the guide has nothing 5' of its spacer, so a flip
and a chance `CCN` cannot be told apart.

## Outputs

--8<-- "src/crispr_analysis_utils/skills/cau-guide-alignment/references/outputs.md:tables"

## What this does not do

**Off-target search.** Matching is exact, so there are no mismatch counts.
Where a guide can also cut with one or two mismatches is a different question
from where it comes from, and it is not answered here.

**Bulges, scores and other nucleases.** No DNA or RNA bulges, no CFD or MIT
scores, and 3' PAMs only.

## Coming from the GEM version

--8<-- "src/crispr_analysis_utils/skills/cau-guide-alignment/references/migration.md:migration"
