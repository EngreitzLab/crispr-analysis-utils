---
name: cau-qc
description: >-
  Draw quality-control figures for a CRISPR screen with crispr-analysis-utils
  (`crispr_analysis_utils.qc`), each published as a PDF, the numbers behind
  it and a caption. Use when a user wants a knee plot or other QC figure from
  an assembled MuData, asks how deep a channel is, wants a figure ready for a
  paper, or asks where the numbers behind a figure are. Covers the three
  files every figure writes, the house style, and why a figure is refused
  without its data.
---

# QC figures (`crispr_analysis_utils.qc`)

```python
import mudata as md
from crispr_analysis_utils import qc

data = md.read_h5mu("results/mudata/mudata.h5mu")
figure, numbers, caption = qc.knee(data)
qc.publish("results/qc/knee_umis_per_barcode", figure, numbers, caption)
```

## Every figure is three files

| File | Holds |
| --- | --- |
| `<stem>.pdf` | the figure at 300 dpi, text kept as text |
| `<stem>.tsv.gz` | the numbers it was drawn from, one row per plotted point |
| `<stem>.txt` | the caption, as it would read under the figure |

**The data file is the point.** A reviewer can redraw the figure and a reader
can find the number behind a claim without rerunning anything. `publish`
refuses an empty table or a blank caption, because a figure without either is
an assertion rather than a result, and failing at write time is cheaper than
finding out in review.

## The style is set once

`qc.style()` puts Okabe-Ito on the colour cycle, `cividis` on images, a grey
grid behind the data, and drops the top and right spines. A figure that names
no colour is already right; reach for `qc.OKABE_ITO` only to override one.
`pdf.fonttype 42` keeps text editable in the PDF, so a figure can be
relabelled without redrawing it.

## The knee plot

`qc.knee(data)` ranks every barcode by total UMIs, both axes logarithmic, one
curve per modality.

- **A barcode a modality does not hold counts as zero**, not as missing. The
  assembled object is the union of the two matrices, and dropping the
  absentees would make a channel look deeper than it is.
- **The caption carries the real counts**, generated with the figure, so it
  cannot drift from what was drawn.
- The knee is where real cells give way to ambient background. Its position
  is what a count threshold has to justify itself against; this step draws
  it and sets no threshold.

## The count floor

`qc.filter_cells(data, min_umis=..., modality="gene")` keeps the barcodes
whose UMI count reaches the floor, returning the filtered copy and the mask.

- **`min_umis` has no default.** Where the knee falls depends on the
  channel's depth, so read it off that channel's own knee plot. A floor
  carried over from another channel is a guess.
- **The floor is applied to one modality**, `gene` by default. A barcode with
  real RNA and no guide counts is kept: an unassigned cell is a result, not a
  failure.
- Fractions computed before filtering describe the whole barcode space, which
  is mostly ambient. A mitochondrial median of 0 on an unfiltered object means
  the median barcode holds almost no UMIs, not that the channel is clean.

## Per-cell fractions

`qc.cell_fractions(data)` returns `total_counts`, `n_genes_by_counts` and a
`pct_counts_<class>` column per class, over `scanpy`'s
`calculate_qc_metrics`, so the names match the IGVF pipeline's own.

`qc.gene_classes(var)` decides the classes:

| Class | From |
| --- | --- |
| `mito` | `gene_chr == "chrM"` |
| `ribo` | an `RPS`/`RPL` symbol on a protein-coding gene |
| `ribo_pseudogene` | the same prefixes on anything else |
| `targeted` | the `readout` column, present only for a targeted panel |

- **Mitochondrial means the contig, not the symbol.** The `MT` prefix also
  names the nuclear metallothioneins and the `MTRNR2L` family, which are not
  mitochondrial. On GENCODE v43 the contig gives 37 genes and the prefix 566.
- **Ribosomal pseudogenes are counted apart**, since they share the `RPS` and
  `RPL` prefixes with the ~80 real ribosomal proteins and are transcribed.
- A cell the modality does not hold gets `NaN`, not zero: it was not
  measured.

## The metric distributions

`qc.metric_distributions(frame, nmads=5)` draws one histogram per column
of a `cell_fractions` frame, with matplotlib directly rather than through
scanpy's plotting.

- **Every panel thresholds on its own metric.** It calls `qc.mad_bounds`,
  draws the result as a dashed black line, and colours by it: blue
  (`qc.KEPT`) within the band, grey (`qc.REMOVED`) beyond it. The line and
  the cut are one computation, so the figure cannot disagree with the
  table.
- **Counts are two-sided and go on a log axis**, fractions are upper-tail
  only and linear (`LOG_METRICS`, overridable with `sides`).
- **The grey bars do not sum to the cells a filter drops.** Each panel
  counts its own rejects, so a cell that two metrics reject is grey twice.
  `mad_outlier` per metric, combined with `any`, is what gives the cells.
- A metric with no spread loses its threshold, not the figure: the panel is
  drawn with no line and the caption says so.
- A cell the modality does not hold is dropped, never binned as zero.

## Two metrics at once

`qc.metric_scatter(frame, x, y, x_threshold=..., y_threshold=...)` draws one
point per cell with the thresholds as dashed lines, each metric's histogram
on its margin, and the quadrant counts in the caption.

- **The thresholds are given, not derived.** A MAD band is one way to reach
  them; a library someone has QCed closely is another, and a fixed number
  chosen that way is a legitimate input.
- **Blue is at or above every threshold**, grey is everything else, so the
  kept quadrant reads at a glance.
- **The points are rasterised** and the text is not. One vector point per
  cell makes a PDF that will not open: the same figure in R came to 39 MB.
- A cell missing either metric is dropped, not drawn at zero.

## Thresholds from the channel itself

`qc.mad_outlier(values, nmads=5, side="both", log=False)` flags the cells
more than `nmads` median absolute deviations from the median. The spread is
the channel's own, so the threshold moves with its depth rather than being
carried over.

- **The MAD is rescaled to a sigma** (`scale="normal"`, a factor of 1.4826),
  so `nmads` reads as standard deviations. A raw MAD would make the same
  number cut a third more tightly.
- **`log=True` for a count**, which is roughly log-normal; off for a
  percentage. The base does not matter: a median and a MAD rescale together,
  so the band in the original units is the same either way.

- **It cannot pick the first cut.** The MAD describes whatever population
  it is given, and on an unfiltered barcode space that population is ambient
  droplets. On one real channel the median barcode held 2 UMIs and the
  5-MAD band on `log10` reached 131, so every actual cell was an outlier.
  The order is: knee for the floor, then MAD on the cells that survive it.
- **`log=True` for a count, always.** On the raw scale the band is symmetric
  around the median and its lower bound goes negative, so only the upper
  tail ever cuts anything.
- **`side="upper"` for a contamination fraction.** A cell with unusually
  little mitochondrial signal is not a problem; one with a lot is.
- 5 MADs is the usual starting point. 3 is aggressive, and on a healthy
  channel cuts cells that are only deep.

## Ask before running

- **Which object?** One channel's MuData, or several concatenated. A knee
  plot over pooled channels hides a channel that failed.
- **Where do the figures go?** Nothing is written outside the stem's folder.
