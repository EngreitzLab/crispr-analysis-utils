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

One pass over a channel: find the floor, measure the cells, cut, write.

```python
import mudata as md
from crispr_analysis_utils import qc

data = md.read_h5mu("results/mudata/mudata.h5mu")
qc.publish("results/qc/knee", *qc.knee(data))
cells, _ = qc.filter_cells(data, min_umis=1000)  # read off that knee
fractions = qc.cell_fractions(cells)
qc.publish("results/qc/metrics", *qc.metric_distributions(fractions))
passes, rules = qc.cell_filters(fractions, nmads=5)
clean = qc.apply_filters(cells, passes, rules, path="results/qc/filtered.h5mu")
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
- **`extra={name: mask}` adds a class of someone else's definition.** A
  threshold set against an older analysis's idea of ribosomal only applies
  to that definition, so measure it rather than assume the two agree.

## The metric distributions

`qc.metric_distributions(frame, nmads=5)` draws one histogram per column
of a `cell_fractions` frame, with matplotlib directly rather than through
scanpy's plotting.

- **Give it the same rules as the filter.** `nmads`, `limits` and `band`
  are `cell_filters`' own arguments, and the panel calls it to get them, so
  the lines are the cuts the file was made with rather than a second
  opinion. Blue (`qc.KEPT`) passes every rule the metric carries, grey
  (`qc.REMOVED`) fails one.
- **Each line is labelled with the rule and the value it drew** -- `5 MAD
  13,063`, `limit  6` -- and each panel is titled with how many cells it
  put beyond them. A metric can carry several lines. A cut that
  removes a handful leaves no grey tall enough to see, so without the count
  a one-colour panel looks like a panel that failed to say anything.
- **Counts are two-sided and go on a log axis**, fractions are upper-tail
  only and linear (`LOG_METRICS`, overridable with `sides`).
- **The grey bars do not sum to the cells a filter drops.** Each panel
  counts its own rejects, so a cell that two metrics reject is grey twice.
  `mad_outlier` per metric, combined with `any`, is what gives the cells.
- A metric no rule covers loses its lines, not the figure: the panel is
  drawn without them and the caption says so.
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
- **A zero cannot sit on a log axis**, and matplotlib drops one without
  saying so. Such a point stays in the data file with `drawn` false, still
  counts in its quadrant, and the caption says how many there were.
- `unit="genes"` makes the same figure a per-gene comparison against a
  reference, thresholds and quadrant counts included.
- **`groups=` colours by a label instead of by the thresholds.** Where
  nothing is being filtered, calling the points below a line "removed" is
  wrong: the lines then only mark where things were expected to fall, and
  the colour is free to carry what the point is. A point with no label is
  grey, see-through and drawn behind the rest; a labelled one is solid
  and thinly outlined, since the alpha that keeps a cloud of ten thousand
  readable makes a handful of points invisible, and Okabe-Ito's yellow is
  a pale fill on a white page.

## Does the panel dominate?

`qc.panel_expression(data)` boxplots the mean expression per gene, in UMIs
per million, for the genes `var["readout"]` marks against every other gene.

- **Gene by gene, not cell by cell.** `pct_counts_targeted` says how much of
  a cell is the panel; this says whether that comes from the panel as a
  whole or from a few genes carrying it.
- `qc.mean_upm(counts)` is the number behind it: each cell scaled to a
  million UMIs, then averaged across cells, computed on the sparse matrix.
  `normalization.counts_per_million` is the same arithmetic for the dense
  guide-by-sample table of a bulk screen and would densify this one.
- One is added before the log, so a gene seen in no cell still appears.

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

## The rules, and the file they make

`qc.cell_filters(frame, nmads=5, limits=...)` turns the metric frame into one
boolean column per rule, and a table of the rules themselves.
`qc.apply_filters(data, passes, rules, path=...)` keeps the cells that pass
every one, records them, and writes the object.

- **Two kinds of rule, kept apart.** A band is measured from this channel's
  own spread; a limit is a number someone decided, `{metric: (low, high)}`
  with `None` for an open side. `total_counts_mad` and `total_counts_limit`
  are separate columns, so a threshold the screen settled on elsewhere sits
  next to the band measured here, and the record says which one cut what.
- **The panel draws these rules.** `metric_distributions` calls
  `cell_filters` for them, so pass both the same arguments and the figure
  cannot disagree with the file.
- **Measure the rules on the object being filtered.** `apply_filters`
  refuses a frame that does not cover exactly its cells: a cell no rule saw
  would otherwise pass every one of them.
- **A `NaN` fails its rule.** The modality did not measure that cell, and a
  metric a cell was never measured for is not a metric it passed.
- **`uns["qc"]` is the record** the file carries: `floor`, written by
  `filter_cells` itself; `rules`, every threshold with its source and how
  many cells it rejected; and the counts in and out. An open bound is an
  infinity, since `None` does not survive HDF5.
- **`band=` names the metrics worth a band**, every column by default. Not
  every metric deserves one: a ribosomal fraction with a long tail loses
  real cells to a band nobody chose.
- `nmads=None` leaves only the limits, for a screen whose thresholds were
  settled elsewhere.

## Ask before running

- **Which object?** One channel's MuData, or several concatenated. A knee
  plot over pooled channels hides a channel that failed.
- **Is this a targeted assay, and is the panel flagged?** `cell_fractions`
  adds `pct_counts_targeted` only when `var` carries the `readout` column,
  and `panel_expression` needs it too. Without the flag a TAP channel is
  QCed as if it were whole-transcriptome, and a `pct_counts_targeted` limit
  is a rule that cannot be written.
- **What is the floor?** `min_umis` has no default; read it off this
  channel's knee.
- **Where do the figures and the filtered object go?** Nothing is written
  outside the stem's folder, and `apply_filters` writes one `.h5mu`.
