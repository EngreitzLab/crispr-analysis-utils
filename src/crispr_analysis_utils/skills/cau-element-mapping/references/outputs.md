# The output of `cau element-mapping`

<!-- --8<-- [start:tables] -->
`--outdir` gets three folders and nothing loose:

```text
elements/        guide_element_mapping.tsv  columns.yaml  qc.tsv
run_metadata/    element_mapping_<UTC timestamp>.json
logs/            element_mapping.log
```

## `guide_element_mapping.tsv`

**One row per site per element it overlaps.** A site inside two overlapping
elements gets two rows. A guide that matches nowhere gets one row with no
site and no element.

| Column | Holds |
| --- | --- |
| `guide_id` | The aligned sequence. With the four site columns, the key back to `guides_alignment_metadata.tsv`. |
| `chr`, `start`, `end`, `strand` | The protospacer, 0-based and half-open. Empty when the guide matches nowhere. |
| `element_id` | The element this site overlaps. Empty unless targeting. |
| `element_category` | The elements BED's fifth column, when it has one. |
| `target_class` | `targeting`, `safe_targeting` or `non_targeting`. |
| `distance_to_tss` | Signed bases from the cut site to the nearest protein-coding TSS, carried from the annotation. |

### The classes are per site

A guide with two sites whose classes differ reports both rows, because the
two sites really do different things:

```text
guide_id   chr    start     element_id               target_class
G...NGG    chr10  86958386  chr10:86958280-86958780  targeting
G...NGG    chr17  79714569                           safe_targeting
```

### `safe_targeting` is not a promise

It says only that the site maps and overlaps none of the elements given. How
far from a gene a guide has to be before it is safe depends on what the
screen reads out, so that is left to `distance_to_tss` rather than fixed by a
threshold here. Threshold it yourself, and say which threshold you used.

It follows that an incomplete element list inflates `safe_targeting`: every
guide for a missing element lands there. Check the count against the library
design before believing it.

## `qc.tsv`

Three columns, `section`, `name`, `value`: how the library covers the
elements, which is the question the mapping table cannot answer by itself.

| Section | Says |
| --- | --- |
| `classes` | How many guides, and how many sites in each class |
| `elements` | Total, how many have a guide, how many have none, and how many are below half the median density |
| `guides_per_element` | Median, min, max, and the median per 100 bases |
| `category` | The same per element category, so a thinly covered kind of element stands out |
| `untiled_elements` | How far the elements with no guide are from the nearest mapped one |

Two of these earn their place:

- **Density, not count.** A wide element collects more guides without being
  better covered, so `median_per_100bp` and `below_half_median_density` are
  what find the elements that are actually thin.
- **Distance for the untiled ones.** An element a few hundred bases from the
  nearest guide was missed at its boundary; one tens of kilobases away was
  never targeted at all. Only the second is a hole in the design.

### `element_id`

Taken from the BED's fourth column when it has one, so an id used elsewhere
keeps its spelling. Otherwise built as `chr:start-end` with a 1-based start,
which is how a locus is read aloud and how most element tables write it.
<!-- --8<-- [end:tables] -->

## Checking a library's own labels

This is the step that tests a design's control sets, so expect disagreement
and report it rather than smoothing it over. Two kinds, which mean opposite
things:

- **a guide the design calls non-targeting that lands inside an element.**
  Read the genomic sequence back before believing either side; an aligner
  that is not exhaustive can miss a perfect site and leave a real cutter
  filed as a control.
- **a guide the design calls safe-targeting whose `distance_to_tss` is
  small.** Usually the design measured safety against its own readout genes
  rather than every promoter. Report the distance, not the label.

A third, quieter one: if the design stripped a 5' G before aligning, it
searched a 19-mer, and a 19-mer matches places the real guide does not. Those
guides look targeting in the old table and map nowhere here.
