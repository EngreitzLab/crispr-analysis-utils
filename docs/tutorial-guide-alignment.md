# Tutorial: guide alignment on a cluster

This tutorial aligns a guide library to hg38 on a SLURM cluster, from a fresh
clone to the summary of the run. The [guide](guide-alignment.md) explains
each step and every output.

## 1. Install

On the cluster's login node, clone the repository and install its default
pixi environment, which holds the package, pysam and GEM3:

```bash
git clone https://github.com/EngreitzLab/crispr-analysis-utils.git ~/crispr-analysis-utils
cd ~/crispr-analysis-utils
pixi install
pixi run gem-mapper --version
```

Install first, on the login node: `pixi run` installs a missing environment
before it runs anything, which a compute node may not have the network for.
The jobs below run `pixi run -m ~/crispr-analysis-utils`, which uses that
environment from any folder.

## 2. Download the reference

The IGVF hg38 reference, accession `IGVFFI0653VCGH`. The alignment needs an
uncompressed FASTA, so decompress it:

```bash
mkdir -p ~/guide-alignment/reference
cd ~/guide-alignment/reference
curl -L "https://api.data.igvf.org/reference-files/IGVFFI0653VCGH/@@download/IGVFFI0653VCGH.fasta.gz" \
    -o IGVFFI0653VCGH.fasta.gz
gunzip IGVFFI0653VCGH.fasta.gz
grep '^>' IGVFFI0653VCGH.fasta | head -3  # the first contig names
```

GEM 3.6 can crash on a read that aligns across the first base of the
reference's first contig (see
[Known limitations](guide-alignment.md#known-limitations)). A FASTA whose
first contig starts with N, as hg38's `chr1` does, cannot trigger it: the
`grep` above shows the first contig, and `sed -n 2p IGVFFI0653VCGH.fasta`
its first bases.

## 3. Prepare the guide table

One guide per line, tab-separated, without a header: an id, then the
sequence, 5' to 3' and without the PAM. For example, with two made-up
guides:

```bash
cd ~/guide-alignment
printf 'guide_0001\tGACCTCAGCAGGCTCTTGCA\nguide_0002\tTTGACCAGGCACACAGTTCC\n' > guides.tsv
```

For another layout, such as a CSV with column names, add
`--sep comma --header` to the command below and name the columns with
`--id-col` and `--spacer-col`.

## 4. Build the index

Indexing the genome is the slow step, done once per reference. Save this as
`~/guide-alignment/index.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=gem-index
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out

set -euo pipefail
cd ~/guide-alignment
pixi run -m ~/crispr-analysis-utils cau guide-alignment index \
    reference/IGVFFI0653VCGH.fasta \
    reference/gem_index/IGVFFI0653VCGH \
    --threads "$SLURM_CPUS_PER_TASK"
```

Adjust the resources to your cluster (and add a `--partition` line if it
needs one), then submit it with `sbatch index.sbatch`. It writes `reference/gem_index/IGVFFI0653VCGH.gem`,
gem-indexer's `.info` and `.log`, and `IGVFFI0653VCGH.gem.json`, which
records the reference's MD5 checksum.

## 5. Align the library

Save this as `~/guide-alignment/align.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=guide-alignment
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out

set -euo pipefail
cd ~/guide-alignment
pixi run -m ~/crispr-analysis-utils cau guide-alignment run \
    --guides guides.tsv \
    --reference reference/IGVFFI0653VCGH.fasta \
    --index reference/gem_index/IGVFFI0653VCGH.gem \
    --outdir results \
    --add-leading-g \
    --threads "$SLURM_CPUS_PER_TASK"
```

Keep `--add-leading-g` only if the guides are expressed with an added 5' G,
as from a U6 promoter; the [options](guide-alignment.md#options) list the
rest. Submit it with `sbatch align.sbatch`, then follow the run:

```bash
squeue -u "$USER"
tail -f results/logs/guide_alignment.log
```

To count only the primary chromosomes, leaving out unplaced contigs and
other sequences, add `--contigs "$contigs"` to the command, with:

```bash
contigs=$(printf 'chr%s,' {1..22} X Y M)
contigs=${contigs%,}  # chr1,chr2,...,chrY,chrM
```

Every name must be in the FASTA, or the run stops with an error naming the
missing ones.

## 6. Read the results

The job's output ends with the class counts, and `results/summary.tsv` holds
them all:

```bash
column -t results/summary.tsv | head -15
```

- **`class`**: how many guides are `unique`, `multi`, `off_target_only` and
  `no_site`. Non-targeting controls are expected to be `off_target_only` at
  the default 3 mismatches: a non-targeting guide that comes out `unique` or
  `multi` matches the genome perfectly.
- **`orientation`**: mostly `ok` for a library given 5' to 3'. Many
  `reversed` guides mean the table holds reverse complements; the run warns
  when they outnumber `ok` ones, and `orientation.tsv` shows each hit's
  flanks.
- **`site_search`**: the hits kept (`n_sites`) and those dropped, by reason.

`guides.tsv` has one row per guide, with the coordinates and cut site of the
unique ones; `guides.bed` and `cut_sites.bed` hold the same as BED, ready for
`bedtools`; `sites.tsv` lists every site kept, with up to 3 mismatches.

## 7. Rerun

Running the same command into the same `--outdir` again reuses the GEM
passes whose reads, options, index and GEM version are unchanged, as
recorded in `results/run.json`, and redoes the rest. Changing
`--max-mismatches`, the PAMs or `--add-leading-g` maps the passes they
affect again; changing only `--contigs` or `--threads` maps nothing again. A run that
stopped can be restarted the same way. `results/logs/guide_alignment.log`
keeps the log of every run.
