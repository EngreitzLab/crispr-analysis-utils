# Tutorial: Guide Alignment QC

This tutorial walks from environment setup to final QC outputs for guide
alignment against hg38.

## 1. Create and activate environment

```bash
conda create -n gem3-map -c conda-forge -c bioconda --strict-channel-priority \
  python=3.11 gem3-mapper pysam pandas -y
conda activate gem3-map
python -m pip install -e .
```

## 2. Download IGVF hg38 reference

Accession: `IGVFFI0653VCGH`

```bash
mkdir -p data
curl -L "https://api.data.igvf.org/reference-files/IGVFFI0653VCGH/@@download/IGVFFI0653VCGH.fasta.gz" \
  -o data/IGVFFI0653VCGH.fasta.gz
gunzip -c data/IGVFFI0653VCGH.fasta.gz > data/IGVFFI0653VCGH.fasta
```

## 3. Prepare input files

- `data/guides.tsv`: two columns (guide id, guide sequence), tab-separated
- `data/hg38.chrom.sizes`: tab-separated chrom sizes used as allowed contig list

## 4. Run the pipeline

Use the single entry script:

```bash
python scripts/run_guide_alignment_qc.py \
  --guides-tsv data/guides.tsv \
  --reference-fasta data/IGVFFI0653VCGH.fasta \
  --chromsizes data/hg38.chrom.sizes \
  --outdir results/guide_alignment_qc \
  --threads 8 \
  --pam NGG \
  --add-leading-g
```

## 5. Final outputs

In `results/guide_alignment_qc/`:

- `logs/pipeline.log`

In `results/guide_alignment_qc/gem_index/`:

- `genome_index.gem` (+ index sidecar files)
- `gem_index.log`
- `index_command.sh`
- `index_inputs.txt`

In `results/guide_alignment_qc/alignment_outputs/`:

- `guides_input.fastq`
- `guides_mapped.sam`
- `guides_mapped.log`

In `results/guide_alignment_qc/guide_alignments_outputs/`:

- `valid_alignments.bed` (protospacer coordinates; excludes PAM)
- `discarded_alignments.tsv`
- `unmapped.tsv`
- `guide_alignment_log.tsv`
- `invalid_alignments.tsv`
- `alignment_summary.tsv`

## 6. Run on SLURM (one Python script)

Use the provided template:

```bash
sbatch scripts/run_guide_alignment_qc.sbatch
```

Template includes:

- conda activation
- cluster resources
- one `python` command for the full workflow
