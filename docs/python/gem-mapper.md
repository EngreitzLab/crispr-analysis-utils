# GEM Mapper

Install GEM3 with conda:

```bash
conda create -n gem3-map -c conda-forge -c bioconda --strict-channel-priority gem3-mapper -y
conda activate gem3-map
```

Download IGVF hg38 reference FASTA:

```bash
mkdir -p ../../annotations/ENCODE/hg38
curl -L "https://api.data.igvf.org/reference-files/IGVFFI0653VCGH/@@download/IGVFFI0653VCGH.fasta.gz" \
  -o ../../annotations/ENCODE/hg38/IGVFFI0653VCGH.fasta.gz
gunzip -c ../../annotations/ENCODE/hg38/IGVFFI0653VCGH.fasta.gz \
  > ../../annotations/ENCODE/hg38/IGVFFI0653VCGH.fasta
```

Reference accession: `IGVFFI0653VCGH`

Run from Python:

```python
import crispr_analysis_utils as cau

# Build index
cau.gem_mapper.build_gem_index(
  "../../annotations/ENCODE/hg38/IGVFFI0653VCGH.fasta",
  "../../annotations/ENCODE/hg38/gem_index/IGVFFI0653VCGH",
  threads=8
)

# Map guides
cau.gem_mapper.map_guides_with_gem(
  "../../annotations/ENCODE/hg38/gem_index/IGVFFI0653VCGH.gem",
  "guides_input.fastq",
  "guides_mapped_hg38.sam",
  mapping_mode="sensitive",
  threads=8,
  sam_compact=False
)
```

By default this writes mapper logs to `guides_mapped_hg38.log` in the same
directory as the SAM output. You can override with `log_path="path/to/log.log"`.

Equivalent command executed by `map_guides_with_gem`:

```bash
gem-mapper \
  -I ../../annotations/ENCODE/hg38/gem_index/IGVFFI0653VCGH.gem \
  -i guides_input.fastq \
  -o guides_mapped_hg38.sam \
  --sam-compact=false \
  --mapping-mode sensitive \
  --threads 8
```

::: crispr_analysis_utils.gem_mapper
