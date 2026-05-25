# Guide QC

Example usage:

```python
import pandas as pd
import crispr_analysis_utils as cau

cau.guide_qc.guides_to_fastq(
  "guides.tsv",
  "guides.fastq",
  pam="NGG",
  add_leading_g=True
)

# DataFrame input with column overrides
guides = pd.DataFrame({"id": ["g1"], "seq": ["ACGT"]})
cau.guide_qc.guides_to_fastq(
  guides,
  "guides.fastq",
  id_col="id",
  sequence_col="seq"
)
```

::: crispr_analysis_utils.guide_qc
