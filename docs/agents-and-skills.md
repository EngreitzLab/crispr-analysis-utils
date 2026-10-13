# Agents and skills

--8<-- "README.md:skills"

## What ships

| Skill | Covers |
| --- | --- |
| `cau-accession` | `cau.accession`: accession ids for SpCas9 guides, encoded from the guide's own sequence |
| `cau-normalization` | `cau.counts_per_million` |
| `cau-guide-alignment` | `cau.guide_alignment` and the `cau guide-alignment` command: finding every guide's genomic site, and reading the outputs |
| `cau-igvf-export` | `cau.igvf_export` and the `cau igvf-export` command: writing a guide design in the IGVF per-guide metadata format |
| `cau-element-mapping` | `cau.element_mapping` and the `cau element-mapping` command: mapping guides to a screen's elements, and checking a design's control sets |
| `cau-guide-annotation` | `cau.guide_annotation` and the `cau guide-annotation` command: building the GENCODE store, placing sites in genes, and reading the output |

`cau install-skills --list` prints the same list from the installed package.

## Skill or agent

Every analysis module has a skill. A skill loads into your conversation, so
you see and steer each step: it suits a module whose use is a matter of
knowing the functions and choosing the parameters.

An agent is added when using a module means a long multi-step job (run a
pipeline, read its logs, summarize the outcome) whose intermediate output
would flood the conversation. The agent works in its own context and returns
a summary. It preloads the module's skill rather than repeating it, so the
module still has exactly one description of how it works.

## Where they live

The skills ship inside the Python package, under
`crispr_analysis_utils/skills/<name>/SKILL.md`, and agents under
`crispr_analysis_utils/agents/`. `cau install-skills` copies them from the
installed package. The plugin serves the same directory straight from the
repository's `main` branch.
