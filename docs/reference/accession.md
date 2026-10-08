# crispr_analysis_utils.accession

A guide's accession comes from its sequence alone, so a guide keeps it across
libraries, labs and renamings, and two names for one sequence share it. It is
computed from the sequence that is aligned: 5' to 3', with the 5' G when the
guide carries one, and without the PAM.

Outside Python, the accession of an uppercase sequence `s` is `GU` followed by
the first 12 hexadecimal digits, in uppercase, of the UUID version 5 of `s` in
the namespace `ee28d900-f51b-5e2b-9ca7-4c5808da45cc`.

Twelve hexadecimal digits are 48 bits. The chance that two of `n` different
sequences share an accession is about `n(n - 1) / 2^49`: 1 in 290,000 for a
library of 43,778 guides, 1 in 56,000 for 100,000 guides, and 1 in 560 for a
million. `guide_accessions` raises an error rather than return a collision, so
pass every library whose accessions must stay distinct in one call.

::: crispr_analysis_utils.accession
    options:
        members_order: source
