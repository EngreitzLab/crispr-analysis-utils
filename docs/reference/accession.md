# crispr_analysis_utils.accession

A guide's accession encodes its **protospacer**, so the same target gets the
same id in every library and table, and two different protospacers can never
share one. `spacer_of` decodes it back.

The protospacer is the key, not the ordered sequence. Two guides differing
only 5' of the protospacer — one lab prepends a G for the U6 promoter, another
does not — are two oligos aimed at one site, so they share an accession and
differ in `guide_alignment`'s `guide_id`. The accession answers "same
target?"; `guide_id` answers "same oligo?".

Nothing in a sequence says where the protospacer ends, so a sequence with its
PAM still attached is accepted and gets a different accession, silently. One
spacer in sixteen ends in `GG`. Decide which bases are the protospacer first;
`guide_genome_annotation.tsv`'s `spacer` column already has.

## The encoding

For an uppercase protospacer:

1. Start from 1, and for each base multiply by 4 and add 0 for A, 1 for C,
   2 for G or 3 for T. The leading 1 keeps `AAA` apart from `AA`.
2. Write the result in Crockford's base 32, alphabet
   `0123456789ABCDEFGHJKMNPQRSTVWXYZ`, after `GU`.

Decoding ignores case and reads I and L as 1 and O as 0. A 20-nt protospacer
gives 9 characters after `GU`; a 19-nt one gives 8.

`ACGT` packs as `0b1_00_01_10_11`, which is 283, which is 8, 27 in base 32:
the accession is `GU8V`.

::: crispr_analysis_utils.accession
    options:
        members_order: source
