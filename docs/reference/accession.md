# crispr_analysis_utils.accession

A guide's accession encodes the guide itself, so it is the same in every
library and table, and two guides can never share one. `guide_sequence`
decodes it back.

The guide is written `G...NGG`, 5' to 3'. `guide_accession` refuses anything
else rather than guess where the PAM ends or whether a G was added: both
answers belong to whoever designed the library. The PAM is always spelled
`NGG`, whatever the genome holds at the N, so a guide keeps one accession
across its genomic sites.

## The encoding

For an uppercase guide:

1. Drop the three PAM characters.
2. Start from 1, and for each remaining base multiply by 4 and add 0 for A,
   1 for C, 2 for G or 3 for T. The leading 1 keeps `GAA` apart from `GA`.
3. Write the result in Crockford's base 32, alphabet
   `0123456789ABCDEFGHJKMNPQRSTVWXYZ`, after `GU`.

Decoding ignores case and reads I and L as 1 and O as 0. A guide with 20 to 22
bases before the PAM gets 9 characters after `GU`; one with 19 gets 8.

`GACGTNGG` drops `NGG` and packs `GACGT` as `0b1_10_00_01_10_11`, which is
1563, which is 1, 16, 27 in base 32: the accession is `GU1GV`.

::: crispr_analysis_utils.accession
    options:
        members_order: source
