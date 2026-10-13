"""The orientation verdicts, and what evidence stands behind one.

A guide table can be written 3' to 5' by mistake, in whole or in part. The
scan answers it by holding both readings of every guide: an entry is
`REVERSE_COMPLEMENT` only when the reversed reading has a site with a PAM and
the sequence as given does not.
"""

from __future__ import annotations

FORWARD = "forward"
"""The sequence as given has a site with a PAM 3' of it."""
REVERSE_COMPLEMENT = "reverse_complement"
"""Only the reversed reading does: this entry is written 3' to 5'."""
NO_PAM = "no_pam"
"""The spacer matches the genome, and no matching site has a PAM either way."""
NO_MATCH = "no_match"
"""The spacer matches nowhere. A non-targeting control belongs here."""
NOT_CHECKED = "not_checked"
"""The orientation check did not run."""
VERDICTS = (FORWARD, REVERSE_COMPLEMENT, NO_PAM, NO_MATCH, NOT_CHECKED)
"""Every orientation verdict, in report order."""

DISTINCT = "distinct"
"""The reversed reading was a different sequence, so the test was independent."""
AMBIGUOUS = "ambiguous"
"""The reversed reading was only the spacer's reverse complement.

It matches the same locus, so ``reverse_complement`` cannot be told from a
``CCN`` that fell 5' of the match by chance, which happens one time in
sixteen. Only a guide carrying no base 5' of its spacer is in this position:
with one, the trim moves to the other end and the reversed spacer is a
different sequence entirely.
"""
EVIDENCE = (DISTINCT, AMBIGUOUS)
"""What `orientation_evidence` can say, for a ``reverse_complement`` entry."""
