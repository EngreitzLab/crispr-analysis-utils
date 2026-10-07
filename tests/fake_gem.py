"""A stand-in for gem-indexer and gem-mapper, for tests that do not need GEM3.

`install` writes ``gem-indexer`` and ``gem-mapper`` wrappers into a folder; each
runs this file with the Python running the tests. The fake indexer writes the
reference's path into the index file. The fake mapper reads the reference back,
finds every placement of each read on both strands with at most
``--complete-search-error`` errors (a read N is always one error), and writes
them as SAM to standard output, the way gem-mapper lays its records out.

Environment variables change what the fakes do:

- ``FAKE_GEM_CALLS``: append each command line, as JSON, to this file.
- ``FAKE_GEM_CRASH``: the mapper writes the SAM header, then a crash report
  naming its first read, and exits with status 1.
- ``FAKE_GEM_SIGNAL_EXIT``: with ``FAKE_GEM_CRASH``, the exit status to use.
- ``FAKE_GEM_OUT_OF_RANGE``: the mapper adds two records per mapped read, one
  past the end of its contig and one at an underflowed position.
- ``FAKE_GEM_INDEX_FAIL``: the indexer writes part of the index and exits 1.
"""

import json
import os
import sys
from pathlib import Path

VERSION = "v0.0.0-fake"
UNDERFLOWED_POSITION = 2**64 - 22
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def install(folder):
    """Write gem-indexer and gem-mapper wrappers into `folder`; return it."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).resolve()
    for tool in ("gem-indexer", "gem-mapper"):
        wrapper = folder / tool
        wrapper.write_text(
            f'#!/bin/sh\nexec "{sys.executable}" "{script}" {tool} "$@"\n',
            encoding="utf-8",
        )
        wrapper.chmod(0o755)
    return folder


def calls(path):
    """The command lines the fakes recorded in `path`, oldest first."""
    if not Path(path).exists():
        return []
    return [json.loads(line) for line in Path(path).read_text().splitlines()]


def main(argv):
    tool, args = argv[0], argv[1:]
    if args == ["--version"]:
        print(VERSION)
        return 0
    record = os.environ.get("FAKE_GEM_CALLS")
    if record:
        with open(record, "a", encoding="utf-8") as handle:
            handle.write(json.dumps([tool, *args]) + "\n")
    if tool == "gem-indexer":
        return _index(args)
    return _map(args)


def _option(args, short, long=None):
    for position, arg in enumerate(args):
        if arg == short:
            return args[position + 1]
        if long and arg.startswith(long + "="):
            return arg.split("=", 1)[1]
    raise SystemExit(f"fake gem: missing {short}")


def _index(args):
    reference = _option(args, "-i")
    prefix = _option(args, "-o")
    scratch = _option(args, "--tmp-folder")
    if not scratch.endswith(os.sep) or not os.path.isdir(scratch):
        print(f"fake gem-indexer: bad --tmp-folder {scratch!r}", file=sys.stderr)
        return 1
    with open(prefix + ".gem", "w", encoding="utf-8") as handle:
        handle.write("" if os.environ.get("FAKE_GEM_INDEX_FAIL") else reference)
    if os.environ.get("FAKE_GEM_INDEX_FAIL"):
        print("GEM::FatalError (fake)", file=sys.stderr)
        return 1
    Path(prefix + ".info").write_text("fake index\n", encoding="utf-8")
    print(f"[GEM Index '{prefix}.gem' was successfully built]", file=sys.stderr)
    return 0


def _read_fasta(path):
    contigs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            contigs[name] = []
        elif name is not None:
            contigs[name].append(line.strip())
    return {name: "".join(parts).upper() for name, parts in contigs.items()}


def _read_fastq(path):
    lines = Path(path).read_text().splitlines()
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 4)]


def _errors(read, window, budget):
    errors = 0
    for read_base, genome_base in zip(read, window, strict=True):
        if read_base == "N" or read_base != genome_base:
            errors += 1
            if errors > budget:
                break
    return errors


def _map(args):
    index = _option(args, "-I")
    fastq = _option(args, "-i")
    budget = int(_option(args, "-E", "--complete-search-error"))
    contigs = _read_fasta(Path(index).read_text().strip())
    out = sys.stdout
    out.write("@HD\tVN:1.4\tSO:unsorted\n")
    for name, sequence in contigs.items():
        out.write(f"@SQ\tSN:{name}\tLN:{len(sequence)}\n")
    out.write(f"@PG\tID:GEM\tPN:gem-mapper\tVN:{VERSION}\n")
    reads = _read_fastq(fastq)
    if os.environ.get("FAKE_GEM_CRASH"):
        out.flush()
        name, sequence = reads[0]
        sys.stderr.write(
            "GEM::Unexpected error occurred. Sorry for the inconvenience\n"
            "GEM::Input.State\n"
            f"Sequence (File '{fastq}' Line '1')\n@{name} \n{sequence}\n+\n"
            ">> GEM.System.Error::Signal raised (no=11)\n"
        )
        return int(os.environ.get("FAKE_GEM_SIGNAL_EXIT", "1"))
    first_contig = next(iter(contigs))
    for name, sequence in reads:
        quality = "I" * len(sequence)
        hits = []
        for contig, genome in contigs.items():
            for strand, read in (("+", sequence), ("-", _revcomp(sequence))):
                for start in range(len(genome) - len(read) + 1):
                    window = genome[start : start + len(read)]
                    errors = _errors(read, window, budget)
                    if errors <= budget:
                        hits.append((errors, contig, start, strand, read))
        hits.sort()
        if not hits:
            out.write(f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}\n")
            continue
        for rank, (errors, contig, start, strand, read) in enumerate(hits):
            flag = (16 if strand == "-" else 0) | (256 if rank else 0)
            seq, qual = (read, quality) if not rank else ("*", "*")
            out.write(
                f"{name}\t{flag}\t{contig}\t{start + 1}\t{0 if rank else 60}\t"
                f"{len(read)}M\t*\t0\t0\t{seq}\t{qual}\tNM:i:{errors}\n"
            )
        if os.environ.get("FAKE_GEM_OUT_OF_RANGE"):
            length = len(contigs[first_contig])
            for position in (length + 1, UNDERFLOWED_POSITION):
                out.write(
                    f"{name}\t272\t{first_contig}\t{position}\t0\t"
                    f"{len(sequence)}M\t*\t0\t0\t*\t*\tNM:i:1\n"
                )
    return 0


def _revcomp(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
