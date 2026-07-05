from pathlib import Path
from typing import Iterable, Iterator, Tuple


ALLOWED_SEQUENCE_CHARS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ*")


def clean_sequence(sequence: Iterable[str]) -> str:
    """Return an uppercase sequence with whitespace, gaps, and digits removed."""
    cleaned = []
    for char in "".join(sequence).upper():
        if char in ALLOWED_SEQUENCE_CHARS:
            cleaned.append(char)
    return "".join(cleaned)


def read_fasta(file_path: str) -> Iterator[Tuple[str, str]]:
    """
    Stream records from a FASTA file.

    Yields:
        (seq_id, sequence), where seq_id is the first token after ">" and
        sequence has been uppercased and cleaned.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Error: 找不到文件 {file_path}")

    with path.open("r", encoding="utf-8") as handle:
        seq_id = None
        seq_lines = []

        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue

            if line.startswith(">"):
                if seq_id is not None:
                    sequence = clean_sequence(seq_lines)
                    if not sequence:
                        raise ValueError(f"FASTA record '{seq_id}' has no sequence")
                    yield seq_id, sequence

                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_no}")

                seq_id = header.split()[0]
                seq_lines = []
                continue

            if seq_id is None:
                raise ValueError(f"Sequence data found before FASTA header at line {line_no}")

            seq_lines.append(line)

        if seq_id is not None:
            sequence = clean_sequence(seq_lines)
            if not sequence:
                raise ValueError(f"FASTA record '{seq_id}' has no sequence")
            yield seq_id, sequence
