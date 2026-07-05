import random
from typing import Optional


def calc_identity(matches: int, align_len: int) -> float:
    """Calculate sequence identity as matches / alignment length."""
    if align_len <= 0:
        return 0.0
    return matches / align_len


def best_ungapped_score(
    query_seq: str,
    subject_seq: str,
    match_score: int = 2,
    mismatch_score: int = -1,
) -> int:
    """Return the best local ungapped alignment score over all diagonals."""
    if not query_seq or not subject_seq:
        return 0

    best_score = 0
    q_len = len(query_seq)
    s_len = len(subject_seq)

    for offset in range(-q_len + 1, s_len):
        q_start = max(0, -offset)
        s_start = max(0, offset)
        diag_len = min(q_len - q_start, s_len - s_start)
        current = 0

        for delta in range(diag_len):
            q_char = query_seq[q_start + delta]
            s_char = subject_seq[s_start + delta]
            current += match_score if q_char == s_char else mismatch_score
            if current < 0:
                current = 0
            elif current > best_score:
                best_score = current

    return best_score


def calc_empirical_pvalue(
    real_score: int,
    query_seq: str,
    subject_seq: Optional[str] = None,
    iterations: int = 50,
    rng: Optional[random.Random] = None,
) -> float:
    """
    Estimate an empirical p-value by shuffling the query and rescoring it.

    The same subject sequence is used as the background target, and each random
    score is the best local ungapped score across all diagonals.
    """
    if real_score <= 0 or iterations <= 0:
        return 1.0

    target_seq = subject_seq if subject_seq is not None else query_seq
    if not query_seq or not target_seq:
        return 1.0

    random_source = rng or random
    query_chars = list(query_seq)
    count = 0

    for _ in range(iterations):
        shuffled = query_chars[:]
        random_source.shuffle(shuffled)
        random_score = best_ungapped_score("".join(shuffled), target_seq)
        if random_score >= real_score:
            count += 1

    return (count + 1) / (iterations + 1)
