import random
from collections import defaultdict
from typing import Dict, Iterable, List, MutableMapping, Optional, Tuple

from .stats import calc_empirical_pvalue, calc_identity


MATCH_SCORE = 2
MISMATCH_SCORE = -1
DEFAULT_X_DROP = 3


def _validate_positive_int(name: str, value: int) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer")


def build_kmer_index(db_records: Iterable[Tuple[str, str]], k: int = 8) -> Dict[str, List[Tuple[str, int]]]:
    """Build a subject k-mer index mapping k-mer -> [(subject_id, position)]."""
    _validate_positive_int("k", k)
    index: MutableMapping[str, List[Tuple[str, int]]] = defaultdict(list)

    for seq_id, seq in db_records:
        if len(seq) < k:
            continue
        for pos in range(len(seq) - k + 1):
            index[seq[pos : pos + k]].append((seq_id, pos))

    return dict(index)


def find_seeds(query_seq: str, kmer_index: Dict[str, List[Tuple[str, int]]], k: int = 8) -> List[Tuple[str, int, int]]:
    """Return seed hits as (subject_id, subject_pos, query_pos)."""
    _validate_positive_int("k", k)
    if len(query_seq) < k:
        return []

    seeds = []
    for query_pos in range(len(query_seq) - k + 1):
        kmer = query_seq[query_pos : query_pos + k]
        for subject_id, subject_pos in kmer_index.get(kmer, []):
            seeds.append((subject_id, subject_pos, query_pos))
    return seeds


def extend_direction(
    query_seq: str,
    subject_seq: str,
    q_pos: int,
    s_pos: int,
    step: int,
    x_drop: int = DEFAULT_X_DROP,
) -> Tuple[int, int, int, int]:
    """Extend one side of a seed with ungapped X-drop scoring."""
    if step not in (-1, 1):
        raise ValueError("step must be -1 or 1")
    if x_drop < 0:
        raise ValueError("x_drop must be >= 0")

    score = 0
    best_score = 0
    best_q = q_pos
    best_s = s_pos
    matches = 0
    best_matches = 0

    curr_q = q_pos + step
    curr_s = s_pos + step

    while 0 <= curr_q < len(query_seq) and 0 <= curr_s < len(subject_seq):
        if query_seq[curr_q] == subject_seq[curr_s]:
            score += MATCH_SCORE
            matches += 1
        else:
            score += MISMATCH_SCORE

        if score > best_score:
            best_score = score
            best_q = curr_q
            best_s = curr_s
            best_matches = matches

        if best_score - score > x_drop:
            break

        curr_q += step
        curr_s += step

    return best_score, best_q, best_s, best_matches


def _extend_seed(
    query_seq: str,
    subject_seq: str,
    q_pos: int,
    s_pos: int,
    k: int,
    x_drop: int,
) -> Tuple[int, int, int, int, int, int, int]:
    left_score, q_start, s_start, left_matches = extend_direction(
        query_seq, subject_seq, q_pos, s_pos, -1, x_drop=x_drop
    )
    right_score, q_end, s_end, right_matches = extend_direction(
        query_seq, subject_seq, q_pos + k - 1, s_pos + k - 1, 1, x_drop=x_drop
    )

    total_score = (k * MATCH_SCORE) + left_score + right_score
    total_matches = k + left_matches + right_matches
    align_len = q_end - q_start + 1
    return total_score, total_matches, q_start, q_end, s_start, s_end, align_len


def _hit_sort_key(hit: Dict[str, object]) -> Tuple[int, int, str, int]:
    return (
        -int(hit["score"]),
        -(int(hit["query_end"]) - int(hit["query_start"]) + 1),
        str(hit["subject_id"]),
        int(hit["subject_start"]),
    )


def _prune_hits(hit_map: Dict[Tuple[str, int, int, int, int], Dict[str, object]], top_n: int) -> None:
    max_cached_hits = max(top_n * 10, top_n + 20)
    if len(hit_map) <= max_cached_hits:
        return

    best_items = sorted(hit_map.items(), key=lambda item: _hit_sort_key(item[1]))[:top_n]
    hit_map.clear()
    hit_map.update(best_items)


def _build_query_index(
    query_records: Iterable[Tuple[str, str]],
    k: int,
) -> Tuple[List[Tuple[str, str]], Dict[str, List[Tuple[int, int]]]]:
    query_list = list(query_records)
    query_index: MutableMapping[str, List[Tuple[int, int]]] = defaultdict(list)

    for query_idx, (_query_id, query_seq) in enumerate(query_list):
        if len(query_seq) < k:
            continue
        for query_pos in range(len(query_seq) - k + 1):
            query_index[query_seq[query_pos : query_pos + k]].append((query_idx, query_pos))

    return query_list, dict(query_index)


def run_search(
    query_records: Iterable[Tuple[str, str]],
    db_records: Iterable[Tuple[str, str]],
    k: int = 8,
    top_n: int = 5,
    pvalue_iterations: int = 50,
    random_seed: Optional[int] = None,
    min_score: Optional[int] = None,
    x_drop: int = DEFAULT_X_DROP,
) -> List[Dict[str, object]]:
    """Run seed-and-extend search and return TSV-ready result dictionaries."""
    _validate_positive_int("k", k)
    _validate_positive_int("top_n", top_n)
    if pvalue_iterations < 0:
        raise ValueError("pvalue_iterations must be >= 0")
    if x_drop < 0:
        raise ValueError("x_drop must be >= 0")

    query_list, query_index = _build_query_index(query_records, k)
    if not query_list or not query_index:
        return []

    rng = random.Random(random_seed) if random_seed is not None else None
    top_hits_by_query: List[Dict[Tuple[str, int, int, int, int], Dict[str, object]]] = [
        {} for _ in query_list
    ]
    all_results: List[Dict[str, object]] = []

    for s_id, s_seq in db_records:
        if len(s_seq) < k:
            continue
        for s_pos in range(len(s_seq) - k + 1):
            kmer = s_seq[s_pos : s_pos + k]
            seed_hits = query_index.get(kmer)
            if not seed_hits:
                continue

            for query_idx, q_pos in seed_hits:
                q_id, q_seq = query_list[query_idx]
                hits_by_span = top_hits_by_query[query_idx]

                score, matches, q_start, q_end, s_start, s_end, align_len = _extend_seed(
                    q_seq, s_seq, q_pos, s_pos, k, x_drop
                )

                if min_score is not None and score < min_score:
                    continue

                span_key = (s_id, q_start, q_end, s_start, s_end)
                identity_val = calc_identity(matches, align_len)
                hit = {
                    "query_id": q_id,
                    "subject_id": s_id,
                    "score": score,
                    "identity": f"{identity_val:.2%}",
                    "query_start": q_start + 1,
                    "query_end": q_end + 1,
                    "subject_start": s_start + 1,
                    "subject_end": s_end + 1,
                    "pvalue": "1.0000",
                    "_subject_seq": s_seq,
                }

                previous = hits_by_span.get(span_key)
                if previous is None or score > previous["score"]:
                    hits_by_span[span_key] = hit
                    _prune_hits(hits_by_span, top_n)

    for query_idx, (_q_id, q_seq) in enumerate(query_list):
        hits = sorted(top_hits_by_query[query_idx].values(), key=_hit_sort_key)[:top_n]

        for hit in hits:
            subject_seq = str(hit.pop("_subject_seq"))
            p_val = calc_empirical_pvalue(
                int(hit["score"]),
                q_seq,
                subject_seq,
                iterations=pvalue_iterations,
                rng=rng,
            )
            hit["pvalue"] = f"{p_val:.4f}"

        all_results.extend(hits)

    return all_results
