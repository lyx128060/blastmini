import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable, Optional

from . import __version__
from .core import DEFAULT_X_DROP, run_search
from .io import read_fasta


RESULT_FIELDS = [
    "query_id",
    "subject_id",
    "score",
    "identity",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "pvalue",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="blastmini: Mini BLAST-style sequence search")
    parser.add_argument("--version", action="version", version=f"blastmini {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)
    search_parser = subparsers.add_parser("search", help="search query FASTA records against a FASTA database")
    search_parser.add_argument("database", help="database FASTA")
    search_parser.add_argument("query", help="query FASTA")
    search_parser.add_argument("-k", "--kmer", type=int, default=8, help="k-mer seed length")
    search_parser.add_argument("-t", "--top", type=int, default=5, help="top hits to keep for each query")
    search_parser.add_argument("-o", "--output", default="results.tsv", help="output TSV path")
    search_parser.add_argument("--min-score", type=int, default=None, help="discard hits below this score")
    search_parser.add_argument("--x-drop", type=int, default=DEFAULT_X_DROP, help="ungapped extension X-drop cutoff")
    search_parser.add_argument(
        "--pvalue-iterations",
        type=int,
        default=50,
        help="query-shuffle iterations for empirical p-value; use 0 to disable",
    )
    search_parser.add_argument("--seed", type=int, default=None, help="random seed for reproducible p-values")
    return parser


def write_results(results: Iterable[dict], output_path: str) -> None:
    path = Path(output_path)
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "search":
        try:
            print(f"[*] 正在加载数据库: {args.database}", file=sys.stderr)
            db_records = read_fasta(args.database)
            query_records = list(read_fasta(args.query))

            print(
                f"[*] 正在进行全序列比对 (K-mer={args.kmer}, Top={args.top})...",
                file=sys.stderr,
            )
            results = run_search(
                query_records,
                db_records,
                k=args.kmer,
                top_n=args.top,
                pvalue_iterations=args.pvalue_iterations,
                random_seed=args.seed,
                min_score=args.min_score,
                x_drop=args.x_drop,
            )
            write_results(results, args.output)
            print(f"[+] 搜索完毕！共找到 {len(results)} 条比对记录，结果已保存至: {args.output}", file=sys.stderr)
            return 0
        except Exception as exc:
            print(f"[!] 发生错误: {exc}", file=sys.stderr)
            return 1

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
