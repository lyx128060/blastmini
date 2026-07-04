import argparse
import sys
import csv
from .io import read_fasta
from .core import run_search

def main():
    parser = argparse.ArgumentParser(description="blastmini: Mini BLAST tool")
    subparsers = parser.add_subparsers(dest="command")

    search_parser = subparsers.add_parser("search")
    search_parser.add_argument("database", help="DB FASTA")
    search_parser.add_argument("query", help="Query FASTA")
    search_parser.add_argument("-k", "--kmer", type=int, default=8)
    search_parser.add_argument("-t", "--top", type=int, default=5)
    search_parser.add_argument("-o", "--output", default="results.tsv")

    args = parser.parse_args()

    if args.command == "search":
        print(f"[*] 正在加载数据库: {args.database}")
        try:
            db_records = list(read_fasta(args.database))
            query_records = list(read_fasta(args.query))
            
            print(f"[*] 正在进行全序列比对 (K-mer={args.kmer}, Top={args.top})...")
            # 这里的 run_search 把所有核心算法都跑完了
            results = run_search(query_records, db_records, k=args.kmer, top_n=args.top)
            
            # 将比对结果写入 TSV 文件
            fields = ["query_id", "subject_id", "score", "identity", "query_start", "query_end", "subject_start", "subject_end", "pvalue"]
            with open(args.output, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(results)
            
            print(f"[+] 搜索完毕！共找到 {len(results)} 条比对记录，结果已保存至: {args.output}")
            
        except Exception as e:
            print(f"[!] 发生错误: {e}")
            sys.exit(1)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()