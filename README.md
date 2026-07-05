# BlastMini: A Lightweight Sequence Search Tool

**BlastMini** 是一个基于纯 Python 开发的轻量级生物序列相似性搜索工具包。本项目作为生物信息学核心算法的教学与工程实践，完整复现了 NCBI BLAST 的核心“种子-延伸（Seed-and-Extend）”启发式比对框架，并提供了标准的命令行接口（CLI）。

---

## 核心特性 (Key Features)

* **高效种子检索**：使用 $k$-mer 字典进行 $O(1)$ 复杂度的种子查找，并对候选命中执行局部延伸。
* **内存流式读取**：底层基于 Python 生成器（Generator）处理 FASTA 文件，支持在计算节点上顺畅解析百兆级（如 UniProt/RefSeq）数据库，有效防止内存溢出（OOM）。
* **双向无间隙延伸 (Ungapped Extension)**：实现精确的种子命中后，执行动态局部打分（Match +2, Mismatch -1），并在分数跌破阈值时自动截断（X-drop 机制）。
* **统计学显著性评估**：内置蒙特卡洛序列打乱机制（Shuffle），动态模拟随机背景分布，并输出严谨的经验 P-value (Empirical P-value)。
* **零外部依赖**：无需安装任何第三方依赖包，原生 Python 环境即可独立运行。

---

## 快速开始 (Quick Start)

### 1. 获取代码
```bash
git clone https://github.com/lyx128060/blastmini.git
cd blastmini
```

### 2. 配置运行环境

本项目没有第三方运行依赖。未安装包时，可以直接把 `src` 加入 `PYTHONPATH`：

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"
```

如果本机有 `pip`，也可以使用可编辑安装：

```bash
python3 -m pip install -e .
```

### 3. 运行序列搜索

工具包提供了标准的命令行调用方式。输入包含目标数据库文件（FASTA）与查询序列文件（FASTA）。

```bash
python3 -m blastmini.cli search <database.fasta> <query.fasta> -k 3 --top 5 -o output/results.tsv
```

例如，使用项目自带的 UniProt 数据库搜索胰高血糖素测试序列：

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"

python3 -m blastmini.cli search \
  data/uniprot_sprot.fasta \
  data/测试序列/seq_1_glucagon.fasta \
  -k 3 \
  --top 5 \
  --seed 1 \
  -o output/seq_1_glucagon_results.tsv
```

运行完成后，结果文件会保存在 `output/` 目录中。

**核心参数说明：**

* `-k`: K-mer 长度（推荐：核酸序列设为 4-6，蛋白质序列设为 3-4）。
* `--top`: 每个 Query 期望输出的最高分匹配记录数。
* `-o`: 输出结果的 TSV 文件路径。
* `--pvalue-iterations`: 经验 P-value 的 Query shuffle 次数；设为 `0` 可关闭统计模拟以便快速调试。
* `--seed`: 固定随机种子，用于复现实验中的 P-value。
* `--x-drop`: 无间隙延伸的 X-drop 截断阈值。

---

## HPC 集群批处理 (Slurm Integration)

针对高性能计算集群环境，本项目支持通过提交脚本进行批量比对。建议在 CPU 计算节点（如 `dcu_htc` 或 `debug` 队列）运行。

示例批处理逻辑：

```bash
# 激活环境并定义路径
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"

# 遍历目录执行比对
for query_file in /path/to/query_dir/*.fasta; do
    filename=$(basename -- "$query_file")
    prefix="${filename%.*}"
    
    python3 -m blastmini.cli search /path/to/db.fasta "$query_file" -k 3 --top 5 -o "output/${prefix}_results.tsv"
done

```

---

## 开发与测试 (Development)

项目测试使用 Python 标准库 `unittest`，无需额外安装测试框架：

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"
python3 -m unittest discover -s tests
python3 -m compileall src tests
```

---

## 工程模块架构 (Architecture)

本项目遵循模块化工程规范，核心代码结构如下：

* `io.py`：负责大文件流式读取与 FASTA 序列清洗（格式校验、去除非法字符）。
* `core.py`：核心检索引擎，包含 $k$-mer 索引构建与双向延伸打分逻辑。
* `stats.py`：统计评估模块，负责计算序列一致性（Identity）与经验 P-value。
* `cli.py`：用户交互接口，负责命令行参数解析与任务流调度。

---

## 示例输出 (Example Output)

运行结束后，系统生成标准的制表符分隔（TSV）结果文件，包含详细的比对坐标与统计学指标：

| query_id | subject_id | score | identity | query_start | query_end | subject_start | subject_end | pvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| query_1 | sp-P0A7V0-RL10_ECOLI | 87 | 95.12% | 1 | 41 | 22 | 62 | 0.0196 |
| query_2 | sp-P01308-INS_HUMAN | 62 | 100.00% | 1 | 31 | 80 | 110 | 0.0196 |

---

## 开发者 (Contributors)

本项目由以下小组成员共同完成，严格遵循模块化分工：

* **李彦熹 (Li Yanxi)** - 负责工程与接口。搭建底层项目包结构、开发基于生成器的 FASTA IO 读取模块、封装 CLI 命令行接口、完成结果导出逻辑及 README 撰写，并搭建基础单元测试框架。
* **杨茗淞 (Yang Mingsong)** - 负责核心算法。实现基于哈希的 k-mer 索引构建、种子搜索 (seed search)、双向无间隙延伸 (ungapped extension) 算法、动态 score 计算与 top hits 排序机制。
* **高睿勤 (Gao Ruiqin)** - 负责实验评估与报告。构建并验证统计指标 (Identity 与经验 p-value)，准备基准测试数据，设计真实与随机序列的比对对照实验，并完成实验结果分析、数据可视化及最终报告撰写。
