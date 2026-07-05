# BlastMini: A Lightweight Sequence Search Tool

**BlastMini** 是一款纯Python实现的轻量级生物序列相似性检索工具，作为生信算法课程实践项目，完整复刻 NCBI BLAST 核心「种子-延伸」启发式比对逻辑，同时提供命令行交互接口与可视化演示 Notebook。
项目仓库地址：https://github.com/lyx128060/blastmini

---

## 核心特性 (Key Features)
* **高效 k-mer 种子检索**：采用哈希字典实现 O(1) 种子匹配，快速筛选潜在同源序列并开展局部双向延伸；
* **流式大文件读取**：基于生成器分批解析 FASTA 文件，适配百兆级 UniProt 数据库，避免内存溢出；
* **无间隙双向延伸机制**：匹配打分规则为匹配+2、错配-1，内置 X-drop 阈值自动截断低分比对；
* **统计学显著性计算**：通过蒙特卡洛序列打乱模拟随机背景分布，输出经验 P 值评估比对可信度；
* **原生 Python 无额外依赖**：不依赖第三方生信库，普通 Python 环境即可运行；
* **配套可视化 Demo**：Jupyter Notebook 一键完成序列比对，自动生成 TSV 结果表与交互式 SVG/HTML 比对可视化图表。

---

## 仓库目录结构
```text
blastmini/
├── README.md                # 项目总说明文档
├── .gitignore               # 缓存、临时文件、大型数据库屏蔽规则
├── environment.yml          # Conda环境配置文件
├── pyproject.toml           # Python项目打包配置
├── src/                     # 全部核心源码模块
│   ├── __init__.py
│   ├── cli.py               # 命令行交互入口
│   ├── core.py              # k-mer索引构建、种子匹配、双向延伸核心算法
│   ├── io.py                # FASTA文件流式读取、序列格式清洗、结果文件输出
│   └── stats.py             # 序列一致性计算、蒙特卡洛打乱、经验P值统计
├── example/
│   └── blastmini_demo.ipynb # 可视化演示Demo，课程复现推荐入口
└── data/
    ├── README_data.md       # 数据集详细说明文档
    └── test_queries/        # 仓库内置5条测试蛋白序列，可直接运行Demo
        ├── seq_1_glucagon.fasta
        ├── seq_2_ubiquitin.fasta
        ├── seq_3_hemoglobin.fasta
        ├── seq_4_gfp.fasta
        └── seq_5_albumin.fasta

```

---

## 快速开始 (Quick Start)

### 1. 克隆项目代码

```bash
git clone [https://github.com/lyx128060/blastmini.git](https://github.com/lyx128060/blastmini.git)
cd blastmini

```

### 2. 环境配置三种方式

**方式1：临时导入源码路径（通用）**

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"

```

**方式2：本地可编辑安装**

```bash
python3 -m pip install -e .

```

**方式3：Conda一键部署**

```bash
conda env create -f environment.yml
conda activate blastmini

```

### 3. 两种项目运行方式

#### 方式 A：Jupyter 可视化 Demo（课程演示首选）

激活环境后启动笔记，直接读取仓库内置 `test_queries` 测试序列，无需下载完整数据库：

```bash
jupyter notebook example/blastmini_demo.ipynb

```

运行全部单元格后自动生成比对结果 TSV 文件与交互式可视化页面，完整展示序列匹配区间与打分信息。

#### 方式 B：命令行 CLI 全库检索

完整 UniProtKB/Swiss-Prot 参考数据库需自行下载，下载方式详见 `data/README_data.md`，数据库文件放置于 `data/uniprot_sprot.fasta` 后执行：

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"

python3 -m blastmini.cli search \
  data/uniprot_sprot.fasta \
  data/test_queries/seq_1_glucagon.fasta \
  -k 3 \
  --top 5 \
  --seed 1 \
  -o output/seq_1_glucagon_results.tsv

```

**命令行核心参数说明**

* `-k`：k-mer种子长度，蛋白序列推荐3~4，核酸序列推荐4~6
* `--top`：单条查询序列保留的最高分匹配条目数量
* `-o`：TSV格式比对结果输出路径
* `--pvalue-iterations`：P值模拟打乱迭代次数，设为0可关闭统计模块加速调试
* `--seed`：随机数种子，固定数值保证实验结果可复现
* `--x-drop`：局部延伸分数衰减截断阈值

---

## HPC集群批处理脚本（Slurm适配）

适用于多条查询序列批量比对，提交至 CPU 计算队列执行：

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"

# 批量遍历所有内置测试序列完成检索
for query_file in data/test_queries/*.fasta; do
    filename=$(basename -- "$query_file")
    prefix="${filename%.*}"
    
    python3 -m blastmini.cli search data/uniprot_sprot.fasta "$query_file" -k 3 --top 5 -o "output/${prefix}_results.tsv"
done

```

---

## 开发与单元测试 (Development)

项目使用 Python 标准库 `unittest` 实现单元测试，无需额外测试框架：

```bash
export PYTHONPATH="$(pwd)/src:$PYTHONPATH"
python3 -m unittest discover -s tests
python3 -m compileall src tests

```

---

## 工程模块架构 (Architecture)

* `io.py`：流式读取 FASTA 文件、序列格式校验、非法字符过滤、比对结果 TSV 文件写入；
* `core.py`：k-mer 哈希索引构建、种子命中检索、双向无间隙序列延伸、匹配打分核心逻辑；
* `stats.py`：序列全局一致性 Identity 计算、序列随机打乱、经验 P 值显著性评估；
* `cli.py`：命令行参数解析、全局任务流程调度、输入输出路径管理。

---

## 标准输出结果示例 (Example Output)

工具运行完成后输出制表符分隔 TSV 文件，包含完整比对坐标、打分与统计指标：

| query_id | subject_id | score | identity | query_start | query_end | subject_start | subject_end | pvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| query_1 | sp-P0A7V0-RL10_ECOLI | 87 | 95.12% | 1 | 41 | 22 | 62 | 0.0196 |
| query_2 | sp-P01308-INS_HUMAN | 62 | 100.00% | 1 | 31 | 80 | 110 | 0.0196 |

---

## 数据集说明

1. **内置测试数据集（仓库自带，直接可用）**
`data/test_queries` 文件夹包含 5 条不同功能的蛋白 FASTA 序列，文件体积极小，克隆仓库后可直接用于 Demo 演示与代码正确性验证。
2. **完整参考数据库（本地自行下载）**
项目采用人工注释、非冗余的 UniProtKB/Swiss-Prot 蛋白数据库作为全量检索库，完整文件约 280MB，受 GitHub 存储限制未上传至仓库，下载命令与部署步骤详见 `data/README_data.md`。

---

## 开发者分工 (Contributors)

本项目由以下小组成员共同完成，严格遵循模块化分工：

* **李彦熹 (Li Yanxi)** - 负责工程与接口。搭建底层项目包结构、开发基于生成器的 FASTA IO 读取模块、封装 CLI 命令行接口、完成结果导出逻辑及 README 撰写，并搭建基础单元测试框架。
* **杨茗淞 (Yang Mingsong)** - 负责核心算法。实现基于哈希的 k-mer 索引构建、种子搜索 (seed search)、双向无间隙延伸 (ungapped extension) 算法、动态 score 计算与 top hits 排序机制。
* **高睿勤 (Gao Ruiqin)** - 负责实验评估与报告。构建并验证统计指标 (Identity 与经验 p-value)，准备基准测试数据，设计真实与随机序列的比对对照实验，并完成实验结果分析、数据可视化及最终报告撰写。
