# 📂 blastmini 数据目录 (Data Directory)

本目录用于存放 `blastmini` 运行所需的测试序列、参考数据库以及比对过程中生成的索引和临时文件。

为了遵守 GitHub 的存储限制并保证代码库的轻量化，**本目录下的所有大型数据库文件和生成的索引文件均已被 `.gitignore` 忽略，请勿强制上传。**

---

## 目录结构说明
完整运行项目后，data目录结构如下：
- `README_data.md` : 当前数据说明文档
- `test_queries/` : **[随代码库提供]** 用于Demo演示的精简FASTA测试蛋白序列，仅KB级大小，已上传仓库，克隆后可直接运行测试
- `uniprot_sprot.fasta` : **[需本地自行下载]** 完整 UniProtKB/Swiss-Prot 手工注释蛋白参考数据库
- `swissprot_db.index` : **[本地运行生成]** blastmini建库流程产生的索引缓存文件
- `output_results.tsv` : **[本地运行生成]** 序列比对完成后输出的结果表格

---

## 数据获取指南
### 1. Test Queries（快速复现Demo）
`test_queries/` 文件夹内置精简测试数据集，无需额外下载，可直接用来验证工具全部比对、可视化功能，满足课程Demo演示需求。

### 2. 完整参考数据库 UniProtKB/Swiss-Prot
本项目采用人工注释、非冗余的Swiss-Prot蛋白数据库作为全量检索库；完整文件约280MB，无法上传GitHub，本地下载命令：
```bash
# 进入data目录执行下载
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# 解压文件
gunzip uniprot_sprot.fasta.gz
