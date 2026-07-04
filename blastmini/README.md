# blastmini

[![Python Version](https://img.shields.io/badge/python-3.8%x2B-blue.svg)]()
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

`blastmini` 是一个用于生物序列相似性搜索的轻量级 Python 工具包。本项目旨在实现类似于 BLAST 的核心流程，包括 k-mer 索引构建、seed 匹配、无间隙延伸以及显著性统计评估。

## 安装指南

建议在虚拟环境中安装：

```bash
git clone <你的 GitHub/Gitee 仓库链接>
cd blastmini
pip install -e .