# skmer-smk2 使用说明

`skmer-smk2` 是一个面向双端测序 FASTQ 数据的系统发育分析流程。它把 FASTQ 质控、可选参考序列过滤、读段合并、按碱基量抽样、Skmer 距离计算、WASTER 建树、Mash 距离热图与 bootstrap 共识树整合到一个 Snakemake 工作流中，并提供统一命令行入口：

```bash
skmer-smk2 run -i FASTQ_DIR -ref REF_FASTA -s 75 -j 48
```

本仓库同时包含一些辅助脚本，用于在正式运行前扫描/修复 FASTQ 文件、处理头部截断 FASTQ、备份并替换已修复数据，以及在 HPC 集群上提交任务。

## 适用场景

该流程适合用于从多个样本的 paired-end FASTQ 文件构建样本间系统发育关系。典型输入是每个样本一对 R1/R2 文件，输出包括：

- 过滤和抽样后的 FASTQ
- 每个样本的 reads/bases 统计表
- Skmer 距离矩阵和系统发育树
- Mash 距离矩阵、系统发育树和距离热图
- WASTER 系统发育树

如果提供 `-ref` 参考序列，流程会先用 `bowtie2` 过滤掉能比对到该参考的 reads，例如去除叶绿体/质体基因组 reads；如果不提供 `-ref`，则跳过参考过滤步骤。

## 仓库结构

```text
.
|-- README.md
|-- skmer_snakemake_project/
|   |-- pyproject.toml
|   |-- README.md
|   |-- raw_data/
|   `-- src/skmer_smk2/
|       |-- cli.py
|       |-- templates/
|       `-- workflow/
|           |-- Snakefile
|           `-- scripts/
|-- repair_head_truncated_fastq.sh
|-- replace_repaired_fastq.sh
|-- scan_repair_fastq.sh
|-- run_skmer.py
|-- skmer_hpc.sh
|-- snakefile
`-- scripts/
```

推荐优先使用 `skmer_snakemake_project/` 中打包好的 `skmer-smk2` 命令。根目录下的 `snakefile`、`run_skmer.py` 和 `scripts/` 是较早的本地运行形式，可作为参考，但日常运行建议使用安装后的命令行工具。

## 软件组成

主程序：

- `skmer-smk2 run`：运行完整 Snakemake 工作流
- `skmer-smk2 doctor`：检查外部依赖是否可用
- `skmer-smk2 init`：导出 HPC/修复脚本模板
- `skmer-smk2 repair-fastq`：复制或运行 FASTQ 扫描修复脚本

核心流程：

1. 自动识别输入目录中的 paired-end FASTQ
2. 使用 `fastp` 做质控和接头处理
3. 可选：使用 `bowtie2` 去除比对到参考序列的 reads
4. 使用 `repair.sh` 修复过滤后 R1/R2 配对关系
5. 使用 `bbmerge.sh` 合并 overlapping reads，并保留未合并 reads
6. 统计每个样本总 bases，根据 `-s` 指定的百分位确定统一抽样深度
7. 为 Skmer、WASTER 和 Mash 准备输入
8. 运行 Skmer、Mash、WASTER 并生成树文件和热图

## 依赖环境

Python 包本身只安装 `skmer-smk2` 命令和 Snakemake 工作流文件；大型生信软件需要安装在当前环境或可被 `PATH` 找到。

必须或常用依赖：

```text
python >= 3.8
snakemake
fastp
bowtie2
bowtie2-build
bbmap: repair.sh, bbmerge.sh
skmer
fastme
raxmlHPC
waster
mash
seqkit
gzip
```

其中多数工具可通过 conda/mamba 安装；`skmer` 和 `waster` 通常需要按各自上游说明手动安装，并保证命令在 `PATH` 中可用。

## 安装方法

### 方法一：从本地源码安装

在仓库根目录执行：

```bash
cd skmer_snakemake_project
python -m pip install .
```

开发调试时可以使用 editable 安装：

```bash
cd skmer_snakemake_project
python -m pip install -e .
```

安装后检查命令：

```bash
skmer-smk2 --version
skmer-smk2 -h
```

### 方法二：从 GitHub 安装

```bash
python -m pip install git+https://github.com/fandunjin/skmer_smk2.git
```

如果 HPC 登录节点访问 GitHub SSL 不稳定，可以用 zip 安装：

```bash
python -m pip install --no-cache-dir --force-reinstall \
  https://github.com/fandunjin/skmer_smk2/archive/refs/heads/main.zip
```

### 方法三：安装或补齐 conda 依赖

如果已有环境，例如 `01bio`：

```bash
source /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh
conda activate 01bio
```

先检查依赖：

```bash
skmer-smk2 doctor
```

严格检查：

```bash
skmer-smk2 doctor --strict
```

自动安装 conda/bioconda 可用的软件：

```bash
skmer-smk2 doctor --install
```

也可以手动安装常见依赖：

```bash
mamba install -c conda-forge -c bioconda \
  snakemake fastp bowtie2 bbmap fastme raxml mash seqkit gzip
```

然后再单独安装 `skmer` 和 `waster`，并确认：

```bash
which skmer
which waster
```

## 输入文件要求

输入目录必须直接包含 FASTQ 文件，不能传入仓库根目录或上一级目录。

支持的双端文件命名：

```text
SampleA_1.fq.gz        SampleA_2.fq.gz
SampleB_R1.fq.gz       SampleB_R2.fq.gz
SampleC.R1.fastq.gz    SampleC.R2.fastq.gz
SampleD-R1.fq.gz       SampleD-R2.fq.gz
```

样本名由 R1/R2 后缀前面的部分决定。例如：

```text
H_asiatica_SAMC1020836_1.fq.gz
H_asiatica_SAMC1020836_2.fq.gz
```

样本名为：

```text
H_asiatica_SAMC1020836
```

运行前建议确认每个样本都有 R1 和 R2：

```bash
ls *_1.fq.gz *_2.fq.gz
```

## 快速运行

不使用参考序列过滤：

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 --printshellcmds
```

使用参考序列过滤：

```bash
skmer-smk2 run \
  -i /path/to/fastq_dir \
  -ref /path/to/refDNA.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

参数说明：

```text
-i, --input              FASTQ 输入目录
-ref, --ref              可选参考基因组 FASTA，用于 bowtie2 过滤
-s, --sample-percentile  按样本总 bases 排序后选取的百分位，默认 75
-j, --jobs               Snakemake 可用核心数/任务数
-b, --bootstraps         bootstrap 重复次数，默认 100
--exclude-samples        跳过指定样本，逗号或空格分隔
--workdir                结果输出目录，默认当前目录
--dry-run                只预演，不实际运行
--printshellcmds         打印 Snakemake 执行的 shell 命令
```

预演命令：

```bash
skmer-smk2 run -i /path/to/fastq_dir -ref /path/to/refDNA.fasta -s 75 -j 1 --dry-run
```

传递额外 Snakemake 参数：

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 -- --keep-going
```

## 示例数据

打包项目中包含小型示例数据：

```text
skmer_snakemake_project/raw_data/raw_data/
skmer_snakemake_project/raw_data/ref.fna
```

示例 dry-run：

```bash
cd skmer_snakemake_project
skmer-smk2 run -i raw_data/raw_data -ref raw_data/ref.fna -s 75 -j 1 -b 2 --dry-run
```

## HPC 集群运行

在集群上，推荐写一个提交脚本，只负责申请资源、激活环境，然后运行同一个 `skmer-smk2 run` 命令。

示例 `jsub` 脚本：

```bash
#!/bin/bash
#JSUB -J skmer
#JSUB -n 48
#JSUB -q normal
#JSUB -o skmer_log.%J
#JSUB -e skmer_err.%J
#JSUB -cwd .

source /hpcfile/users/92024286/anaconda3/etc/profile.d/conda.sh
conda activate 01bio

skmer-smk2 run \
  -i /hpcfile/users/92024286/Huperzia \
  -ref /path/to/refDNA.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

提交：

```bash
jsub < skmer_hpc.sh
```

如果集群使用 `sbatch` 或 `qsub`，只需要替换脚本头部的调度器参数，主体命令不变。

## FASTQ 检查和修复

正式运行前建议先检查 FASTQ 是否完整、是否有 malformed record、R1/R2 是否能配对。

### 通用扫描修复脚本

复制模板：

```bash
skmer-smk2 repair-fastq --workdir . --copy-only
```

运行：

```bash
bash scan_repair_fastq.sh /path/to/fastq_dir
```

输出：

```text
repaired_fastq/fastq_repair_report.tsv
repaired_fastq/input_for_skmer/
```

报告字段中的常见状态：

```text
UNCHANGED      原始文件可用，未生成修复副本
REPAIRED       seqkit sana 修复后记录数发生变化
NEED_REUPLOAD  gzip 流损坏或截断，需要重新下载/上传
FAILED         修复失败，需要查看日志
```

后续可以把 `repaired_fastq/input_for_skmer/` 作为 `skmer-smk2 run -i` 的输入目录。

### 头部截断 FASTQ 的成对修复

本仓库根目录提供了专门脚本：

```bash
repair_head_truncated_fastq.sh
```

它用于处理某些文件末端或记录头截断的 FASTQ。脚本会：

1. 对每个样本分别运行 `seqkit sana`
2. 使用 `seqkit pair` 重新提取成对 reads
3. 检查输出 gzip 是否正常
4. 检查 R1/R2 reads 是否大于 0 且完全相等
5. 将可用输出写入 `repaired_for_skmer/`
6. 将空输出或坏输出隔离到 `failed_repair_outputs/`
7. 写出 `fastq_repair_report.tsv`

运行：

```bash
bash repair_head_truncated_fastq.sh
```

或提交作业：

```bash
jsub < repair_head_truncated_fastq.sh
```

检查报告：

```bash
cat fastq_repair_report.tsv
ll repaired_for_skmer
ll failed_repair_outputs
```

状态解释：

```text
REPAIRED          修复成功，可用于后续分析
EMPTY_OUTPUT      R1/R2 无可配对 reads，不能用于 paired-end 分析
UNPAIRED_OUTPUT   R1/R2 reads 数不相等，不能直接用于 paired-end 分析
BAD_OUTPUT        gzip 或 FASTQ 统计异常
MISSING_INPUT     输入文件不存在或为空
```

如果 `seqkit sana` 后 reads 很多，但 `seqkit pair` 得到 0 pairs，通常说明 R1/R2 的 read ID 没有交集，可能是 mate 文件下载错、上传错或来自不同数据集。这种情况不能靠普通 repair 恢复，应该回源重新下载对应 FASTQ。

### 用修复成功文件替换原始文件

对于确认 `REPAIRED` 的样本，可以使用：

```bash
replace_repaired_fastq.sh
```

该脚本会：

1. 只替换脚本中列出的成功样本
2. 替换前检查 gzip 和 reads 配对数
3. 创建备份目录 `original_fastq_backup_YYYYMMDD_HHMMSS`
4. 将原始 FASTQ 移入备份目录
5. 将 `repaired_for_skmer/` 中的修复文件复制回原文件名

运行：

```bash
bash replace_repaired_fastq.sh
```

检查：

```bash
ll original_fastq_backup_*
ll *.fq.gz
```

## 主要输出结果

运行完成后，结果通常位于当前工作目录的 `results/` 下。

统计结果：

```text
results/stats/post_filter_summary.sorted.tsv
results/stats/head_summary.sorted.tsv
results/stats/head_base_cutoff.txt
```

每个样本的中间结果：

```text
results/<sample>/clean/
results/<sample>/nDNA/
results/<sample>/nDNAOK/<sample>.fq
results/<sample>/stats/
```

Skmer 输出：

```text
results/skmer/dimtrx_main.txt
results/skmer/dimtrx_main_cor_.txt
results/skmer/tree.direct.tre
results/skmer/tree.bootstrap.tre
results/skmer/tree.merged.tre
```

WASTER 输出：

```text
results/waster/input.tsv
results/waster/waster.tree
```

Mash 输出：

```text
results/mash/distances.tsv
results/mash/distances.phy
results/mash/distance_heatmap.svg
results/mash/tree.direct.tre
results/mash/tree.bootstrap.tre
results/mash/tree.merged.tre
```

## 常见问题

### 1. `ADDR2LINE: unbound variable`

原因通常是脚本启用了 `set -u`，而 conda 的某些 `activate.d` 脚本访问了未定义变量。解决方法是在 `conda activate` 前后临时关闭 nounset：

```bash
set +u
source /path/to/conda.sh
conda activate 01bio
set -u
```

本仓库中的修复脚本已经包含该处理。

### 2. `No paired FASTQ files found`

检查 `-i` 是否指向真正存放 FASTQ 的目录，并确认命名符合：

```text
sample_1.fq.gz / sample_2.fq.gz
sample_R1.fq.gz / sample_R2.fq.gz
sample.R1.fq.gz / sample.R2.fq.gz
sample-R1.fq.gz / sample-R2.fq.gz
```

### 3. `seqkit pair` 输出 0 paired-end reads

如果 R1/R2 文件都有大量 reads，但配对结果为 0，说明两端 read ID 不匹配。常见原因：

- R1/R2 不是同一个样本
- 下载或上传时 mate 文件错配
- 数据源本身给的是非配对文件
- read name 被预处理工具改写且两端不一致

这种情况不能直接作为 paired-end 数据运行，应重新获取正确的 R1/R2。

### 4. 结果目录已有旧文件

Snakemake 会复用已完成结果。若需要完全重跑，建议先新建工作目录，或谨慎清理旧的 `results/` 和 `.skmer_smk2_workflow/`。不要在不确认路径的情况下批量删除数据。

### 5. 升级后仍使用旧 Snakefile

如果 `--printshellcmds` 中还出现类似 `python scripts/fastq_stats_and_sample.py` 的旧路径，说明旧 workflow cache 仍在使用。可以删除缓存：

```bash
rm -rf .skmer_smk2_workflow
```

然后重新运行：

```bash
skmer-smk2 run -i /path/to/fastq_dir -s 75 -j 48 --printshellcmds
```

## 推荐完整流程

1. 激活环境并检查依赖：

```bash
conda activate 01bio
skmer-smk2 doctor
```

2. 检查并修复 FASTQ：

```bash
bash repair_head_truncated_fastq.sh
cat fastq_repair_report.tsv
```

3. 只保留 `REPAIRED` 或原本正常的 paired-end 样本。

4. 运行主流程：

```bash
skmer-smk2 run \
  -i /path/to/fastq_dir \
  -ref /path/to/refDNA.fasta \
  -s 75 \
  -j 48 \
  --printshellcmds
```

5. 查看最终树和统计结果：

```bash
ls results/skmer/*.tre
ls results/mash/*.tre
ls results/waster/*.tree
cat results/stats/head_summary.sorted.tsv
```
