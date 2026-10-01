# 家系 Trio WES/WGS 基因组分析：技能 + 流水线 + 容器

本仓库包含三部分：

| 目录 / 文件 | 内容 |
|---|---|
| `skill/prenatal-trio-genome-analysis/` | Claude 技能（编排层）：六阶段流程、两处强制核对关卡、疑难病例扩大挖掘。打包成 zip 即可上传到 claude.ai |
| `pipeline/scripts/` | 生信流水线脚本（00 环境 → 01 参考库 → 02 比对 → 03 call → 04 注释 → 05 遗传模型筛选 → 98 量级断言 → 99 溯源清单），含 quad / singleton / 线粒体分支 |
| `pipeline/docs/` | SOP：家系 WES 病例分析标准流程与踩坑记录 |
| `pipeline/README.md` | 流水线详细说明：硬件门槛、运行顺序、双机分工 |
| `environment.yml` | 工具版本锁定（本机 conda 与容器共用同一份） |
| `Dockerfile` / `run_container.sh` | 容器化：镜像只装工具，参考库/数据/结果全部以卷挂载 |

> ⛔ **本仓库不含任何病例数据**：fastq/BAM/CRAM/VCF、报告、病例 config、HPO 文件、参考库都不进仓库（见 `.gitignore`）。
> 病例配置请从 `pipeline/scripts/config.example.sh` 复制为 `config.sh` 后修改，`config.sh` 已被忽略。

## 一、本机直跑（macOS / Linux，conda）

```bash
cd pipeline/scripts
cp config.example.sh config.sh          # 填 CASE_ID、FAMILY_MODE、样本 ID、fastq 路径、性别、panel/HPO
export SSD_ROOT=/你的外置盘/genomics      # 参考库在 $SSD_ROOT/refs，中间产物/结果按 CASE_ID 隔离
bash 00_setup_env.sh                    # 首次：装 miniforge + conda 环境 trio（按 environment.yml）
bash 01_download_refs.sh                # 首次：参考库（数小时）
bash 01b_download_snpeff.sh
bash 01c_download_spliceai.sh
bash 00a_preflight.sh                   # 每病例第一个跑：配置断言
bash 00b_detect_datatype.sh             # WES / WGS 判定
bash 02_align.sh fetus                  # 各样本分别比对（father / mother / fetus）
bash 02b_qc_somalier.sh
bash 03_call_trio.sh
bash 03b_denovo_refine.sh
bash 04_annotate.sh
bash 05a_make_allgenes_panel.sh
bash 05_inheritance_filter.sh
bash 98_sanity_check.sh                 # 出报告前必跑：量级断言
bash 99_manifest.sh
```

quad / singleton 的运行顺序见 `pipeline/docs/SOP-家系WES病例分析流程.md` 末节。

## 二、容器运行（Linux 服务器 / 环境复现）

构建镜像（支持 linux/amd64 与 linux/arm64）：

```bash
docker build -t trio-wes:latest .
```

多架构一次构建：

```bash
docker buildx build --platform linux/amd64,linux/arm64 -t trio-wes:latest .
```

跑某一步（参考库、数据、结果均为挂载卷）：

```bash
REFS=/data/genomics/refs DATA=/data/seqdata/CASE0001/raw OUT=/data/genomics CASE_ID=CASE0001_trio CONFIG=$PWD/config.sh ./run_container.sh 02_align.sh father
```

| 变量 | 挂载到 | 说明 |
|---|---|---|
| `REFS` | `/refs`（只读） | 参考库 |
| `DATA` | `/data`（只读） | 本病例原始 fastq |
| `OUT` | `/out` | 其下自动建 `work/` `results/` `logs/`，按 `CASE_ID` 隔离 |
| `TOOLS` | `/tools`（只读，可选） | ANNOVAR（需自行注册下载，不进镜像） |
| `CASE` | `/case`（只读，可选） | 本病例 HPO / panel 文件；config 里写 `/case/...` |
| `CONFIG` | `/pipeline/scripts/config.sh`（可选） | 不给则用模板 |

`EXTRA_ENV="THREADS MAXMEM_GB"` 可把当前 shell 的同名变量透传进容器。路径里有空格时先建软链接再传入。

## 三、注意事项

- **内存门槛**：比对与 call（00b/02/03）需要 ≥32GB 物理内存，不足时脚本会主动拦下。Docker Desktop 需在设置里把 VM 内存调到 ≥32GB。
- **参考版本统一 GRCh38**，全链（基因组 / known-sites / snpEff / ClinVar / dbNSFP / gnomAD）必须同一版本。
- **不装 VEP**：ensembl-vep 在 arm64 段错误，且其 perl 依赖会把 samtools 拖回 0.1.19；注释改用 snpEff+SnpSift + ANNOVAR-dbNSFP + SpliceAI 预计算 VCF。
- **CRAM 版本**：samtools 1.22+ 默认写 CRAM 3.1，GATK 4.6 只认 3.0，02 已显式指定 `version=3.0`。
- **脚本保持 bash 3.2 兼容**（macOS 自带 bash），不使用关联数组。
- 所有变异结论须经有资质实验室复核与遗传咨询后方可用于临床决策。
