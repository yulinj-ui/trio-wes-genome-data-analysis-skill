# 家系 WES/WGS 遗传分析流水线

FASTQ → 比对 → 联合 call → 注释 → 表型驱动的 AD/AR-复合杂合/XL 筛选 → 排序 → ACMG 致病性判读 → 中文报告。
支持 trio（父+母+先证者）、quad（父+母+两个子代，复发病例）、singleton（仅先证者）三种家系规模（`FAMILY_MODE`）。

## 分工：Claude = 编排者 + 判读者（不亲自算比对）
- **计算层（01–05）**：标准生信工具（bwa-mem2/GATK/**snpEff+SnpSift**/ANNOVAR/slivar/somalier），Claude 写脚本+调度+报错自愈
  - 注：注释原设计用 VEP，但 **VEP 在 Apple Silicon(arm64) 启动段错误，已彻底弃用**，改用 arm64 原生的 snpEff+SnpSift（基因后果+ClinVar）+ **ANNOVAR-dbNSFP 4.7a**（CADD/REVEL/AlphaMissense）+ **SpliceAI 预计算 VCF**（剪接）
- **判读层（06）**：Claude + 技能 `gene-variant-pathogenicity` / `franklin-variant-lookup` / `prenatal-report-interpreter` 做排序、ACMG、报告

## ⛔ 硬件门槛与双机分工（换机器前必读）

**比对与变异检测需要 ≥32GB 物理内存，这是硬门槛。** bwa-mem2 的 GRCh38 索引本身就要 **16.3GB** 常驻
（`.bwt.2bit.64` 10.1G + `.0123` 6.2G），叠加 fastp/sort/GATK 堆后 32GB 才够用。
**内存不足不会报错，只会持续 swap → 表现为「极慢甚至跑不动」**，调 `THREADS`/`MAXMEM_GB` 解决不了。
> 📌 2026-07-28 在一台 **16GB Mac mini** 上实测：conda 环境能装完，比对根本无法开始。已据此加内存守卫。

| 阶段 | 脚本 | 内存要求 |
|---|---|---|
| 比对 / 变异检测 | `00b` · `02` · `03` | **≥32GB（必需，脚本会主动拦下）** |
| 注释 / 遗传模型筛选 | `04` · `05` | 8–16GB（ANNOVAR 吃磁盘 I/O 不吃内存） |
| 阶段⑤⑥ 判读 | 无脚本（模型层） | 任何机器（输入仅几 KB–20MB） |

**双机分工**：参考库和中间产物都在共享外置盘上，两台机器可以接力——
**计算机**（≥32GB）跑 `00`→`05`；**判读机**（16GB 够用）只读盘上的 `results/candidates/*.tsv` 和
`trio.annot.vcf.gz` 做阶段⑤⑥，**不需要装 conda 环境，也不需要参考库**。

**计算资源已自适应**：`config.sh` 按 `sysctl hw.memsize`（macOS）或 `/proc/meminfo`（Linux/容器）与核数自动算 `THREADS`、`MAXMEM_GB`、
`samtools sort` 内存，不再写死某台机器。要手动覆盖就在跑脚本前 `export THREADS=4 MAXMEM_GB=8`。
强行在低内存机上跑：`ALLOW_LOW_MEM=1 bash 02_align.sh father`（后果自负）。

## 存储

参考库 + 中间产物 + 原始测序数据建议全部放外置盘（APFS/ext4，**不要 exFAT**）：
- `$SSD_ROOT/{refs,work,results,logs,tools}` —— refs 为**多病例共享**参考库（~214GB）
- 原始 fastq 放 `$SEQ_DATA_ROOT/<CASE_ID>/raw/`
- 容器内固定挂载为 `/refs`、`/data`、`/out`（见仓库根 README）

## 运行顺序
```bash
cd scripts
cp config.example.sh config.sh   # 按病例修改（config.sh 不进仓库）
bash 00a_preflight.sh         # 开跑前配置断言（每病例第一个跑）
bash 00_setup_env.sh          # 装 miniforge + conda 环境（~15min，唯一装在本机的东西）
bash 01_download_refs.sh      # GRCh38 + bwa-mem2 索引 + known-sites（一次性，数小时）
bash 01b_download_snpeff.sh   # snpEff 库 + ClinVar + slivar-gnomAD（一次性）
bash 01c_download_spliceai.sh # somalier sites + SpliceAI 预计算 VCF（一次性；注意 chr 前缀转换）
bash 00b_detect_datatype.sh   # 抽样判定 WES/WGS → 回填 config.sh 的 ASSAY/CAPTURE_BED
bash 02_align.sh father        # 三样本各跑（≥64GB 才建议并行开三终端，索引占用会×3）
bash 02_align.sh mother
bash 02_align.sh fetus
bash 02b_qc_somalier.sh       # somalier 亲缘/性别/指纹 + VerifyBamID2 母源污染(MCC)
bash 03_call_trio.sh          # 联合 call → trio VCF + PED
bash 03b_denovo_refine.sh     # de novo 后验精修 → hiConf/loConfDeNovo 分层
bash 04_annotate.sh           # snpEff + ANNOVAR-dbNSFP + ClinVar + SpliceAI
bash 05_inheritance_filter.sh # AD/de novo · AR纯合 · AR复合杂合 · XL 候选表
bash 99_manifest.sh           # 生成溯源清单 manifest.json（每次跑完/重跑后都要更新）
# 06：回到 Claude 对话，按 06_interpret_with_claude.md 做排序+ACMG+报告
# 07_benchmark_giab.sh：GIAB HG002 trio 基准，仅换机/改流程后做临床级验证时跑
```
> **换新病例只改 `config.sh` 的 CASE_ID / FAMILY_MODE 与第 2–6 节**：数据路径、样本 ID、性别、ASSAY、panel/HPO。
> 第 0 节（计算资源）自适应、第 1 节（参考库路径）跨病例共用，用环境变量覆盖即可。

## 新病例开跑前的四件确认事
1. **这台机器内存够不够**（`sysctl -n hw.memsize`）——<32GB 就别从 fastq 开始跑，见上面「硬件门槛」。
2. **数据↔个体对应**——每份 fastq 分别是谁、哪一次妊娠，**必须书面确认，不能靠文件夹名推断**
   （历史教训：目录名写着"第一胎"，却被默认当成当前妊娠，大半程分析对象搞错）。
3. **WES vs WGS**（`00b`）——决定磁盘方案与是否需 CAPTURE_BED；顺带回填 `ASSAY`。
4. **胎儿性别**——`02b` 的 somalier 判定（或 `03` 后看 chrY 覆盖）回填 `SEX_FETUS`，XL 分析必需。
   另：**捕获 kit 的 target.bed** 若测序公司能提供，填入 `CAPTURE_BED` 让 WES on-target 更准。

## 预估耗时（WES trio，M1 Max/64GB 实测）
环境+参考下载 半天（一次性）｜比对 3×约1–2h｜call+注释+筛选 约1–2h。首次全程约 1 个工作日。
> 32GB 机器请串行跑 `02`，总时长约 ×1.5；<32GB 机器跑不了这一段（见「硬件门槛」）。

## 合规与免责
本流程为科研/临床辅助分析，所有变异结论须经**有资质实验室**复核与遗传咨询后方可用于临床决策；
致病性判读遵循 ACMG/AMP 指南，报告须显式声明分析盲区（见 `06` 表）。
