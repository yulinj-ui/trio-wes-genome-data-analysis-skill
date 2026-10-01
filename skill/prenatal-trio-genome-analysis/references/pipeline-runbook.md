# 生信流水线运行手册(阶段③④用)

> 配套脚本在本仓库 `pipeline/scripts/`(配置模板 `config.example.sh` → 复制为 `config.sh`,后者含病例信息,已 gitignore)。
> 容器化运行见文末「容器化运行」节。
> 本手册是运行顺序 + **已踩坑排错表**。所有坑都已在脚本里修复,新病例直接复用;但重跑/新机器/报错时对照本表能快速定位。

---

## 运行顺序

| 步骤 | 脚本 | 作用 | 一次性? |
|---|---|---|---|
| 0 | `00_setup_env.sh` | 装 miniforge + conda 环境 `trio` | ✅ 一次性 |
| 1 | `01_download_refs.sh` | GRCh38 + bwa-mem2 索引 + known-sites + CDS bed | ✅ 一次性 |
| 1b | `01b_download_snpeff.sh` | snpEff 库 + ClinVar + slivar-gnomAD(+dbNSFP 需手填链接) | ✅ 一次性 |
| 0b | `00b_detect_datatype.sh` | 抽样比对判 WES/WGS | 每病例 |
| 2 | `02_align.sh father\|mother\|fetus` | fastp→bwa-mem2→标记重复→BQSR→CRAM | 每病例(每样本) |
| 3 | `03_call_trio.sh` | HaplotypeCaller×3→联合定型→硬过滤→PED | 每病例 |
| 4 | `04_annotate.sh` | snpEff + SnpSift(dbNSFP/ClinVar) | 每病例 |
| 5 | `05_inheritance_filter.sh` | AD/de novo·AR纯合·AR复合杂合·XL 筛选 + panel 分层 | 每病例 |

**新病例只改 `config.sh`**:CASE_ID、FAMILY_MODE、样本路径(FQ_*)、样本 ID、SEX_*、CAPTURE_BED、PANEL_SCOPE/PANEL_FILE/PANEL_RATIONALE、HPO_FILE。REFDIR/SSD_ROOT/WORKDIR/RESULTDIR 均为 `${VAR:-默认}`,本机用 export 指定、容器内由 run_container.sh 注入。环境与参考库不动。

**长步骤(01 下载、02 比对、03 call)放后台跑 + 挂 Monitor**,零 token;完成/报错自动通知。

---

## 已踩坑排错表(按报错信号查)

| 报错信号 | 根因 | 修复 |
|---|---|---|
| `CMAKE_PREFIX_PATH: unbound variable` | conda 激活脚本与 `set -u` 冲突 | 激活用 `set +u; conda activate trio; set -u` 包裹(脚本已改) |
| `line NN: dbsnp: unbound variable` 之类 | **macOS 系统 bash 3.2 不支持关联数组 `declare -A`** | 改用普通数组(脚本已改)。这是 macOS 通病,写任何脚本都避开 declare -A |
| `samtools ... unrecognized command 'dict'` / 各种旧式报错 | **conda 把 samtools 装成了古董版 0.1.19**(被 VEP 相关 perl 包 `perl-bio-samtools` 硬依赖拖回) | `conda remove perl-bio-samtools ensembl-vep perl-bio-db-hts` 后 `conda install "samtools>=1.20"` → 1.22.x |
| `samtools dict` 不存在 | 旧 samtools 无该子命令 | 改用 `gatk CreateSequenceDictionary` |
| VEP 启动 `exit 139` 段错误 | **ensembl-vep 在 Apple Silicon(arm64) bioconda 版段错误** | **弃用 VEP,改 snpEff+SnpSift**(arm64 原生 Java,稳)。注释链见 04 脚本 |
| GATK HaplotypeCaller `RuntimeException: CRAM version 3.1 is not supported` | **samtools 1.22+ 默认写 CRAM 3.1,GATK 4.6 的 htsjdk 只认 3.0** | `samtools view -C --output-fmt-option version=3.0 ...`(02 脚本已加)。已生成的 CRAM 转码即可,不必重新比对 |
| 下 known-sites `curl (56) error 403` | 旧 GCS 桶 `genomics-public-data` 已失效 | 正确桶 `https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/`,用 `.gz`+`.tbi` |
| slivar gnomAD 库 404 | URL 错 | `https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip`(~4GB) |
| snpEff `OutOfMemoryError: Java heap space` | 默认堆太小,加载 GRCh38.p14 库时 OOM | `export _JAVA_OPTIONS="-Xmx${MAXMEM_GB}g"`(04 脚本已加) |
| SnpSift extractFields `INFO field 'dbNSFP_xxx' not found in VCF header` | 请求了未注入的字段(dbNSFP 未就绪时) | 05 脚本已改为**动态探测**字段是否存在,缺失则不请求 |
| dbNSFP 下载反复 `curl (18)`/`(56)`/`(33)` 中断 | 50GB 学术版单连接服务器不稳、断点续传支持时好时坏 | `curl --http1.1 -C - --retry 100 --retry-all-errors --speed-time 60 --speed-limit 200`;多次重试;下完 `bgzip -t` 验完整性再用。dbNSFP 是可选增强,缺它 04 会降级为 snpEff+ClinVar |
| mosdepth 读 CRAM `ERROR: specify a reference file` | CRAM 解码需参考 | mosdepth 加 `-f "$REF_FASTA"`(02 脚本已加) |

---

## 若资产不在位(新机器 / 重建)

按 `00→01→01b` 顺序重建。关键约束:
- **参考版本统一 GRCh38**,全链(基因组/known-sites/snpEff库/ClinVar/dbNSFP/slivar-gnomAD)必须同一版本。
- conda 用 miniforge(00 脚本按 uname 自动选 macOS/Linux、arm64/x86_64 安装包);工具版本统一锁定在仓库根 `environment.yml`,本机与容器同一份。
- 存储:参考库 ~40GB + 中间产物,建议外置盘(APFS 格式化,支持符号链接/POSIX 权限,exFAT 会出诡异错误)。做成多病例共享参考库,新病例只改 WORK/RESULT 指向。
- dbNSFP(学术版,提供 CADD/REVEL/AlphaMissense/MetaRNN + gnomAD4.1 频率)需去 dbnsfp.org 申请下载链接,填 `DBNSFP_URL` 环境变量。字段名对齐 5.3.1a 实测表头:**无 SpliceAI**;gnomAD 频率字段是 `gnomAD4.1_joint_AF`(非 `gnomAD_exomes_AF`)。

---

## 关键实操提醒(阶段⑤会用到)

- **VCF 样本列顺序 = GATK 按样本名字母排序**,不一定是"父母子"。解读基因型前先 `bcftools view -h ... | tail -1` 核实列顺序,再对应到 PED。
- **胎儿/先证者性别判定**:`samtools coverage -r chrX/chrY` 看平均深度,chrX≈常染色体基线且 chrY≈0 → 女(XX);chrX≈半 + chrY 有覆盖 → 男(XY)。回填 config.sh 的 SEX_ 和 PED。
- **既往报告身份核验**:报告变异若为 GRCh37,用 `https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38/<NM>:<c.变异>/all` 按 cDNA 命名拿 GRCh38 坐标(比 liftOver chain 可靠),再 `bcftools view -r chr:pos` 去 VCF 核对基因型。
- **多等位位点(GT 含 2)** 的遗传模型判断不可全信工具,人工看原始 reads 复核(HLA 区尤其是比对伪影热点)。

---

## 容器化运行(Linux 服务器 / 换机器时)

- 构建:`docker build -t trio-wes:latest .`(多架构:`docker buildx build --platform linux/amd64,linux/arm64 ...`)。构建期会自检 java/samtools/bcftools/bwa-mem2/gatk/snpEff 等,任何一个起不来就构建失败。
- 运行:`REFS=… DATA=… OUT=… CASE_ID=… CONFIG=…/config.sh ./run_container.sh 02_align.sh father`
  - 挂载:`REFS→/refs`(只读)、`DATA→/data`(只读)、`OUT→/out`(其下按 CASE_ID 建 work/results/logs)、可选 `TOOLS→/tools`(ANNOVAR)、`CASE→/case`(本例 HPO/panel)
  - 容器内 `TRIO_IN_CONTAINER=1`:`activate_env` 跳过 conda 激活,工具已在 PATH
  - `config.sh` 若写了 HPO_FILE/PANEL_FILE,容器里要写成 `/case/...` 路径
- **参考库不进镜像**(~200GB),ANNOVAR 需自行注册下载,也不进镜像。
- ⚠️ 内存门槛不变:比对/call 仍需 ≥32GB;Docker Desktop(macOS)要在设置里把 VM 内存调到 ≥32GB,否则 config 读到的是 VM 内存,02/03 会被守卫拦下。
- ⚠️ macOS 上 Docker 跑 linux/arm64 镜像有虚拟化开销且 CRAM/fastq 走 bind mount I/O 较慢;本机 M 系列 Mac 仍建议直接用 conda 跑,容器主要用于 Linux 服务器与环境复现。
