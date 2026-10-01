#!/usr/bin/env bash
# ============================================================================
# config.example.sh —— 配置模板。用法：cp config.example.sh config.sh 后按病例修改。
# config.sh 已在 .gitignore 中：它含病例数据路径与样本信息，绝不提交到仓库。
# 所有路径均可指向外置硬盘/云盘，"存储待定"时先跑起来后随时改 WORK/REF 即可。
# ============================================================================

# ---- 0. 计算资源（按本机内存/核数自适应，不再写死某台机器）----------------
# 想手动覆盖：跑脚本前 export THREADS=4 MAXMEM_GB=8 即可（下面用 ${VAR:-默认} 尊重外部值）。
# macOS 用 sysctl，Linux/容器用 nproc + /proc/meminfo；都取不到才回落到 8 核/16GB。
_ncpu=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 8)
_memb=$(sysctl -n hw.memsize 2>/dev/null \
        || awk '/^MemTotal:/{printf "%d", $2*1024}' /proc/meminfo 2>/dev/null || echo 17179869184)
[ -n "$_memb" ] || _memb=17179869184
export HOST_MEM_GB=$(( _memb / 1073741824 ))
export HOST_NCPU=$_ncpu
export THREADS=${THREADS:-$(( _ncpu > 2 ? _ncpu - 2 : 1 ))}     # 留 2 核给系统
# GATK/snpEff Java 堆：物理内存的 ~40%，夹在 [4,32]GB。绝不允许超过物理内存（否则全程 swap）。
_heap=$(( HOST_MEM_GB * 2 / 5 )); [ "$_heap" -lt 4 ] && _heap=4; [ "$_heap" -gt 32 ] && _heap=32
export MAXMEM_GB=${MAXMEM_GB:-$_heap}
# samtools sort：线程数 × 每线程内存，总量控制在物理内存 1/8 左右
export SORT_THREADS=${SORT_THREADS:-4}
_sortmem=$(( HOST_MEM_GB / 8 / SORT_THREADS )); [ "$_sortmem" -lt 1 ] && _sortmem=1
export SORT_MEM_G=${SORT_MEM_G:-$_sortmem}

# ---- 0b. 内存门槛守卫（★换机器必读）--------------------------------------
# bwa-mem2 的 GRCh38 索引本身就需 16.3GB 常驻内存（.bwt.2bit.64 10.1G + .0123 6.2G），
# 叠加 fastp / samtools sort / GATK 堆后，比对与 call 阶段需要 ≥32GB 物理内存。
# 低于门槛的机器（如 16GB Mac mini，2026-07-28 实测）不会报错，只会持续 swap——
# 表现为「极慢甚至跑不动」，这是硬门槛，调 THREADS/MAXMEM 也解决不了。
# 但 04/05 注释筛选与阶段⑤⑥判读只读几十 MB 的 VCF/候选表，任何机器都能跑（见 README「双机分工」）。
export MIN_MEM_GB_COMPUTE=${MIN_MEM_GB_COMPUTE:-32}
require_compute_mem() {
  local step="${1:-本步骤}"                       # 全部用 ${VAR:-默认}，避免调用方 set -u 下炸掉
  local need="${MIN_MEM_GB_COMPUTE:-32}" have="${HOST_MEM_GB:-0}"
  [ "$have" -ge "$need" ] && return 0
  {
    echo ""
    echo "❌ [内存守卫] 本机物理内存 ${have}GB < ${need}GB，拒绝执行：${step}"
    echo "   原因：bwa-mem2 GRCh38 索引需 16.3GB 常驻，GATK 堆另需 ${MAXMEM_GB:-?}GB。"
    echo "   内存不足不会报错，只会持续 swap → 表现为「极慢/卡死」，不是配置问题。"
    echo "   建议：00b/02/03 在 ≥32GB 的机器上跑；本机可直接用 \$RESULTDIR 里已算好的"
    echo "         trio.annot.vcf.gz / candidates.tsv 做 04/05 与阶段⑤⑥判读。"
    echo "   确要强行跑：ALLOW_LOW_MEM=1 bash <脚本>（后果自负，预计数十倍耗时）"
    echo ""
  } >&2
  [ "${ALLOW_LOW_MEM:-0}" = "1" ] || exit 1
  echo "⚠️  ALLOW_LOW_MEM=1 已设置，继续执行（预计极慢）。" >&2
}

# ---- 1. 目录（★存储决策点：把 WORK/REF 指到空间大的盘）--------------------
# WES trio 约需 60-90GB 峰值；WGS 更大。参考库+中间产物建议放外置盘。
# 容器内（run_container.sh）固定挂载为 /refs /data /work /results /logs /tools，
#   对应变量由镜像的环境变量给出，下面的 ${VAR:-默认} 会自动采用。
# 本机直跑：在 shell 里 export SSD_ROOT=/你的外置盘/genomics（或单独 export REFDIR 等）。
# ⛔ 规则：凡"换病例要改"或"调试时可能想临时指向别处"的变量，一律 ${VAR:-默认}，
#    绝不用裸赋值 —— 否则命令行传参（如 `RESULTDIR=/tmp/x bash 98_...`）会被静默覆盖。
export PIPE_ROOT="${TRIO_PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export SSD_ROOT="${SSD_ROOT:-$HOME/genomics}"
export REFDIR="${REFDIR:-$SSD_ROOT/refs}"          # 参考基因组+索引+注释缓存（多病例共用）
export TOOLS_DIR="${TOOLS_DIR:-$SSD_ROOT/tools}"   # ANNOVAR 等需自行注册下载的工具
export CASE_ID="${CASE_ID:-CASE0001_trio}"          # ★换病例必改（建议 送检号_家系规模，不写姓名）
export WORKDIR="${WORKDIR:-$SSD_ROOT/work/$CASE_ID}"       # BAM/CRAM/GVCF 中间产物
export RESULTDIR="${RESULTDIR:-$SSD_ROOT/results/$CASE_ID}"   # 最终 VCF / 候选变异表 / 报告
export LOGDIR="${LOGDIR:-$SSD_ROOT/logs/$CASE_ID}"
mkdir -p "$WORKDIR" "$RESULTDIR" "$LOGDIR"

# ---- 1a. 运行环境（conda 本机 / 容器）------------------------------------
export ENVNAME="${ENVNAME:-trio}"                          # conda 环境名
export CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
# 独立 JDK 环境（本机若 trio 环境里的 java 不可用时用；容器内 java 已在 PATH）
export JAVA_ENV_BIN="${JAVA_ENV_BIN:-$CONDA_ROOT/envs/java/lib/jvm/bin}"
activate_env() {
  # 容器镜像里工具已在 PATH，无需 conda activate
  if [ "${TRIO_IN_CONTAINER:-0}" = "1" ]; then return 0; fi
  # shellcheck disable=SC1091
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  set +u; conda activate "$ENVNAME"; set -u   # conda 激活脚本非 set -u 洁净，需临时关闭
}
# 可用磁盘（GB，整数）—— 用 POSIX `df -Pk`，macOS 与 Linux 通用（BSD 的 df -g 在 Linux 不存在）
df_avail_gb() { df -Pk "$1" 2>/dev/null | awk 'NR==2{printf "%d", $4/1048576}'; }

# ---- 1b. ★ FAMILY_MODE —— 本病例【声明的家系规模】（与 PANEL_SCOPE 同为意图声明）----
#   trio = 父+母+先证者（三样本，经典流程）
#   quad = 父+母+两个子代（四样本）。用于【复发病例】：同一对父母的两次异常妊娠/两个患儿，
#          可做「两胎共有而父母不共有」的交集分析 —— 这是复发病例最核心的判据，
#          三样本流程做不出来。
#   ⛔ 作用同 PANEL_SCOPE：把"本例是几口人"变成机器可核对的声明。00a/98 会硬断言
#      FAMILY_MODE 与实际配置的样本数、fastq 数、VCF 列数一致，防止"配了四个人却跑了三人流程"
#      这类内部自洽、任何交叉校验都不触发的静默失效。
#   singleton = 只有先证者一人（无父母数据）。⚠️ 代价必须在报告里写明：
#          de novo 无法本地验证（PS2/PM6 一律不得赋）、复合杂合无法定相（PM3 降级为"疑似"）、
#          亲缘门与 MCC 无对照样本可跑 —— 这些不是"跑一下就有"的东西，缺样本就是缺证据。
export FAMILY_MODE="${FAMILY_MODE:-trio}"   # ★换病例必改：trio / quad / singleton（须经用户书面确认）
# ★ VCF_PREFIX —— 联合 call 产物的文件名前缀。singleton 用 proband，避免把单人 VCF
#   叫成 trio.*（历史上"文件名与实际内容不符"正是本流水线静默失效的温床）。
if [ "${FAMILY_MODE:-trio}" = "singleton" ]; then
  export VCF_PREFIX="${VCF_PREFIX:-proband}"
else
  export VCF_PREFIX="${VCF_PREFIX:-trio}"
fi


# ============================================================================
# ★★★ 换新病例：从这里往下到第 3 节，是唯一需要改的部分。上面的环境/参考库不动。★★★
# ============================================================================
# ---- 2. 原始数据 -----------------------------------------------------------
# ★ 每份 fastq 分别是谁、哪一次妊娠，必须书面确认，不能靠文件夹名推断。
export SEQ_DATA_ROOT="${SEQ_DATA_ROOT:-$SSD_ROOT/seqdata}"
export DATA_ROOT="${DATA_ROOT:-$SEQ_DATA_ROOT/$CASE_ID/raw}"
# ★ 样本 ID（用于 BAM read-group 与 VCF 列名）
#   ⚠️ 字母序必须 = 先证者 < 母 < 父 < 第二子代（05 按 GEN[0]=先证者 解析），故用 P1_/P2_/P3_/P4_ 前缀。
#   不要用姓名做 ID。
export ID_FETUS="${ID_FETUS:-P1_PROBAND}"   # 先证者（胎儿或患儿）
export ID_MOTHER="${ID_MOTHER:-P2_MOTHER}"  # singleton 时必须置空
export ID_FATHER="${ID_FATHER:-P3_FATHER}"  # singleton 时必须置空
export ID_FETUS1="${ID_FETUS1:-}"           # 仅 quad：第二个子代（如 P4_SIB）
# fastq 路径（R1/R2）
export FQ_FETUS_R1="${FQ_FETUS_R1:-$DATA_ROOT/proband_1.fq.gz}"
export FQ_FETUS_R2="${FQ_FETUS_R2:-$DATA_ROOT/proband_2.fq.gz}"
export FQ_MOTHER_R1="${FQ_MOTHER_R1:-$DATA_ROOT/mother_1.fq.gz}"
export FQ_MOTHER_R2="${FQ_MOTHER_R2:-$DATA_ROOT/mother_2.fq.gz}"
export FQ_FATHER_R1="${FQ_FATHER_R1:-$DATA_ROOT/father_1.fq.gz}"
export FQ_FATHER_R2="${FQ_FATHER_R2:-$DATA_ROOT/father_2.fq.gz}"
export FQ_FETUS1_R1="${FQ_FETUS1_R1:-}"; export FQ_FETUS1_R2="${FQ_FETUS1_R2:-}"
# ★ 线粒体专项文库（仅当厂商另交付时；供 smt_call_mito.sh 用，不进核 WES 流程）
export FQ_MITO_R1="${FQ_MITO_R1:-}"
export FQ_MITO_R2="${FQ_MITO_R2:-}"
export ID_MITO="${ID_MITO:-P1_PROBAND_MT}"

# 原始数据在位检查（只有 00b/02 需要 fastq；供这两个脚本调用，避免跑到一半才发现路径错）
require_fastq() {
  local missing=0 f
  local list="$FQ_FATHER_R1 $FQ_FATHER_R2 $FQ_MOTHER_R1 $FQ_MOTHER_R2 $FQ_FETUS_R1 $FQ_FETUS_R2"
  # quad 模式额外要求第二个子代的两个 fastq
  if [ "${FAMILY_MODE:-trio}" = "quad" ]; then
    list="$list ${FQ_FETUS1_R1:-} ${FQ_FETUS1_R2:-}"
  fi
  # singleton 模式：只有先证者一人（FQ_FETUS_* 槽位），不要求父母 fastq
  if [ "${FAMILY_MODE:-trio}" = "singleton" ]; then
    list="$FQ_FETUS_R1 $FQ_FETUS_R2"
  fi
  for f in $list; do
    [ -f "$f" ] || { echo "❌ [数据守卫] fastq 不存在: $f" >&2; missing=1; }
  done
  if [ "$missing" = 1 ]; then
    echo "   请确认外置盘已挂载（当前 SEQ_DATA_ROOT=$SEQ_DATA_ROOT），并核对第 2 节的路径。" >&2
    echo "   当前 FAMILY_MODE=${FAMILY_MODE:-trio}（quad 需要 8 个 fastq，trio 需要 6 个）。" >&2
    exit 1
  fi
}

# 本病例参与分析的全部样本 ID（trio 三个 / quad 四个），供各脚本统一遍历
all_sample_ids() {
  if [ "${FAMILY_MODE:-trio}" = "quad" ]; then
    echo "$ID_FETUS $ID_MOTHER $ID_FATHER ${ID_FETUS1:-}"
  elif [ "${FAMILY_MODE:-trio}" = "singleton" ]; then
    echo "$ID_FETUS"
  else
    echo "$ID_FETUS $ID_MOTHER $ID_FATHER"
  fi
}

# ---- 3. 家系与性别（PED 编码：1=男 2=女 0=未知/无数据）-------------------
# 先证者性别由 02b 的 somalier（或 chrX/chrY 覆盖比）复核后回填，XL 分析必需。
export SEX_FATHER="${SEX_FATHER:-1}"
export SEX_MOTHER="${SEX_MOTHER:-2}"
export SEX_FETUS="${SEX_FETUS:-0}"    # ★比对后回填
export SEX_FETUS1="${SEX_FETUS1:-0}"  # 仅 quad

# ---- 4. 数据类型（★由 00b_detect_datatype 判定后回填）----------------------
export ASSAY="${ASSAY:-WES}"            # WES / WGS
export CAPTURE_BED="${CAPTURE_BED:-}"   # 厂商 target.bed；未提供则退化为 CDS±10bp 通用区间

# ---- 5. 参考基因组版本 ----------------------------------------------------
export GENOME_BUILD="GRCh38"
export REF_FASTA="$REFDIR/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# ---- 5b. ANNOVAR + dbNSFP（功能预测注释；参考侧共享，一次性）-------------
# 用户放的是 ANNOVAR 格式 dbNSFP 4.7a（hg38_dbnsfp47a.txt[.idx]），故注释走 ANNOVAR 而非 SnpSift。
export ANNOVAR_DIR="${ANNOVAR_DIR:-$TOOLS_DIR/annovar}"        # table_annovar.pl 所在（ANNOVAR 需自行注册下载，不进镜像）
export ANNOVAR_HUMANDB="$REFDIR/dbsnfp"                        # hg38_dbnsfp47a.txt + .idx 所在目录
export DBNSFP_PROTOCOL="dbnsfp47a"                             # 对应 hg38_${DBNSFP_PROTOCOL}.txt
export ANNOVAR_BUILD="hg38"                                    # GRCh38 对应 ANNOVAR 的 hg38

# ---- 5c. v2 增强：MCC / somalier / SpliceAI（参考侧共享 + 每病例阈值）------
export MCC_THRESHOLD=0.05        # 母源细胞污染预警阈值（VerifyBamID2 FREEMIX 或信息位点 alt 比例）
export SOMALIER_SITES="$REFDIR/GRCh38/somalier.sites.GRCh38.vcf.gz"   # 01c 下载（小，~几MB）
export SPLICEAI_SNV_VCF="$REFDIR/GRCh38/spliceai_scores.masked.snv.hg38.vcf.gz"     # 01c 下载（~29GB，可选）
export SPLICEAI_INDEL_VCF="$REFDIR/GRCh38/spliceai_scores.masked.indel.hg38.vcf.gz" # 01c 下载（~1GB，可选）

# ---- 6. 表型驱动分析 ------------------------------------------------------
# ⚠️ 2026-08-06 修（某病例踩到）：原为无条件赋值，命令行 `PANEL_FILE=xxx bash 05_...` 传入的值
#    会在 source config.sh 时被静默覆盖 → 05 用了上个病例的 panel，输出一批与本例表型无关的基因，
#    且不报错。凡"换病例要改"的变量一律改成 ${VAR:-默认} 形式，尊重外部覆盖。
export PANEL_FILE="${PANEL_FILE:-$PIPE_ROOT/scripts/panel_ALLGENES.txt}"   # 候选基因清单
# ★ PANEL_SCOPE —— 本病例【声明的扫描范围】，换病例必须显式确认。
#   allgenes = 全外显子不分层（PANEL_FILE 须为全基因清单，>5000 基因）
#   panel    = 表型定向 panel（PANEL_FILE 为小 panel，须与本病例表型对应）
#   作用：把"用哪个 panel"这一意图变成机器可核对的声明。
#   起因（2026-08-07）：某病例误用上一例的 CAKUT+先心 panel，产出的 5 条候选
#   与该 panel 完全自洽，任何交叉校验都不会触发——只有把"本例应扫全外显子"这个
#   意图写下来，脚本才有办法判定它跑错了。00a/98 会据此硬断言。
export PANEL_SCOPE="${PANEL_SCOPE:-allgenes}"
# ★ PANEL_RATIONALE —— 一句话写明为何选这个范围（换病例必改，供人工复核与报告留痕）
export PANEL_RATIONALE="${PANEL_RATIONALE:-<示例> 一句话写明本例为何选 allgenes / panel（表型谱、既往检测、漏诊代价等）}"
export HPO_FILE="${HPO_FILE:-$PIPE_ROOT/scripts/phenotype_HPO.example.txt}"   # ★换病例必改（建议放在病例目录，勿提交到仓库）
export MAX_AF=0.001    # 罕见变异群体频率阈值（gnomAD popmax），AD/AR 可分别覆写

echo "[config] ASSAY=$ASSAY  本机=${HOST_NCPU}核/${HOST_MEM_GB}GB  THREADS=$THREADS  MAXMEM_GB=$MAXMEM_GB  sort=${SORT_THREADS}×${SORT_MEM_G}G"
echo "[config] FAMILY_MODE=$FAMILY_MODE  PANEL_SCOPE=${PANEL_SCOPE:-<未声明>}"
echo "[config] REFDIR=$REFDIR  WORKDIR=$WORKDIR  DATA_ROOT=$DATA_ROOT"
# 注意：本文件被 `set -e` 的脚本 source，最后一条命令必须返回 0，否则会连累调用方退出。
if [ "${HOST_MEM_GB:-0}" -lt "$MIN_MEM_GB_COMPUTE" ]; then
  echo "[config] ⚠️ 本机内存 ${HOST_MEM_GB}GB < ${MIN_MEM_GB_COMPUTE}GB：00b/02/03（比对与 call）会被守卫拦下；04/05/判读可正常跑。"
fi
true
