#!/usr/bin/env bash
# ============================================================================
# 02_align.sh —— 单样本: fastp质控 → bwa-mem2比对 → 排序 → 标记重复 → BQSR → CRAM
# 用法: bash 02_align.sh <father|mother|fetus|fetus1>
#       trio 跑前三个；FAMILY_MODE=quad 时再跑 fetus1（第二个子代）。可并行，见内存说明。
# ★ 全程边算边清中间文件、监控磁盘。
# ★ 内存门槛：本步骤需 ≥32GB 物理内存（bwa-mem2 索引常驻 16.3GB），低于门槛会被 config.sh
#   的 require_compute_mem 拦下。并行开三个终端会让索引占用×3，仅在 ≥64GB 机器上这么做。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
require_compute_mem "02_align.sh（bwa-mem2 比对）"   # 内存不足直接拦下，见 config.sh 第 0b 节
require_fastq                                        # 原始 fastq 路径在位检查
activate_env
who="${1:?用法: bash 02_align.sh father|mother|fetus|fetus1}"

case "$who" in
  father) SID=$ID_FATHER; R1=$FQ_FATHER_R1; R2=$FQ_FATHER_R2;;
  mother) SID=$ID_MOTHER; R1=$FQ_MOTHER_R1; R2=$FQ_MOTHER_R2;;
  fetus)  SID=$ID_FETUS;  R1=$FQ_FETUS_R1;  R2=$FQ_FETUS_R2;;
  fetus1) # 第二个子代，仅 FAMILY_MODE=quad 时可用
          [ "${FAMILY_MODE:-trio}" = "quad" ] || { echo "❌ FAMILY_MODE=${FAMILY_MODE:-trio}，未声明 quad，拒绝比对 fetus1。
   若本例确为四人家系，请在 config.sh 设 FAMILY_MODE=quad 并填 ID_FETUS1/FQ_FETUS1_*。" >&2; exit 1; }
          SID=${ID_FETUS1:?FAMILY_MODE=quad 但 ID_FETUS1 未设}
          R1=${FQ_FETUS1_R1:?FQ_FETUS1_R1 未设}; R2=${FQ_FETUS1_R2:?FQ_FETUS1_R2 未设};;
  *) echo "未知样本 $who"; exit 1;;
esac
O="$WORKDIR/$SID"; mkdir -p "$O"; LOG="$LOGDIR/02_$SID.log"; exec > >(tee -a "$LOG") 2>&1
echo "===== [$SID] $(date) 开始比对 ====="
dfree(){ echo "[磁盘] $(df_avail_gb "$WORKDIR")GB可用"; }; dfree

# 1) fastp —— MGI/DNBSEQ 数据，PE 重叠自动检测接头；输出 QC 报告
CLEAN1="$O/clean_1.fq.gz"; CLEAN2="$O/clean_2.fq.gz"
if [ ! -f "$CLEAN1" ]; then
  fastp -i "$R1" -I "$R2" -o "$CLEAN1" -O "$CLEAN2" \
    --detect_adapter_for_pe --thread "$THREADS" --qualified_quality_phred 15 \
    --length_required 50 --json "$O/${SID}.fastp.json" --html "$O/${SID}.fastp.html"
fi

# 2) 比对 + 排序（read-group 里写样本名，trio 联合 call 必需）
RG="@RG\tID:${SID}\tSM:${SID}\tPL:DNBSEQ\tLB:${SID}"
SORT="$O/${SID}.sorted.bam"
if [ ! -f "$SORT" ]; then
  bwa-mem2 mem -t "$THREADS" -R "$RG" "$REF_FASTA" "$CLEAN1" "$CLEAN2" \
    | samtools sort -@ "$SORT_THREADS" -m "${SORT_MEM_G}G" -o "$SORT" -
  samtools index "$SORT"
  rm -f "$CLEAN1" "$CLEAN2"          # 清理 fastp 中间文件
fi; dfree

# 3) 标记重复
MD="$O/${SID}.md.bam"
if [ ! -f "$MD" ]; then
  gatk --java-options "-Xmx${MAXMEM_GB}g" MarkDuplicates \
    -I "$SORT" -O "$MD" -M "$O/${SID}.dupmetrics.txt" --CREATE_INDEX true
  rm -f "$SORT" "$SORT.bai"
fi; dfree

# 4) BQSR
REC="$O/${SID}.recal.table"; DBSNP=$(ls "$REFDIR/GRCh38/"*dbsnp*.vcf.gz 2>/dev/null|head -1)
MILLS=$(ls "$REFDIR/GRCh38/"Mills*.vcf.gz 2>/dev/null|head -1)
FINAL_CRAM="$O/${SID}.final.cram"
if [ ! -f "$FINAL_CRAM" ]; then
  gatk --java-options "-Xmx${MAXMEM_GB}g" BaseRecalibrator \
    -I "$MD" -R "$REF_FASTA" --known-sites "$DBSNP" --known-sites "$MILLS" -O "$REC"
  gatk --java-options "-Xmx${MAXMEM_GB}g" ApplyBQSR \
    -I "$MD" -R "$REF_FASTA" --bqsr-recal-file "$REC" -O "$O/${SID}.bqsr.bam"
  # 转 CRAM 省空间（比 BAM 小 ~40%）
  # ★ 显式写 CRAM 3.0：samtools 1.22+ 默认写 3.1，GATK 4.6 的 htsjdk 只认到 3.0（不指定会导致 03 步 HaplotypeCaller 报错）
  samtools view -@4 -C --output-fmt-option version=3.0 -T "$REF_FASTA" -o "$FINAL_CRAM" "$O/${SID}.bqsr.bam"
  samtools index "$FINAL_CRAM"
  rm -f "$MD" "$O/${SID}.md.bai" "$O/${SID}.bqsr.bam" "$O/${SID}.bqsr.bai"
fi

# 5) 覆盖度 QC（同时给 00b 用来判定 WES/WGS）
mosdepth -t "$THREADS" -n -f "$REF_FASTA" --by "${CAPTURE_BED:-$REFDIR/GRCh38/exons_cds.bed}" \
  "$O/${SID}" "$FINAL_CRAM" 2>/dev/null || true
echo "===== [$SID] ✅ 完成 → $FINAL_CRAM ====="; dfree
