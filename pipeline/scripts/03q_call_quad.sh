#!/usr/bin/env bash
# ============================================================================
# 03q_call_quad.sh —— 四人家系联合变异检测（父+母+两胎）
# 本例特有：复发病例家系有两次异常妊娠的胎儿数据，四人联合 call 才能做
#   「两胎共有 / 两胎均不共有」的交集分析（复发病例的核心判据）。
# 产物：
#   quad.filtered.vcf.gz (P1_F2,P2_MO,P3_FA,P4_F1) —— 交集分析 / P4 基因型查询用
#   trio.filtered.vcf.gz (P1_F2,P2_MO,P3_FA)       —— 由 quad 抽样本得到，喂给既有 03b/04/05 链
#   quad.ped / trio.ped
# ★ 按 runbook 2026-08-22 教训：HaplotypeCaller 用 -Xmx8g（比 25g 快一倍）且四样本并行；
#   先写 tmp_hc/ 并做记录数断言，通过后才 mv 入正式路径（杜绝半成品被下游当成品）。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
require_compute_mem "03q_call_quad.sh"
activate_env
mkdir -p "$RESULTDIR"; LOG="$LOGDIR/03q_call.log"; exec > >(tee -a "$LOG") 2>&1

L=""; [ -n "${CAPTURE_BED:-}" ] && L="-L $CAPTURE_BED" || \
      { [ "$ASSAY" = "WES" ] && L="-L $REFDIR/GRCh38/exons_cds.bed"; }
ALL="$ID_FETUS $ID_MOTHER $ID_FATHER $ID_FETUS1"

# ---- 1) 每样本 GVCF（并行，小堆）------------------------------------------
pids=""
for SID in $ALL; do
  G="$WORKDIR/$SID/${SID}.g.vcf.gz"
  [ -f "$G" ] && { echo "[call] $SID GVCF 已存在，跳过"; continue; }
  ( set -e
    T="$WORKDIR/$SID/tmp_hc"; mkdir -p "$T"
    gatk --java-options "-Xmx8g" HaplotypeCaller \
      -R "$REF_FASTA" -I "$WORKDIR/$SID/${SID}.final.cram" \
      -O "$T/${SID}.g.vcf.gz" -ERC GVCF $L
    n=$(bcftools view -H "$T/${SID}.g.vcf.gz" | wc -l | tr -d ' ')
    if [ "$n" -lt 100000 ]; then echo "❌ [$SID] GVCF 记录数仅 $n (<10万)，判为半成品，保留在 $T" >&2; exit 1; fi
    mv "$T/${SID}.g.vcf.gz" "$G"; mv "$T/${SID}.g.vcf.gz.tbi" "$G.tbi"; rmdir "$T"
    echo "[call] ✅ $SID GVCF 完成，$n 条记录"
  ) > "$LOGDIR/03q_hc_$SID.log" 2>&1 &
  pids="$pids $!"
done
fail=0; for p in $pids; do wait "$p" || fail=1; done
[ "$fail" = 0 ] || { echo "❌ 有 HaplotypeCaller 子任务失败，见 $LOGDIR/03q_hc_*.log"; exit 1; }

# ---- 2) 四人合并 + 联合定型 ------------------------------------------------
COMBI="$WORKDIR/quad.combined.g.vcf.gz"; RAW="$WORKDIR/quad.raw.vcf.gz"
if [ ! -f "$RAW" ]; then
  VARGS=""; for SID in $ALL; do VARGS="$VARGS -V $WORKDIR/$SID/${SID}.g.vcf.gz"; done
  gatk --java-options "-Xmx${MAXMEM_GB}g" CombineGVCFs -R "$REF_FASTA" $VARGS -O "$COMBI"
  gatk --java-options "-Xmx${MAXMEM_GB}g" GenotypeGVCFs -R "$REF_FASTA" -V "$COMBI" -O "$RAW"
fi

# ---- 3) 硬过滤 -------------------------------------------------------------
QFILT="$RESULTDIR/quad.filtered.vcf.gz"
gatk VariantFiltration -R "$REF_FASTA" -V "$RAW" -O "$WORKDIR/quad.tagged.vcf.gz" \
  --filter-expression "QD<2.0" --filter-name QD2 \
  --filter-expression "FS>60.0" --filter-name FS60 \
  --filter-expression "MQ<40.0" --filter-name MQ40 \
  --filter-expression "SOR>3.0" --filter-name SOR3 \
  --filter-expression "MQRankSum<-12.5" --filter-name MQRS \
  --filter-expression "ReadPosRankSum<-8.0" --filter-name RPRS
bcftools view -f PASS,. "$WORKDIR/quad.tagged.vcf.gz" -Oz -o "$QFILT"; tabix -f -p vcf "$QFILT"

# ---- 4) 由 quad 抽出标准 trio VCF（喂既有 03b/04/05 链，样本列顺序 P1<P2<P3）----
TFILT="$RESULTDIR/trio.filtered.vcf.gz"
bcftools view -s "$ID_FETUS,$ID_MOTHER,$ID_FATHER" "$QFILT" -Oz -o "$TFILT"; tabix -f -p vcf "$TFILT"

# ---- 5) PED ---------------------------------------------------------------
{ echo -e "#FID\tIID\tPID\tMID\tSEX\tPHENO"
  echo -e "FAM1\t$ID_FATHER\t0\t0\t$SEX_FATHER\t1"
  echo -e "FAM1\t$ID_MOTHER\t0\t0\t$SEX_MOTHER\t1"
  echo -e "FAM1\t$ID_FETUS\t$ID_FATHER\t$ID_MOTHER\t$SEX_FETUS\t2"; } > "$RESULTDIR/trio.ped"
{ cat "$RESULTDIR/trio.ped"
  echo -e "FAM1\t$ID_FETUS1\t$ID_FATHER\t$ID_MOTHER\t$SEX_FETUS1\t2"; } > "$RESULTDIR/quad.ped"

echo "[call] ✅ quad: $QFILT  ($(bcftools view -H "$QFILT" | wc -l | tr -d ' ') 条)"
echo "[call] ✅ trio: $TFILT  ($(bcftools view -H "$TFILT" | wc -l | tr -d ' ') 条)"
echo "[call] 样本列顺序: $(bcftools view -h "$QFILT" | tail -1 | cut -f10-)"
echo "下一步: bash 03b_denovo_refine.sh → bash 04_annotate.sh → bash 05_inheritance_filter.sh"
