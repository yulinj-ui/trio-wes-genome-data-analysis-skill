#!/usr/bin/env bash
# ============================================================================
# s03_call_single.sh —— 单人（singleton）变异检测
# HaplotypeCaller(GVCF) → GenotypeGVCFs → 硬过滤 → $RESULTDIR/${VCF_PREFIX}.filtered.vcf.gz
#
# ⛔ 与 03_call_trio.sh 的本质区别（不是"少两个样本"这么简单，报告必须写明）：
#    · 没有父母基因型 → de novo 不可判定，03b 的 CalculateGenotypePosteriors/PossibleDeNovo
#      全部不适用（它们以 PED 三人为前提），故本脚本【不】产 hiConf/loConfDeNovo 标签。
#    · 复合杂合只能"同基因两个杂合"提示，无法定相 → 一律标"疑似、相位未验证"。
#    · 亲缘门无对照样本 → 只能做性别 + 自身指纹，报告顶部须标"家系关系未验证"。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
[ "${FAMILY_MODE:-trio}" = "singleton" ] || {
  echo "❌ FAMILY_MODE=${FAMILY_MODE:-trio}，本脚本仅用于 singleton。三/四人家系请跑 03_call_trio.sh / 03q_call_quad.sh" >&2; exit 1; }
require_compute_mem "s03_call_single.sh（GATK HaplotypeCaller，Java 堆 ${MAXMEM_GB}g）"
activate_env
mkdir -p "$RESULTDIR"; LOG="$LOGDIR/s03_call.log"; exec > >(tee -a "$LOG") 2>&1
VP="${VCF_PREFIX:-proband}"; SID="$ID_FETUS"

L=""; [ -n "${CAPTURE_BED:-}" ] && L="-L $CAPTURE_BED" || \
      { [ "$ASSAY" = "WES" ] && L="-L $REFDIR/GRCh38/exons_cds.bed"; }

CRAM="$WORKDIR/$SID/${SID}.final.cram"
[ -f "$CRAM" ] || { echo "❌ 找不到 $CRAM，先跑 02_align.sh fetus" >&2; exit 1; }

# 1) GVCF（小堆更快：runbook 2026-08-22 教训，-Xmx8g 比 -Xmx25g 快一倍）
G="$WORKDIR/$SID/${SID}.g.vcf.gz"
if [ ! -f "$G" ]; then
  T="$WORKDIR/$SID/tmp_hc"; mkdir -p "$T"
  gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R "$REF_FASTA" -I "$CRAM" -O "$T/${SID}.g.vcf.gz" -ERC GVCF $L
  n=$(bcftools view -H "$T/${SID}.g.vcf.gz" | wc -l | tr -d ' ')
  [ "$n" -ge 100000 ] || { echo "❌ GVCF 记录数仅 $n (<10万)，判为半成品，保留在 $T" >&2; exit 1; }
  mv "$T/${SID}.g.vcf.gz" "$G"; mv "$T/${SID}.g.vcf.gz.tbi" "$G.tbi"; rmdir "$T"
  echo "[call] ✅ GVCF 完成，$n 条记录"
fi

# 2) 定型
RAW="$RESULTDIR/${VP}.raw.vcf.gz"
[ -f "$RAW" ] || gatk --java-options "-Xmx${MAXMEM_GB}g" GenotypeGVCFs -R "$REF_FASTA" -V "$G" -O "$RAW"

# 3) 硬过滤（与 03_call_trio.sh 同一套阈值，保持跨病例可比）
FILT="$RESULTDIR/${VP}.filtered.vcf.gz"
gatk VariantFiltration -R "$REF_FASTA" -V "$RAW" -O "$WORKDIR/${VP}.tagged.vcf.gz" \
  --filter-expression "QD<2.0" --filter-name QD2 \
  --filter-expression "FS>60.0" --filter-name FS60 \
  --filter-expression "MQ<40.0" --filter-name MQ40 \
  --filter-expression "SOR>3.0" --filter-name SOR3 \
  --filter-expression "MQRankSum<-12.5" --filter-name MQRS \
  --filter-expression "ReadPosRankSum<-8.0" --filter-name RPRS
bcftools view -f PASS,. "$WORKDIR/${VP}.tagged.vcf.gz" -Oz -o "$FILT"; tabix -p vcf "$FILT"

NV=$(bcftools view -H "$FILT" | wc -l | tr -d ' ')
echo "[call] ✅ 单人 VCF: $FILT （$NV 个 PASS 变异）"
[ "$NV" -ge 20000 ] || echo "⚠️ 变异数 $NV 低于 WES 经验区间下界(2万)，先怀疑静默失效再下结论" >&2
echo "[call] ⛔ 本例无 PED、无 de novo 判定 —— 跳过 03b_denovo_refine.sh（它以三人 PED 为前提）"
echo "下一步: bash 04_annotate.sh"
