#!/usr/bin/env bash
# ============================================================================
# 03_call_trio.sh —— 家系联合变异检测（trio joint calling）
# HaplotypeCaller(GVCF)×3 → CombineGVCFs → GenotypeGVCFs → 过滤 → 生成 PED
# trio 联合 call 的价值: 直接支持 de novo 检测与孟德尔错误(ME)标注。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
require_compute_mem "03_call_trio.sh（GATK HaplotypeCaller，Java 堆 ${MAXMEM_GB}g）"
activate_env
mkdir -p "$RESULTDIR"; LOG="$LOGDIR/03_call.log"; exec > >(tee -a "$LOG") 2>&1

# WES 限定在捕获/CDS 区间可大幅提速
L=""; [ -n "${CAPTURE_BED:-}" ] && L="-L $CAPTURE_BED" || \
      { [ "$ASSAY" = "WES" ] && L="-L $REFDIR/GRCh38/exons_cds.bed"; }

# 1) 每样本 GVCF —— 并行 + 小堆
# ★ runbook 2026-08-22 教训：HaplotypeCaller 用 -Xmx8g 比 -Xmx25g 快一倍（大堆 → G1GC 停顿变长、
#   缓存局部性变差；实测 RSS 仅约 1.2GB）；且 arm64 上 Intel GKL 原生库加载失败、PairHMM 退化为
#   单线程纯 Java，故三样本必须并行，否则白白浪费核。
# ★ 先写 tmp_hc/ 并做记录数断言，通过后才 mv 入正式路径 —— 杜绝"半成品 GVCF 被下游当成品"。
pids=""
for who in father mother fetus; do
  case "$who" in father)SID=$ID_FATHER;;mother)SID=$ID_MOTHER;;fetus)SID=$ID_FETUS;;esac
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
  ) > "$LOGDIR/03_hc_$SID.log" 2>&1 &
  pids="$pids $!"
done
fail=0; for p in $pids; do wait "$p" || fail=1; done
[ "$fail" = 0 ] || { echo "❌ 有 HaplotypeCaller 子任务失败，见 $LOGDIR/03_hc_*.log"; exit 1; }

# 2) 合并 + 联合定型
COMBI="$WORKDIR/trio.combined.g.vcf.gz"; RAW="$RESULTDIR/trio.raw.vcf.gz"
gatk --java-options "-Xmx${MAXMEM_GB}g" CombineGVCFs -R "$REF_FASTA" \
  -V "$WORKDIR/$ID_FATHER/${ID_FATHER}.g.vcf.gz" \
  -V "$WORKDIR/$ID_MOTHER/${ID_MOTHER}.g.vcf.gz" \
  -V "$WORKDIR/$ID_FETUS/${ID_FETUS}.g.vcf.gz" -O "$COMBI"
gatk --java-options "-Xmx${MAXMEM_GB}g" GenotypeGVCFs \
  -R "$REF_FASTA" -V "$COMBI" -O "$RAW"

# 3) 过滤（trio/WES 规模用硬过滤即可，稳健）
FILT="$RESULTDIR/trio.filtered.vcf.gz"
gatk VariantFiltration -R "$REF_FASTA" -V "$RAW" -O "$WORKDIR/trio.tagged.vcf.gz" \
  --filter-expression "QD<2.0" --filter-name QD2 \
  --filter-expression "FS>60.0" --filter-name FS60 \
  --filter-expression "MQ<40.0" --filter-name MQ40 \
  --filter-expression "SOR>3.0" --filter-name SOR3 \
  --filter-expression "MQRankSum<-12.5" --filter-name MQRS \
  --filter-expression "ReadPosRankSum<-8.0" --filter-name RPRS
bcftools view -f PASS,. "$WORKDIR/trio.tagged.vcf.gz" -Oz -o "$FILT"; tabix -p vcf "$FILT"

# 4) PED（家系文件；胎儿性别若未回填=0，后续 XL 分析前需确定）
PED="$RESULTDIR/trio.ped"
{ echo -e "#FID\tIID\tPID\tMID\tSEX\tPHENO"
  echo -e "FAM1\t$ID_FATHER\t0\t0\t$SEX_FATHER\t1"
  echo -e "FAM1\t$ID_MOTHER\t0\t0\t$SEX_MOTHER\t1"
  echo -e "FAM1\t$ID_FETUS\t$ID_FATHER\t$ID_MOTHER\t$SEX_FETUS\t2"; } > "$PED"

echo "[call] ✅ trio VCF: $FILT ；PED: $PED"
[ "$SEX_FETUS" = "0" ] && echo "⚠️ 胎儿性别未定，XL 分析前先: samtools idxstats 看 chrY 覆盖，回填 config.sh 的 SEX_FETUS"
echo "下一步: bash 04_annotate.sh"
