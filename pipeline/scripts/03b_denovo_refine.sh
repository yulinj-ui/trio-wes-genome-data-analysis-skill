#!/usr/bin/env bash
# ============================================================================
# 03b_denovo_refine.sh —— de novo 后验精修（压假阳性）
# 裸孟德尔过滤(0/0,0/0,0/1)会因父母某位点低GQ/低覆盖造出假 de novo。
# GATK CalculateGenotypePosteriors(带PED做家系后验) → VariantAnnotator PossibleDeNovo,
# 产物带 hiConfDeNovo / loConfDeNovo INFO 标签,供 05 优先取 hiConf。
# 在 03_call_trio 之后、04_annotate 之前运行。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
IN="$RESULTDIR/trio.filtered.vcf.gz"; PED="$RESULTDIR/trio.ped"
PP="$WORKDIR/trio.pp.vcf.gz"; OUT="$RESULTDIR/trio.dn.vcf.gz"

[ -f "$IN" ] || { echo "❌ 缺 $IN，请先跑 03_call_trio.sh"; exit 1; }
[ -f "$PED" ] || { echo "❌ 缺 PED"; exit 1; }

# 1) 家系后验重定型（用 trio 结构 + 群体先验修正基因型，压低质量假阳）
echo "[03b] CalculateGenotypePosteriors …"
gatk --java-options "-Xmx${MAXMEM_GB}g" CalculateGenotypePosteriors \
  -V "$IN" -ped "$PED" --skip-population-priors -O "$PP"

# 2) 标注可能的 de novo（PossibleDeNovo 依 GQ/后验给 hiConf/loConf 分层）
echo "[03b] VariantAnnotator PossibleDeNovo …"
gatk --java-options "-Xmx${MAXMEM_GB}g" VariantAnnotator \
  -V "$PP" -A PossibleDeNovo -ped "$PED" -O "$OUT"
tabix -f -p vcf "$OUT" 2>/dev/null || true

echo "[03b] ✅ 完成: $OUT"
echo "     hiConf de novo 数: $(bcftools view -H -i 'INFO/hiConfDeNovo!="."' "$OUT" 2>/dev/null | wc -l | tr -d ' ')"
echo "     loConf de novo 数: $(bcftools view -H -i 'INFO/loConfDeNovo!="."' "$OUT" 2>/dev/null | wc -l | tr -d ' ')"
echo "下一步: bash 04_annotate.sh（04 会优先注释本步的 trio.dn.vcf.gz）"
