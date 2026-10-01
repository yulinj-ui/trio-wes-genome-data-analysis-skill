#!/usr/bin/env bash
# ============================================================================
# 02c_king_kinship.sh —— KING-robust 亲缘系数（somalier 亲缘门的独立复核）
# 背景(2026-08-24 某四人家系实测)：somalier 的 relatedness 在本数据上有系统性偏移
#   —— 四对亲子实测 0.72–0.78(应 0.5)、父母对 0.587(应 0)，触发 98 的亲缘门误报。
# KING-robust 不依赖等位基因频率面板、对群体结构稳健，是该门的独立仲裁手段。
# 用法: bash 02c_king_kinship.sh [VCF]   默认用 quad.filtered.vcf.gz
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
VCF="${1:-$RESULTDIR/quad.filtered.vcf.gz}"
[ -f "$VCF" ] || { echo "❌ 缺 VCF: $VCF"; exit 1; }
OUT="$RESULTDIR/qc/king_kinship.txt"; mkdir -p "$RESULTDIR/qc"
TMP="$WORKDIR/_king_gt.tsv"
bcftools view -m2 -M2 -v snps -t ^chrX,chrY,chrM "$VCF" -Ou \
  | bcftools query -f '%CHROM\t%POS[\t%GT\t%DP\t%GQ]\n' > "$TMP"
python3 "$DIR/king_kinship.py" "$TMP" "$(bcftools query -l "$VCF" | tr '\n' ',')" | tee "$OUT"
rm -f "$TMP"
echo "[02c] ✅ 结果已落盘: $OUT"
