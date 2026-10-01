#!/usr/bin/env bash
# ============================================================================
# 05q_quad_shared.sh —— 四人家系专项分析（在 05 之后跑）
# 复发病例的核心判据：两个子代共有而父母不共有的变异。三样本流程做不出来。
# 依次执行：
#   ① quad_shared_models.py  六种交集模型（X连锁/共有de novo/AR纯合/复合杂合/母源/父源显性）
#   ② pheno_reverse_search.py 表型优先反向搜索 + 全基因组 ClinVar P/LP + HIGH impact 扫描
#   ③ carrier_scan_3C.py      3-C 夫妇共同携带（不看子代基因型的独立扫描）
# ⛔ 三者的输出都只是【待核实清单】，写进报告前必须逐条回原始 BAM 过伪影筛（见 runbook）。
# 用法: bash 05q_quad_shared.sh [表型定向panel文件]
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env

[ "${FAMILY_MODE:-trio}" = "quad" ] || { echo "❌ FAMILY_MODE=${FAMILY_MODE:-trio}，本脚本仅用于 quad。" >&2; exit 1; }
QUAD="$RESULTDIR/quad.filtered.vcf.gz"
FULL="$WORKDIR/all_variants_full.tsv"
OUT="$RESULTDIR/candidates"; mkdir -p "$OUT"
PHENO_PANEL="${1:-$DIR/panel_pheno_${CASE_ID}.txt}"

[ -f "$QUAD" ] || { echo "❌ 缺 $QUAD，请先跑 03q_call_quad.sh" >&2; exit 1; }
[ -f "$FULL" ] || { echo "❌ 缺 $FULL，请先跑 05_inheritance_filter.sh" >&2; exit 1; }

# 样本列顺序硬断言：quad_shared_models.py 按 先证者/母/父/第二子代 解析四列
EXPECT="$(printf '%s\n' "$ID_FETUS" "$ID_MOTHER" "$ID_FATHER" "$ID_FETUS1" | sort | tr '\n' ' ')"
ACTUAL="$(bcftools query -l "$QUAD" | tr '\n' ' ')"
[ "$EXPECT" = "$ACTUAL" ] || { echo "❌ quad VCF 样本列顺序 [$ACTUAL] 与预期字母序 [$EXPECT] 不符 —— 基因型会被错位解读。" >&2; exit 1; }
echo "[05q] 样本列顺序核实通过: $ACTUAL"

GT="$WORKDIR/quad_gt.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%AD:%DP:%GQ]\n' "$QUAD" > "$GT"
n=$(wc -l < "$GT"); [ "$n" -gt 1000 ] || { echo "❌ quad_gt.tsv 仅 $n 行，疑静默失效" >&2; exit 1; }
echo "[05q] quad 基因型表: $n 行"

echo "══════ ① 四人交集六模型 ══════"
python3 "$DIR/quad_shared_models.py" "$FULL" "$GT" "$OUT/quad_shared" | tail -5

if [ -f "$PHENO_PANEL" ]; then
  echo "══════ ② 表型优先反向搜索（panel: $(basename "$PHENO_PANEL")，$(grep -cv '^#' "$PHENO_PANEL") 基因）══════"
  python3 "$DIR/pheno_reverse_search.py" "$FULL" "$GT" "$PHENO_PANEL" > "$OUT/pheno_reverse_search.txt"
  grep -c '^  \[' "$OUT/pheno_reverse_search.txt" | xargs echo "  命中条目数:"
else
  echo "⚠️ 未找到表型定向 panel（$PHENO_PANEL），跳过②。表型优先反向搜索是疑难病例的关键步骤，"
  echo "   建议按 HPO 官方 annotation API 生成本例的 panel 后重跑。"
fi

echo "══════ ③ 3-C 夫妇共同携带扫描（不看子代基因型）══════"
python3 "$DIR/carrier_scan_3C.py" "$FULL" "$GT" > "$OUT/carrier_scan_3C.txt"
grep -E '^  → ' "$OUT/carrier_scan_3C.txt" | sed 's/^/  /'

echo
echo "[05q] ✅ 完成。产物: $OUT/quad_shared.tsv, pheno_reverse_search.txt, carrier_scan_3C.txt"
echo "⛔ 下一步：逐条回原始 BAM 过伪影筛（bcftools mpileup 取真实 AD、查均聚物/MAPQ/孟德尔一致性），"
echo "   再进入阶段⑥ ACMG 判读。未过伪影筛的条目不得写进三类变异清单。"
