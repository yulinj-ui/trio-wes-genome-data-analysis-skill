#!/usr/bin/env bash
# ============================================================================
# 05a_make_allgenes_panel.sh —— 生成"全基因清单"(= 取消 panel 分层)
#
# 当用户选择"全外显子、无 panel 分层"时，05 的 path B(遗传来源标注扫描)仍需要一个
# 基因清单文件。此脚本从【本病例】的注释 VCF 现场生成，避免复用上一病例的清单。
#
# ⚠️ 必须每病例重新生成：panel_ALLGENES.txt 是从某个病例的 annot VCF 抽出来的，
#    跨病例复用会漏掉本病例特有的基因（且不会报错，只会少几条候选）。
#    00a_preflight.sh 会检查该文件是否比本病例的 annot VCF 旧并告警。
#
# 用法: bash 05a_make_allgenes_panel.sh
#       之后: PANEL_FILE=<输出路径> bash 05_inheritance_filter.sh
#       （config.sh 的 PANEL_FILE 已改为 ${PANEL_FILE:-...} 形式，尊重外部覆盖）
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env

ANN="$RESULTDIR/trio.annot.vcf.gz"
[ -f "$ANN" ] || { echo "❌ 需先跑 04_annotate.sh 产出 $ANN" >&2; exit 1; }

OUT="$RESULTDIR/panel_ALLGENES_${CASE_ID:-case}.txt"
# snpEff ANN 字段格式: ALLELE|EFFECT|IMPACT|GENE|GENEID|... → 第 4 列是基因名
bcftools query -f '%INFO/ANN\n' "$ANN" 2>/dev/null | cut -d'|' -f4 | sort -u | grep -v '^$' > "$OUT"

n=$(wc -l < "$OUT" | tr -d ' ')
if [ "$n" -lt 5000 ]; then
  echo "❌ 只抽出 $n 个基因，远低于全外显子应有的量级(约 1.2 万)。" >&2
  echo "   可能是 ANN 字段未注入(查 JDK/04 是否静默失败)或字段分隔位置变化。" >&2
  exit 1
fi
echo "✅ 全基因清单已生成: $OUT （$n 个基因，来源: 本病例 $ANN）"
echo "下一步: PANEL_FILE=\"$OUT\" bash 05_inheritance_filter.sh"
