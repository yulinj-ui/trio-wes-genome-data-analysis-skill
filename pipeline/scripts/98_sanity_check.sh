#!/usr/bin/env bash
# ============================================================================
# 98_sanity_check.sh —— 跑完之后的量级断言（★出报告之前必跑）
#
# 【为什么需要它】
# 本流水线的失败模式不是崩溃，而是【沉默】：退出码 0、文件照样产出、格式完好，
# 但内容是空的、或来自上一个病例。这类错误没有退出码可查，唯一可靠的发现手段
# 就是【核对每一步产物的数量级是否落在合理区间】。
#
# 2026-08-06 某病例的三次实际教训：
#   · on-target 比例算出来是空字符串（awk 取错列）—— 靠"值为空"发现
#   · 亲缘门 somalier relate 抛异常但脚本报"✅ 完成" —— 靠"pairs.tsv 不存在"发现
#   · 05 用了上个病例的 panel，输出 5 条候选（应为约 500 条）—— 靠【数量级】发现
# 第三条最危险：它产出的是一张真实基因、带完整评分、格式无懈可击的表。
# 若不核对数量级，它会被直接当成本病例的候选写进报告。
#
# 用法: bash 98_sanity_check.sh
# 退出码: 0=全部通过  1=存在致命异常（禁止据此出报告）
# ============================================================================
set -uo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env 2>/dev/null || true

ERR=0; WARN=0
fail(){ echo "❌ $*" >&2; ERR=$((ERR+1)); }
warn(){ echo "⚠️  $*" >&2; WARN=$((WARN+1)); }
ok(){   printf "✅ %-46s %s\n" "$1" "${2:-}"; }
# 断言某个数落在 [lo, hi]；超出即视为"疑似静默失效"，而不是"这个病例就是这样"
rng(){ # rng <名称> <实测值> <下限> <上限> <超范围时的提示>
  local name="$1" v="$2" lo="$3" hi="$4" hint="$5"
  if ! [ "${v:-}" -eq "${v:-}" ] 2>/dev/null; then fail "$name: 取不到数值（很可能是静默失效）。$hint"; return; fi
  if [ "$v" -lt "$lo" ] || [ "$v" -gt "$hi" ]; then
    fail "$name = $v，超出合理区间 [$lo, $hi]。$hint"
  else ok "$name" "$v （区间 $lo–$hi）"; fi
}

echo "════════════════════════════════════════════════════════════"
echo " 产物量级断言  CASE_ID=${CASE_ID:-<未设>}  $(date '+%F %H:%M')"
echo "════════════════════════════════════════════════════════════"

VCF="$RESULTDIR/trio.filtered.vcf.gz"; ANN="$RESULTDIR/trio.annot.vcf.gz"
QUADVCF="$RESULTDIR/quad.filtered.vcf.gz"
KING="$RESULTDIR/qc/king_kinship.txt"
FM="${FAMILY_MODE:-trio}"
if [ "$FM" = "quad" ]; then NEXP_PAIR=6; NSAMP=4; else NEXP_PAIR=3; NSAMP=3; fi
echo "  [声明] FAMILY_MODE=$FM （期望 $NSAMP 样本 / $NEXP_PAIR 个配对）"
PED="$RESULTDIR/trio.ped"; PAIRS="$RESULTDIR/qc/somalier/trio.pairs.tsv"
[ -f "$PAIRS" ] || PAIRS="$RESULTDIR/qc/somalier.pairs.tsv"

# ---- 1. 样本身份：VCF 里的样本必须就是 config 里的三个人 -------------------
if [ -f "$VCF" ]; then
  cols=$(bcftools query -l "$VCF" | tr '\n' ' ')
  for s in "$ID_FETUS" "$ID_MOTHER" "$ID_FATHER"; do
    case " $cols " in *" $s "*) :;; *) fail "VCF 中找不到样本 $s。实际样本: [$cols]。config 与结果不是同一个病例！";; esac
  done
  [ "$ERR" = 0 ] && ok "VCF 样本身份与 config 一致" "$cols"
  # 列顺序：05 硬编码 GEN[0]=先证者 GEN[1]=母 GEN[2]=父，错位不报错但结论全反
  first=$(bcftools query -l "$VCF" | head -1)
  [ "$first" = "$ID_FETUS" ] || fail "VCF 第一列是 $first 而非先证者 $ID_FETUS。
     05 脚本硬编码 GEN[0]=先证者，此处错位会让父母基因型张冠李戴且【不报错】。
     修复：改样本 ID 使字母序为 先证者<母<父（如 P1_/P2_/P3_ 前缀），重跑 03。"
else fail "缺 $VCF"; fi

# ---- 2. 变异总数 -----------------------------------------------------------
if [ -f "$VCF" ]; then
  n=$(bcftools view -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
  if [ "${ASSAY:-WES}" = "WES" ]; then
    rng "PASS 变异总数(WES)" "$n" 20000 120000 "WES trio 通常 3-6 万。过低→捕获区间或过滤有问题；过高→可能没限定在 CDS。"
  else
    rng "PASS 变异总数(WGS)" "$n" 3000000 8000000 "WGS trio 通常 400-600 万。"
  fi
fi

# ---- 3. 亲缘门：somalier 初判 + KING 仲裁 ---------------------------------
# ⚠️ 2026-08-24 某四人家系实测：somalier 的 relatedness 会整体偏移（四对亲子 0.72–0.78、
#   无血缘的父母对 0.587），使本断言误报。偏移是【系统性】的（所有配对同向抬高约 +0.25，
#   含亲子对），而真近亲只抬高父母对、亲子对基本不动 —— 形态本身就能区分。
#   故本节改为两段式：somalier 不过 → 自动查 KING-robust 仲裁结果（02c_king_kinship.sh 产出，
#   不依赖等位基因频率面板、对群体结构稳健）。KING 通过则降级为 WARN 并要求报告写明仲裁过程；
#   KING 也不过、或压根没跑 KING → 维持 FAIL。⛔ 任何情况下都不允许手动放行。
king_verdict(){   # 返回 0=KING 确认家系成立  1=KING 判定不成立  2=没有 KING 结果
  [ -f "$KING" ] || return 2
  # 亲子/同胞对的 kinship 应 >0.177（KING 一级亲阈值）；统计一级亲对数
  local n1
  n1=$(awk '$0 ~ /一级亲/ {c++} END{print c+0}' "$KING")
  [ "$n1" -ge "$((NSAMP-1))" ] && return 0 || return 1
}
if [ -f "$PAIRS" ]; then
  bad=$(awk -F'\t' 'NR>1 && $17==0.5 && ($3<0.40 || $3>0.60){c++} END{print c+0}' "$PAIRS")
  npair=$(awk 'NR>1' "$PAIRS" | wc -l | tr -d ' ')
  [ "$npair" -ge 3 ] || fail "somalier pairs 只有 $npair 行（trio 应 3 行 / quad 应 6 行）"
  if [ "$bad" = 0 ]; then
    ok "亲缘门：somalier 所有期望亲子对的 relatedness 在 0.40–0.60"
  else
    king_verdict; kv=$?
    case $kv in
      0) warn "亲缘门：somalier 有 $bad 对亲子落在 0.40–0.60 之外，但 KING-robust 独立重算确认家系成立
     （见 $KING）。判定为 somalier 估计量偏移，非家系问题。
     ⛔ 报告中必须写明【somalier 触发断言失败 → 经 KING 仲裁通过】的完整过程，不得只写\"亲缘门通过\"。" ;;
      1) fail "亲缘门：somalier 有 $bad 对异常，且 KING-robust 仲裁【也未确认】一级亲关系（见 $KING）。
     禁止继续家系推理 —— 高度怀疑样本互换或非亲生。" ;;
      2) fail "亲缘门：somalier 有 $bad 对期望为亲子的关系落在 0.40–0.60 之外，且【未跑 KING 仲裁】。
     ⛔ 不得手动放行。请先跑: bash 02c_king_kinship.sh
     （KING-robust 不依赖等位基因频率面板，是本门的独立仲裁手段）" ;;
    esac
  fi
else
  fail "缺 somalier pairs.tsv → 【亲缘门未通过】。注意 02b 曾有 bug：relate 失败后仍打印'✅ QC 完成'，
     不要仅凭日志里的对勾就认为亲缘门跑过了。"
fi

# ---- 3b. quad 专项产物断言 ------------------------------------------------
if [ "$FM" = "quad" ]; then
  if [ -f "$QUADVCF" ]; then
    qs=$(bcftools query -l "$QUADVCF" 2>/dev/null | wc -l | tr -d ' ')
    [ "$qs" = 4 ] && ok "quad VCF 样本数" "4" \
      || fail "FAMILY_MODE=quad 但 quad.filtered.vcf.gz 只有 $qs 个样本 —— 第二子代未进入联合 call。"
    for s4 in $(all_sample_ids); do
      bcftools query -l "$QUADVCF" 2>/dev/null | grep -qx "$s4" \
        || fail "quad VCF 中找不到样本 $s4"
    done
  else
    fail "FAMILY_MODE=quad 但缺 $QUADVCF —— 四人联合分析未跑（应跑 03q_call_quad.sh 而非 03_call_trio.sh）。"
  fi
  [ -f "$RESULTDIR/quad.ped" ] && ok "quad PED 在位" || fail "缺 quad.ped"
  QS="$RESULTDIR/candidates/quad_shared.tsv"
  if [ -f "$QS" ]; then
    nqs=$(( $(wc -l < "$QS") - 1 ))
    rng "四人交集候选(quad_shared)" "$nqs" 5 3000 "六种交集模型合计。为 0 → 交集脚本静默失效（多为样本列顺序或 join 键错位）。"
  else
    warn "未见 candidates/quad_shared.tsv —— 四人交集分析（复发病例的核心判据）尚未跑，见 quad_shared_models.py"
  fi
fi

[ -f "$PED" ] && ok "PED 在位" || fail "缺 $PED"

# ---- 4. 注释链：每一环都可能"静默注释 0 条" -------------------------------
if [ -f "$ANN" ]; then
  na=$(bcftools view -H "$ANN" 2>/dev/null | grep -c 'ANN='    || echo 0)
  nc=$(bcftools view -H "$ANN" 2>/dev/null | grep -c 'CLNSIG=' || echo 0)
  ns=$(bcftools view -H "$ANN" 2>/dev/null | grep -c 'SpliceAI='|| echo 0)
  [ "$na" -gt 0 ] && ok "snpEff ANN=" "$na 条" || fail "snpEff 注释 0 条 —— 典型无 JDK 症状（wrapper 无 Java 时正常退出、退出码 0）。"
  [ "$nc" -gt 0 ] && ok "ClinVar CLNSIG=" "$nc 条" || fail "ClinVar 注释 0 条。注意 ClinVar VCF 不带 chr 前缀但 SnpSift 能正确匹配，
     不要误以为是染色体命名问题去'修'一个不存在的 bug —— 先查 JDK 与 ClinVar 文件本身。"
  [ "$ns" -gt 0 ] && ok "SpliceAI" "$ns 条" || warn "SpliceAI 注释 0 条。若 SpliceAI VCF 在位却注释 0 条，
     几乎一定是染色体命名不一致（Ensembl 版用 '14'，本流水线用 'chr14'）→ 需 bcftools annotate --rename-chrs。"
  bcftools view -h "$ANN" 2>/dev/null | grep -q "ID=REVEL_score" \
    && ok "dbNSFP 字段已注入" || warn "无 REVEL_score → PP3/BP4 不可用，报告须显式声明降级。"
else fail "缺 $ANN"; fi

# ---- 5. 候选表量级：本次最危险的一类（panel 张冠李戴）---------------------
T1="$RESULTDIR/candidates/tier1_panel_hits.tsv"
CAND="$RESULTDIR/candidates/candidates.tsv"
if [ -f "$T1" ]; then
  nrow=$(( $(wc -l < "$T1") - 1 ))
  ngene=$(awk -F'\t' 'NR>1{print $5}' "$T1" | sort -u | wc -l | tr -d ' ')
  npanel=$(grep -vc '^#' "$PANEL_FILE" 2>/dev/null || echo 0)
  echo "   [参考] PANEL_FILE=$(basename "$PANEL_FILE") 含 $npanel 基因"
  echo "   [声明] PANEL_SCOPE=${PANEL_SCOPE:-<未声明>}  理由: ${PANEL_RATIONALE:-<未填>}"

  # ★ 声明与实际必须一致 —— 这是能确定性抓住"误用上一病例 panel"的唯一判据。
  #   错误 panel 的产物是内部自洽的（候选全在该 panel 内），交叉校验抓不到；
  #   只有把"本例应扫什么范围"声明出来，才有对照物。
  case "${PANEL_SCOPE:-}" in
    allgenes)
      [ "$npanel" -gt 5000 ] \
        && ok "扫描范围声明一致" "allgenes ↔ $npanel 基因" \
        || fail "PANEL_SCOPE 声明为 allgenes（全外显子不分层），实际 PANEL_FILE 只有 $npanel 个基因
     → 几乎一定是误用了其他病例的小 panel。它产出的候选表会与该 panel 完全自洽、
     看不出任何异常，但基因与本病例表型无关。请跑 05a_make_allgenes_panel.sh 后重跑 05。" ;;
    panel)
      [ "$npanel" -le 5000 ] \
        && ok "扫描范围声明一致" "panel ↔ $npanel 基因" \
        || fail "PANEL_SCOPE 声明为 panel，实际却是 $npanel 个基因的全基因清单。" ;;
    *)
      fail "PANEL_SCOPE 未声明（应为 allgenes 或 panel）。换病例时必须显式声明扫描范围，
     否则无法判定 PANEL_FILE 是本病例的还是上一病例的残留。" ;;
  esac

  if [ "$npanel" -gt 5000 ]; then
    # 全外显子无分层：先证者携带的罕见中高危害变异，正常在数百条量级
    rng "候选表行数(全基因扫描)" "$nrow" 100 3000 \
      "全外显子扫描却只有个位/十位数候选，几乎一定是 PANEL_FILE 被静默换成了小 panel。
     实测教训：某病例因 config 无条件赋值覆盖了命令行传参，用了上个病例的 CAKUT+先心 panel，
     只出 5 条候选（应为 496 条），且这 5 条是真实基因、带完整评分、看不出异常。"
  else
    ok "候选表行数(小 panel)" "$nrow 行 / $ngene 基因"
    warn "当前用的是小 panel（$npanel 基因）。请确认这是本病例表型对应的 panel，
     而不是上个病例的残留 —— 小 panel 出少量候选看起来完全正常，不会触发任何告警。"
  fi
  # 候选基因必须是 panel 的子集，否则说明两者根本不是一次运行的产物
  notin=$(awk -F'\t' 'NR>1{print $5}' "$T1" | sort -u | grep -vxF -f <(grep -v '^#' "$PANEL_FILE" | cut -f1) 2>/dev/null | head -5)
  [ -z "$notin" ] && ok "候选基因均在 PANEL_FILE 内" || fail "候选表中出现 PANEL_FILE 之外的基因: $(echo $notin)
     → 候选表与当前 PANEL_FILE 不是同一次运行的产物，请重跑 05。"
else warn "缺 $T1（未跑 05？）"; fi

if [ -f "$CAND" ]; then
  nc2=$(( $(wc -l < "$CAND") - 1 ))
  rng "严格模型候选(de novo/AR/XL)" "$nc2" 1 500 "为 0 说明 slivar 模型全落空，先查 PED 与样本列顺序。"
fi

# ---- 6. de novo 数量：WES trio 的经验区间 ----------------------------------
if [ -f "$ANN" ]; then
  ndn=$(bcftools view -H "$ANN" 2>/dev/null | grep -c 'hiConfDeNovo' || echo 0)
  rng "hiConfDeNovo 数" "$ndn" 1 300 "真实 de novo 编码区通常个位数，其余多为 HLA/黏蛋白/重复区伪影；
     为 0 → 查 PED 与 03b 是否跑过；过高 → 查是否样本关系错配。"
fi

echo "────────────────────────────────────────────────────────────"
if [ "$ERR" -gt 0 ]; then
  echo "❌ 量级断言未通过：$ERR 个致命异常，$WARN 个警告。" >&2
  echo "   ⛔ 在查清之前，禁止把这些产物当作本病例的分析结果写进报告。" >&2
  exit 1
fi
echo "✅ 量级断言通过（$WARN 个警告）。可进入阶段⑤⑥判读。"
exit 0
