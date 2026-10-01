#!/usr/bin/env bash
# ============================================================================
# 00a_preflight.sh —— 开跑前配置守卫（★新病例第一个要跑的脚本）
#
# 【为什么需要它】
# 本流水线历史上的失败，几乎都不是"报错崩掉"，而是"退出码 0、结论却是错的"：
#   · 00b 取错 mosdepth 列 → on-target 恒为空，判定不出 WES/WGS（2026-08-06 实际踩到）
#   · 02b 依赖 03 才生成的 PED → 亲缘门实际没跑，却打印"✅ QC 完成"（同上）
#   · config 无条件覆盖 PANEL_FILE → 05 用了上个病例的 panel，产出一张格式完好、
#     基因却完全无关的候选表（同上，最危险：不核对数量级就发现不了）
# 三者的共同点是【沉默】。本脚本把"换病例时最容易残留上一例配置"的部分，
# 在消耗任何算力之前显式打印出来并逐项断言。
#
# 用法: bash 00a_preflight.sh          # 只检查
#       STRICT=0 bash 00a_preflight.sh # 把致命错降级为警告（不建议）
# ============================================================================
set -uo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
STRICT="${STRICT:-1}"
ERR=0; WARN=0
fail(){ echo "❌ $*" >&2; ERR=$((ERR+1)); }
warn(){ echo "⚠️  $*" >&2; WARN=$((WARN+1)); }
ok(){   echo "✅ $*"; }

echo "════════════════════════════════════════════════════════════════"
echo " 本次病例配置确认（请逐行核对，确认是当前病例而非上一个病例）"
echo "════════════════════════════════════════════════════════════════"
printf "  CASE_ID     : %s\n" "${CASE_ID:-<未设置>}"
printf "  先证者      : %s\n" "$ID_FETUS"
printf "  母          : %s\n" "$ID_MOTHER"
printf "  父          : %s\n" "$ID_FATHER"
[ "${FAMILY_MODE:-trio}" = "quad" ] && printf "  第二子代    : %s  ★FAMILY_MODE=quad\n" "${ID_FETUS1:-<未设>}"
printf "  性别(先/母/父): %s / %s / %s   (1=男 2=女 0=未知)\n" "$SEX_FETUS" "$SEX_MOTHER" "$SEX_FATHER"
printf "  ASSAY       : %s\n" "${ASSAY:-<未定>}"
printf "  DATA_ROOT   : %s\n" "$DATA_ROOT"
printf "  WORKDIR     : %s\n" "$WORKDIR"
printf "  RESULTDIR   : %s\n" "$RESULTDIR"
printf "  PANEL_FILE  : %s\n" "$PANEL_FILE"
printf "  HPO_FILE    : %s\n" "$HPO_FILE"
printf "  本机        : %s核 / %sGB\n" "$HOST_NCPU" "$HOST_MEM_GB"
echo "────────────────────────────────────────────────────────────────"

# ---- 1. 病例隔离：WORK/RESULT 必须落在 CASE_ID 专属子目录 -------------------
# 防止新病例覆盖既往病例结果（refs 是共享的，work/results 不是）。
if [ -z "${CASE_ID:-}" ]; then
  fail "CASE_ID 未设置。新病例必须设 CASE_ID，否则 work/results 会写到根目录并覆盖既往病例。"
else
  case "$WORKDIR"   in *"$CASE_ID"*) ok "WORKDIR 已隔离到病例子目录";;   *) fail "WORKDIR 不含 CASE_ID($CASE_ID)，有覆盖既往病例的风险: $WORKDIR";; esac
  case "$RESULTDIR" in *"$CASE_ID"*) ok "RESULTDIR 已隔离到病例子目录";; *) fail "RESULTDIR 不含 CASE_ID($CASE_ID): $RESULTDIR";; esac
fi

# ---- 2. 残留检测：结果目录里是否已有【别的样本】的产物 ---------------------
# 换病例时忘了改 CASE_ID 的典型症状：RESULTDIR 里躺着上一例的 VCF。
if [ -f "$RESULTDIR/trio.filtered.vcf.gz" ]; then
  if command -v bcftools >/dev/null 2>&1; then
    prev=$(bcftools query -l "$RESULTDIR/trio.filtered.vcf.gz" 2>/dev/null | tr '\n' ',' | sed 's/,$//')
    case ",$prev," in
      *",$ID_FETUS,"*) ok "结果目录中已有本病例的 VCF（续跑场景，正常）";;
      *) fail "结果目录已存在【其他样本】的 VCF：[$prev]，与本次 [$ID_FETUS] 不符。
     很可能是换病例时没改 CASE_ID。请改正或换目录，不要在此基础上续跑。";;
    esac
  fi
fi

# ---- 3. fastq 六件套在位（require_fastq 的前置友好版）---------------------
miss=0
FQ_LIST="$FQ_FETUS_R1 $FQ_FETUS_R2 $FQ_MOTHER_R1 $FQ_MOTHER_R2 $FQ_FATHER_R1 $FQ_FATHER_R2"
NFQ=6; R1_LIST="$FQ_FETUS_R1 $FQ_MOTHER_R1 $FQ_FATHER_R1"
if [ "${FAMILY_MODE:-trio}" = "singleton" ]; then
  FQ_LIST="$FQ_FETUS_R1 $FQ_FETUS_R2"; NFQ=2; R1_LIST="$FQ_FETUS_R1"
fi
if [ "${FAMILY_MODE:-trio}" = "quad" ]; then
  FQ_LIST="$FQ_LIST ${FQ_FETUS1_R1:-} ${FQ_FETUS1_R2:-}"; NFQ=8
  R1_LIST="$R1_LIST ${FQ_FETUS1_R1:-}"
fi
for f in $FQ_LIST; do
  [ -f "$f" ] || { fail "fastq 缺失: $f"; miss=1; }
done
[ "$miss" = 0 ] && ok "${NFQ} 个 fastq 全部在位（FAMILY_MODE=${FAMILY_MODE:-trio}）"

# ★ FAMILY_MODE 意图声明与实际配置的一致性硬断言（同 PANEL_SCOPE 的设计动机）
case "${FAMILY_MODE:-}" in
  singleton)
    # 单人模式：父母槽位必须留空，防止上一个病例的父母 fastq 残留被悄悄带进来
    [ -z "${ID_MOTHER:-}" ] && [ -z "${ID_FATHER:-}" ] \
      || fail "FAMILY_MODE=singleton 但 ID_MOTHER/ID_FATHER 仍有值（${ID_MOTHER:-}/${ID_FATHER:-}）—— 疑为上个病例残留。单人模式必须置空。"
    [ -z "${FQ_MOTHER_R1:-}" ] && [ -z "${FQ_FATHER_R1:-}" ] \
      || fail "FAMILY_MODE=singleton 但 FQ_MOTHER_R1/FQ_FATHER_R1 仍有值 —— 疑为上个病例残留。"
    ok "singleton 模式：仅先证者 ${ID_FETUS} 一人，父母槽位已置空"
    warn "singleton 代价（报告须写明）：de novo 无法本地验证（PS2/PM6 禁赋）、复合杂合无法定相（PM3 只能写\"疑似\"）、亲缘门/MCC 无对照样本。" ;;
  trio)
    [ -z "${ID_FETUS1:-}" ] || warn "FAMILY_MODE=trio 但 ID_FETUS1 已设为 ${ID_FETUS1} —— 疑为上个 quad 病例的残留，第二子代不会被分析。" ;;
  quad)
    [ -n "${ID_FETUS1:-}" ] || fail "FAMILY_MODE=quad 但 ID_FETUS1 未设置。"
    [ -n "${FQ_FETUS1_R1:-}" ] && [ -n "${FQ_FETUS1_R2:-}" ] || fail "FAMILY_MODE=quad 但 FQ_FETUS1_R1/R2 未设置。"
    [ "${SEX_FETUS1:-0}" != "0" ] || warn "FAMILY_MODE=quad 但 SEX_FETUS1 未回填（比对后用 chrX/chrY 覆盖度或 somalier 确定）。"
    # 四个样本 ID 必须互不重复
    ndup=$(printf '%s\n' "$ID_FETUS" "$ID_MOTHER" "$ID_FATHER" "${ID_FETUS1:-}" | sort | uniq -d | wc -l | tr -d ' ')
    [ "$ndup" = 0 ] && ok "四个样本 ID 互不重复" || fail "四个样本 ID 有重复（config 漏改）。" ;;
  *)
    fail "FAMILY_MODE 未声明（应为 singleton / trio / quad）。换病例必须显式声明家系规模 —— 配了四个人却跑三人流程是内部自洽、不会报错的静默失效。" ;;
esac

dup=$(printf '%s\n' $R1_LIST | sort | uniq -d)
[ -n "$dup" ] && fail "三个样本的 R1 路径有重复（config 漏改）: $dup" || ok "各样本 fastq 路径互不重复"

# ---- 4. PANEL_FILE / HPO_FILE：本次最危险的一类 ---------------------------
# 教训：config 若用无条件赋值，命令行传参会被静默覆盖 → 用上个病例的 panel 出一张
#       格式完好但基因无关的候选表。此处强制显式确认。
echo "  [声明] PANEL_SCOPE=${PANEL_SCOPE:-<未声明>}"
echo "  [理由] ${PANEL_RATIONALE:-<未填>}"
if [ ! -s "$PANEL_FILE" ]; then
  fail "PANEL_FILE 不存在或为空: $PANEL_FILE"
else
  ngene=$(grep -vc '^#' "$PANEL_FILE" 2>/dev/null || echo 0)
  ok "PANEL_FILE 可读，含 $ngene 个基因"
  # 声明 vs 实际（与 98_sanity_check.sh 同一判据，此处提前到开跑前）
  case "${PANEL_SCOPE:-}" in
    allgenes) [ "$ngene" -gt 5000 ] || fail "PANEL_SCOPE=allgenes 但 PANEL_FILE 仅 $ngene 基因 —— 疑为其他病例的小 panel 残留。";;
    panel)    [ "$ngene" -le 5000 ] || fail "PANEL_SCOPE=panel 但 PANEL_FILE 有 $ngene 基因（全基因清单）。";;
    *)        fail "PANEL_SCOPE 未声明（应为 allgenes 或 panel）。换病例必须显式声明扫描范围。";;
  esac
  case "${PANEL_RATIONALE:-}" in
    *"<示例>"*) warn "PANEL_RATIONALE 仍是 config.example.sh 的示例占位，换病例时请改写为本例的理由。";;
    ""|"<未填>") warn "PANEL_RATIONALE 为空，建议一句话写明本例选此范围的依据（供报告留痕）。";;
  esac
  case "$(basename "$PANEL_FILE")" in
    panel_ALLGENES.txt)
      # 全基因清单是从【某个病例的】annot VCF 生成的，跨病例复用会漏基因。
      if [ -f "$RESULTDIR/trio.annot.vcf.gz" ] && [ "$PANEL_FILE" -ot "$RESULTDIR/trio.annot.vcf.gz" ]; then
        warn "panel_ALLGENES.txt 比本病例的 trio.annot.vcf.gz 旧 —— 它可能是上个病例生成的。
     跑 05 前请重新生成：bash 05a_make_allgenes_panel.sh"
      fi
      ;;
    panel_CAKUT_CHD.txt)
      warn "当前使用 CAKUT+先心 panel。若本病例不是肾/心表型，这几乎肯定是上个病例的残留 ——
     它不会报错，只会给你一张与表型无关却看起来很正常的候选表。"
      ;;
  esac
fi
[ -s "$HPO_FILE" ] && ok "HPO_FILE 可读" || warn "HPO_FILE 不存在或为空: $HPO_FILE"

# ---- 5. 注释链依赖：JDK（缺它 snpEff/SnpSift 静默注释 0 条）---------------
JAVA_BIN="$JAVA_ENV_BIN/java"
if [ -x "$JAVA_BIN" ]; then
  ok "JDK 在位: $("$JAVA_BIN" -version 2>&1 | head -1)"
elif command -v java >/dev/null 2>&1; then
  ok "JDK 在位(PATH): $(java -version 2>&1 | head -1)"
else
  fail "找不到 JDK。snpEff/SnpSift 的 conda 包只是 python wrapper，无 Java 时它们
     把错误打到 stderr 后【正常退出】，04 退出码 0、VCF 照样产出，但 ANN=/CLNSIG= 全为 0 条。
     修复: conda create -y -n java 'openjdk>=17'
     并把 \$JAVA_ENV_BIN 加到 PATH 最前（该包不建 bin/java 软链）。"
fi

# ---- 6. 参考资产 -----------------------------------------------------------
for f in "$REF_FASTA" "$REFDIR/GRCh38/exons_cds.bed" "$SOMALIER_SITES" \
         "$REFDIR/annot/clinvar_GRCh38.vcf.gz" "$REFDIR/GRCh38/gnomad.hg38.genomes.v3.fix.zip"; do
  [ -e "$f" ] || fail "参考资产缺失: $f"
done
[ -e "$REF_FASTA.bwt.2bit.64" ] || fail "bwa-mem2 索引缺失: $REF_FASTA.bwt.2bit.64"
# 可选增强：缺失只警告，但必须显式说出来（否则报告会以为做了）
[ -e "${SPLICEAI_SNV_VCF:-}" ]   || warn "SpliceAI SNV VCF 缺失 → 剪接注释降级，须写入报告的未覆盖清单"
[ -e "${SPLICEAI_INDEL_VCF:-}" ] || warn "SpliceAI indel VCF 缺失 → indel 剪接影响无法评估，须写入未覆盖清单"
[ -e "$ANNOVAR_HUMANDB/hg38_${DBNSFP_PROTOCOL}.txt" ] || warn "dbNSFP 明文缺失 → 无 REVEL/CADD/AlphaMissense，
     ACMG 的 PP3/BP4 将【不可用】，须在报告中显式声明降级。"

# ---- 7. 内存门槛 -----------------------------------------------------------
if [ "${HOST_MEM_GB:-0}" -ge "${MIN_MEM_GB_COMPUTE:-32}" ]; then
  ok "物理内存 ${HOST_MEM_GB}GB ≥ ${MIN_MEM_GB_COMPUTE}GB，可跑 00b/02/03"
else
  warn "物理内存 ${HOST_MEM_GB}GB < ${MIN_MEM_GB_COMPUTE}GB：比对/call 会被守卫拦下。
     请走双机分工（大内存机跑 00→05，本机做阶段⑤⑥判读）。"
fi

# ---- 8. 磁盘 ---------------------------------------------------------------
avail=$(df_avail_gb "$WORKDIR")
[ -n "${avail:-}" ] && { [ "$avail" -ge 150 ] && ok "工作盘可用 ${avail}GB" || warn "工作盘仅剩 ${avail}GB，WES trio 峰值需 60-90GB"; }

echo "────────────────────────────────────────────────────────────────"
if [ "$ERR" -gt 0 ]; then
  echo "❌ preflight 未通过：$ERR 个致命问题，$WARN 个警告。修正后再跑 00b/02。" >&2
  [ "$STRICT" = "1" ] && exit 1
  echo "⚠️  STRICT=0，强行继续（后果自负）。" >&2
else
  echo "✅ preflight 通过（$WARN 个警告）。下一步: bash 00b_detect_datatype.sh"
fi
exit 0
