#!/usr/bin/env bash
# ============================================================================
# 04_annotate.sh —— 注释（snpEff + ANNOVAR-dbNSFP + SnpSift-ClinVar）
# 链路: snpEff(基因/后果) → ANNOVAR dbnsfp(CADD/REVEL/AlphaMissense/MetaRNN…) → SnpSift annotate(ClinVar)
# dbNSFP(ANNOVAR格式)缺失/未解压时自动降级为 snpEff+ClinVar（频率过滤改由 05 的 slivar gnotate 承担）。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
export _JAVA_OPTIONS="-Xmx${MAXMEM_GB}g"   # snpEff/SnpSift 默认堆太小，加载 GRCh38.p14 会 OOM

# ---- 0) Java 预检（★静默失效防护，2026-07-30 某病例实际踩到）------
# snpEff / SnpSift 的 bioconda 包只是 python wrapper：环境里没有 JDK 时，它们把
# "Unable to locate a Java Runtime" 打到 stderr 后**以退出码 0 正常退出**。
# 后果：本脚本退出码 0、VCF 文件照常产出，但 ANN=/CLNSIG=/SpliceAI= 全部 0 条，
# 下游 05 与阶段⑥判读完全无从察觉。故此处硬性拦下，绝不带着无注释的 VCF 往下走。
# 独立 java 环境若在标准位置，自动补进 PATH（每次跑不必手动 export）。
_JAVA_BIN="$JAVA_ENV_BIN"
[ -x "$_JAVA_BIN/java" ] && case ":$PATH:" in *":$_JAVA_BIN:"*) ;; *) export PATH="$_JAVA_BIN:$PATH";; esac
if ! java -version >/dev/null 2>&1; then
  cat >&2 <<'EOF'

❌ [annot] 本环境找不到可用的 Java 运行时（java -version 失败），拒绝执行注释。
   原因：snpEff / SnpSift 的 bioconda 包只是 python wrapper —— 无 JDK 时它们会把
   "Unable to locate a Java Runtime" 打到 stderr 后**以退出码 0 正常退出**，
   于是脚本"成功"产出 VCF，而 ANN= / CLNSIG= / SpliceAI= 一条都没有（静默失效）。
   （2026-07-30 某病例实际踩到：04 退出码 0，注释条数全 0。）

   修复 —— 装到独立环境，然后把 lib/jvm/bin 加进 PATH：
     conda create -y -n java "openjdk>=17"
     export PATH="$JAVA_ENV_BIN:$PATH"
     java -version        # 确认能打印版本号
   ⚠️ conda 的 openjdk 包**不会**创建 $PREFIX/bin/java 符号链接 ——
      PATH 必须指向 envs/java/lib/jvm/bin，只加 envs/java/bin 无效（这是坑点）。
      本脚本会自动探测上述标准路径，故按此命令装完即可直接重跑 04。
   ⚠️ **不要**往 trio 环境里 conda install openjdk：该环境有过 perl-bio-samtools
      把 samtools 连带降级成 0.1.19 的事故史，动它的依赖解算风险太高。

EOF
  exit 1
fi
echo "[annot] java: $(java -version 2>&1 | head -1)"

# ---- 0b) 注释条数断言（同上，防一切"退出码 0 但没注上"的静默失效）---------
# 只数数据行（header 里的 ##INFO=<ID=ANN,...> 不算命中）。
count_tag(){ grep -v '^#' "$1" 2>/dev/null | grep -c "$2" || true; }
assert_annotated(){   # $1=vcf $2=INFO标签正则 $3=步骤名
  local n; n=$(count_tag "$1" "$2")
  if [ "${n:-0}" -eq 0 ]; then
    if [ "${ALLOW_EMPTY_ANNOT:-0}" = "1" ]; then
      echo "⚠️  [annot] ${3}: ${2} 命中 0 条，但 ALLOW_EMPTY_ANNOT=1 已设置，继续（后果自负）。" >&2
      return 0
    fi
    {
      echo ""
      echo "❌ [annot] ${3} 完成后 ${2} 命中 0 条 —— 该步实际未生效，拒绝把无注释的 VCF 交给 05。"
      echo "   最常见根因：环境内无 JDK（snpEff/SnpSift wrapper 静默退出，见本脚本开头预检）。"
      echo "   其次：注释源缺失/未建索引，或注释源与 VCF 染色体命名不一致（chr14 vs 14）。"
      echo "   排查：直接看上一步的 stderr；并确认 $1 里确有变异行（bcftools view -H | head）。"
      echo "   若确认确为真阴性（如极小 panel VCF 无 ClinVar 命中）：ALLOW_EMPTY_ANNOT=1 bash 04_annotate.sh"
      echo ""
    } >&2
    exit 1
  fi
  echo "[annot] ✔ ${3}: ${2} 命中 ${n} 条"
}
AN="$REFDIR/annot"; export SNPEFF_DB="${SNPEFF_DB:-GRCh38.p14}"
# 输入优先用 03b 的 de novo 精修产物(带 hiConf/loConfDeNovo)，无则退回 03 的过滤产物
# ★ VCF_PREFIX 由 config.sh 按 FAMILY_MODE 给出（trio/quad=trio，singleton=proband）。
VP="${VCF_PREFIX:-trio}"
IN="$RESULTDIR/${VP}.filtered.vcf.gz"; [ -f "$RESULTDIR/${VP}.dn.vcf.gz" ] && IN="$RESULTDIR/${VP}.dn.vcf.gz"
OUT="$RESULTDIR/${VP}.annot.vcf.gz"; echo "[annot] 输入: $IN"
[ -f "$IN" ] || { echo "❌ [annot] 输入 VCF 不存在: $IN（FAMILY_MODE=${FAMILY_MODE:-trio}）" >&2; exit 1; }
CLINVAR="$AN/clinvar_GRCh38.vcf.gz"
# ANNOVAR dbNSFP：需要解压后的 hg38_${DBNSFP_PROTOCOL}.txt + .idx（见 config.sh 5b 段）
DBNSFP_TXT="$ANNOVAR_HUMANDB/${ANNOVAR_BUILD}_${DBNSFP_PROTOCOL}.txt"
tmp1="$WORKDIR/ann.snpeff.vcf"; tmp2="$WORKDIR/ann.dbnsfp.vcf"

# 1) snpEff：基因/转录本/变异后果（HIGH/MODERATE… + HGVS），写 ANN 字段
echo "[annot] snpEff ($SNPEFF_DB) …"
snpEff -dataDir "$AN/snpeff_data" -hgvs -canon -noStats \
  "$SNPEFF_DB" "$IN" > "$tmp1"
assert_annotated "$tmp1" 'ANN=' "snpEff"          # 硬要求：无基因/后果注释则 05 无从筛选
N_ANN=$(count_tag "$tmp1" 'ANN=')

# 2) dbNSFP：功能预测 + gnomAD 频率（ANNOVAR table_annovar；需解压后的 .txt 与 .idx 齐备）
if [ -f "$DBNSFP_TXT" ] && [ -f "$DBNSFP_TXT.idx" ] && [ -x "$ANNOVAR_DIR/table_annovar.pl" ]; then
  echo "[annot] ANNOVAR dbNSFP ($DBNSFP_PROTOCOL) …"
  perl "$ANNOVAR_DIR/table_annovar.pl" "$tmp1" "$ANNOVAR_HUMANDB" \
    -buildver "$ANNOVAR_BUILD" -out "$WORKDIR/annovar_dbnsfp" \
    -protocol "$DBNSFP_PROTOCOL" -operation f -nastring . -vcfinput -remove -thread "$THREADS"
  # ANNOVAR -vcfinput 产出 {out}.${ANNOVAR_BUILD}_multianno.vcf，dbNSFP列已写入INFO并保留原ANN
  mv "$WORKDIR/annovar_dbnsfp.${ANNOVAR_BUILD}_multianno.vcf" "$tmp2"
  rm -f "$WORKDIR/annovar_dbnsfp"*.avinput "$WORKDIR/annovar_dbnsfp"*.txt 2>/dev/null || true
  # 可选增强：为 0 不退出，但必须醒目告警 + 在 manifest 标降级（PP3/BP4 不可用）
  N_DBNSFP=$(count_tag "$tmp2" 'REVEL_score=[^.;[:space:]]')
  if [ "${N_DBNSFP:-0}" -eq 0 ]; then
    echo "⚠️⚠️ [annot] dbNSFP 跑完但 REVEL_score 实值 0 条 —— 功能预测**未实际注上**（降级运行）。" >&2
    echo "        后果：ACMG 的 PP3/BP4 不可用，判读强度下降；已在 manifest 标注降级。" >&2
    echo "        排查：hg38_${DBNSFP_PROTOCOL}.txt/.idx 是否完整、染色体命名是否一致。" >&2
  else
    echo "[annot] ✔ dbNSFP: REVEL_score 实值 ${N_DBNSFP} 条"
  fi
else
  echo "⚠️⚠️ [annot] 无可用 ANNOVAR dbNSFP（.txt/.idx 未就绪），跳过功能预测 → **降级运行**" >&2
  echo "        PP3/BP4 不可用；频率过滤改由 05 的 slivar gnotate 承担。已在 manifest 标注。" >&2
  cp "$tmp1" "$tmp2"
  N_DBNSFP=0
fi

# 3) ClinVar：临床意义
tmp3="$WORKDIR/ann.clinvar.vcf"
if [ ! -f "$CLINVAR" ]; then
  # ClinVar 不是可选增强：缺它等于判读没有临床意义证据，不允许静默跳过
  echo "❌ [annot] ClinVar VCF 不存在: $CLINVAR" >&2
  echo "   跑 01_download_refs.sh 补齐后重跑 04（缺 ClinVar 的注释结果不可用于判读）。" >&2
  exit 1
fi
echo "[annot] SnpSift annotate ClinVar …"
SnpSift annotate -info CLNSIG,CLNREVSTAT,CLNDN "$CLINVAR" "$tmp2" > "$tmp3"
assert_annotated "$tmp3" 'CLNSIG=' "ClinVar"      # 硬要求：同上，0 条即 SnpSift 没真跑
N_CLNSIG=$(count_tag "$tmp3" 'CLNSIG=')

# 4) SpliceAI：剪接影响预测（有哪个预计算 VCF 就注哪个；SNV/indel 各自可选，二者皆无才跳过）
#    注意 SpliceAI VCF 染色体命名须与参考一致(带 chr)——Ensembl 版需先 bcftools --rename-chrs 转换。
if [ -f "$SPLICEAI_SNV_VCF" ] || [ -f "$SPLICEAI_INDEL_VCF" ]; then
  echo "[annot] SpliceAI 注释（SNV=$([ -f "$SPLICEAI_SNV_VCF" ]&&echo 有||echo 无) indel=$([ -f "$SPLICEAI_INDEL_VCF" ]&&echo 有||echo 无)）…"
  cp "$tmp3" "$WORKDIR/ann.sp.vcf"
  if [ -f "$SPLICEAI_SNV_VCF" ]; then
    SnpSift annotate -info SpliceAI "$SPLICEAI_SNV_VCF" "$WORKDIR/ann.sp.vcf" > "$WORKDIR/ann.sp2.vcf" && mv "$WORKDIR/ann.sp2.vcf" "$WORKDIR/ann.sp.vcf"
  fi
  if [ -f "$SPLICEAI_INDEL_VCF" ]; then
    SnpSift annotate -info SpliceAI "$SPLICEAI_INDEL_VCF" "$WORKDIR/ann.sp.vcf" > "$WORKDIR/ann.sp2.vcf" && mv "$WORKDIR/ann.sp2.vcf" "$WORKDIR/ann.sp.vcf"
  fi
  # 可选增强：为 0 不退出，但必须醒目告警 + 在 manifest 标降级
  N_SPLICEAI=$(count_tag "$WORKDIR/ann.sp.vcf" 'SpliceAI=')
  if [ "${N_SPLICEAI:-0}" -eq 0 ]; then
    echo "⚠️⚠️ [annot] SpliceAI 跑完但 SpliceAI= 0 条 —— 剪接注释**未实际注上**（降级运行）。" >&2
    echo "        典型根因：注释源染色体命名不一致（Ensembl 版 '14' vs 流水线 'chr14'）；" >&2
    echo "        或无 JDK 致 SnpSift 静默退出。验证：tabix \$SPLICEAI_SNV_VCF chr14:53951900-53952000" >&2
  else
    echo "[annot] ✔ SpliceAI: SpliceAI= 命中 ${N_SPLICEAI} 条"
  fi
  bgzip -@ "$THREADS" -c "$WORKDIR/ann.sp.vcf" > "$OUT"; rm -f "$WORKDIR/ann.sp.vcf"
else
  echo "⚠️⚠️ [annot] 无 SpliceAI 预计算 VCF（跑 01c 下载），跳过剪接注释 → **降级运行**" >&2
  echo "        深内含子/剪接挖掘将无 SpliceAI 分。已在 manifest 标注。" >&2
  bgzip -@ "$THREADS" -c "$tmp3" > "$OUT"
  N_SPLICEAI=0
fi
tabix -p vcf "$OUT"; rm -f "$tmp1" "$tmp2" "$tmp3"

# ---- 注释条数落盘，供 99_manifest.sh 写入 manifest.json（事后追溯是否降级运行）----
_DEG=""; [ "${N_DBNSFP:-0}" -eq 0 ] && _DEG="dbnsfp"
[ "${N_SPLICEAI:-0}" -eq 0 ] && _DEG="${_DEG:+$_DEG,}spliceai"
cat > "$RESULTDIR/annot_counts.env" <<ENV
N_ANN=${N_ANN:-0}
N_CLNSIG=${N_CLNSIG:-0}
N_DBNSFP=${N_DBNSFP:-0}
N_SPLICEAI=${N_SPLICEAI:-0}
ANNOT_DEGRADED="${_DEG}"
ANNOT_JAVA="$(java -version 2>&1 | head -1 | tr -d '"')"
ENV
echo "[annot] 注释条数: ANN=${N_ANN:-0} CLNSIG=${N_CLNSIG:-0} dbNSFP=${N_DBNSFP:-0} SpliceAI=${N_SPLICEAI:-0}${_DEG:+  ⚠️降级: $_DEG}"
echo "[annot] ✅ 注释完成: $OUT"
echo "下一步: bash 05_inheritance_filter.sh"
