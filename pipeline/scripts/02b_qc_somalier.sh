#!/usr/bin/env bash
# ============================================================================
# 02b_qc_somalier.sh —— 一体化 trio QC + 母源细胞污染(MCC)评估
# 产前 trio 头号 QC：MCC 会以中间 VAF 污染变异检出，摧毁 de novo 判定与遗传模型。
# somalier(亲缘/性别/指纹) + VerifyBamID2(FREEMIX 估污染，自带 GRCh38 资源，开箱即用)。
# somalier 只需 CRAM(可 02 后跑)；VerifyBamID2 只需胎儿 CRAM。均 arm64 原生已装。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
QC="$RESULTDIR/qc"; mkdir -p "$QC"; PED="$RESULTDIR/trio.ped"
# ⚠️ 2026-08-05 修（实际病例踩到）：PED 原由 03_call_trio.sh 末尾生成，但 02b 跑在 03 之前，
#    首次跑新病例时 PED 必然不存在。而第 23 行的 ${PED:+...} 判的是"变量非空"（路径字符串恒非空），
#    不是"文件存在"，于是 --ped 照传 → somalier 抛 IOError → 被 `|| true` 吞掉 → 仍打印"✅ QC 完成"。
#    典型的"亲缘门实际没跑，却报告跑过了"。此处改为缺失即自行生成（与 03 的写法保持一致）。
if [ ! -f "$PED" ]; then
  echo "[02b] PED 不存在，按 config 的家系与性别自动生成: $PED"
  { echo -e "#FID\tIID\tPID\tMID\tSEX\tPHENO"
    echo -e "FAM1\t$ID_FATHER\t0\t0\t$SEX_FATHER\t1"
    echo -e "FAM1\t$ID_MOTHER\t0\t0\t$SEX_MOTHER\t1"
    echo -e "FAM1\t$ID_FETUS\t$ID_FATHER\t$ID_MOTHER\t$SEX_FETUS\t2"; } > "$PED"
fi
CRAM_F="$WORKDIR/$ID_FATHER/${ID_FATHER}.final.cram"
CRAM_M="$WORKDIR/$ID_MOTHER/${ID_MOTHER}.final.cram"
CRAM_C="$WORKDIR/$ID_FETUS/${ID_FETUS}.final.cram"

# ---- 1) somalier：亲缘系数 / 遗传性别 / 跨样本指纹（需 SOMALIER_SITES，01c 下载）----
if [ -f "$SOMALIER_SITES" ]; then
  echo "[02b] somalier extract ×3 → relate …"
  rm -rf "$QC/somalier"; mkdir -p "$QC/somalier"
  for cram in "$CRAM_F" "$CRAM_M" "$CRAM_C"; do
    somalier extract -d "$QC/somalier/" --sites "$SOMALIER_SITES" -f "$REF_FASTA" "$cram"
  done
  somalier relate --ped "$PED" -o "$QC/somalier/trio" "$QC/somalier/"*.somalier 2>&1 | tail -2 || true
  # ⛔ 亲缘门是硬门控：relate 没产出结果就不许把本步当"通过"，否则后续家系推理建立在未验证的假设上。
  if [ ! -f "$QC/somalier/trio.pairs.tsv" ]; then
    echo "❌ [02b] somalier relate 未产出 trio.pairs.tsv —— 亲缘门未通过，禁止继续家系推理。" >&2
    exit 1
  fi
  { echo "  亲缘(亲子应≈0.5):"; awk -F'\t' 'NR>1{print "   "$1" <-> "$2": "$3}' "$QC/somalier/trio.pairs.tsv"; }
  [ -f "$QC/somalier/trio.samples.tsv" ] && echo "  推断性别见 $QC/somalier/trio.samples.tsv"
else
  echo "[02b] ⚠️ 无 SOMALIER_SITES（跑 01c 下载它，很小），跳过 somalier；亲缘/性别用手动兜底(见 runbook)"
fi

# ---- 2) MCC —— VerifyBamID2 FREEMIX（主方法，自带 GRCh38 资源，无需下载）-----
# WES 用外显子专属 10k 面板；WGS 用全基因组 100k 面板。资源随 conda 包分发。
VBDIR=$(ls -d "$CONDA_PREFIX/share/"verifybamid2-* 2>/dev/null | head -1)
if [ "${ASSAY:-WES}" = "WES" ]; then
  SVD="${VERIFYBAMID_RES:-$VBDIR/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat}"
else
  SVD="${VERIFYBAMID_RES:-$VBDIR/resource/1000g.phase3.100k.b38.vcf.gz.dat}"
fi
if command -v verifybamid2 >/dev/null && [ -f "$SVD.UD" ]; then
  echo "[02b] VerifyBamID2 估胎儿 FREEMIX（母源污染分数）… 面板: $(basename "$SVD")"
  ( cd "$QC" && verifybamid2 --SVDPrefix "$SVD" --Reference "$REF_FASTA" \
      --BamFile "$CRAM_C" --NumThread "$THREADS" --Output "$QC/verifybamid_fetus" >/dev/null 2>&1 ) || \
      echo "[02b] ⚠️ VerifyBamID2 运行异常，检查 $QC/"
  if [ -f "$QC/verifybamid_fetus.selfSM" ]; then
    FREEMIX=$(awk 'NR==2{print $7}' "$QC/verifybamid_fetus.selfSM")
    echo "[02b] 胎儿 FREEMIX(污染估计) = $FREEMIX （阈值 $MCC_THRESHOLD）"
    awk -v e="$FREEMIX" -v t="$MCC_THRESHOLD" 'BEGIN{
      if(e==""){print "[02b] ⚠️ 未取到 FREEMIX，检查 selfSM"; exit}
      if(e+0 > t+0) printf "[02b] 🚨 WARN: FREEMIX(%.4f) 超阈值(%.4f)！母源污染风险高——胎儿 de novo/低VAF结果不可靠，须重采/重测或在报告显式标注不可靠，不得据此下阴性结论\n", e, t;
      else printf "[02b] ✅ FREEMIX(%.4f) 在阈值(%.4f)内，母源污染风险低\n", e, t }'
  fi
else
  echo "[02b] ⚠️ VerifyBamID2 或其资源不可用，MCC 未评估——产前样本务必设法补此项(重装 verifybamid2 或指定 VERIFYBAMID_RES)"
fi

# 结构化 QC 汇总(供 manifest / 阶段⑥报告)
{ echo -e "metric\tvalue"
  [ -f "$QC/verifybamid_fetus.selfSM" ] && echo -e "fetus_FREEMIX\t$(awk 'NR==2{print $7}' "$QC/verifybamid_fetus.selfSM")"
  echo -e "MCC_THRESHOLD\t$MCC_THRESHOLD"
  [ -f "$QC/somalier/trio.pairs.tsv" ] && awk -F'\t' 'NR>1{print "relatedness_"$1"_"$2"\t"$3}' "$QC/somalier/trio.pairs.tsv"
} > "$QC/qc_summary.tsv"
echo "[02b] ✅ QC 完成 → $QC/qc_summary.tsv"
