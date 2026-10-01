#!/usr/bin/env bash
# ============================================================================
# 02bq_qc_quad.sh —— 四人版一体化 QC（父+母+两胎）
# 本例特有：两次异常妊娠的胎儿数据都在，亲缘门与 MCC 必须【两个胎儿都做】。
# ⛔ 亲缘门硬门控：relate 未产出结果即 exit 1，不得当"通过"。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
QC="$RESULTDIR/qc"; mkdir -p "$QC"; PED="$RESULTDIR/quad.ped"

# quad PED（父母为 unaffected，两胎均 affected）
{ echo -e "#FID\tIID\tPID\tMID\tSEX\tPHENO"
  echo -e "FAM1\t$ID_FATHER\t0\t0\t$SEX_FATHER\t1"
  echo -e "FAM1\t$ID_MOTHER\t0\t0\t$SEX_MOTHER\t1"
  echo -e "FAM1\t$ID_FETUS\t$ID_FATHER\t$ID_MOTHER\t$SEX_FETUS\t2"
  echo -e "FAM1\t$ID_FETUS1\t$ID_FATHER\t$ID_MOTHER\t$SEX_FETUS1\t2"; } > "$PED"
echo "[02bq] PED: $PED"; cat "$PED"

C_F="$WORKDIR/$ID_FATHER/${ID_FATHER}.final.cram"
C_M="$WORKDIR/$ID_MOTHER/${ID_MOTHER}.final.cram"
C_K2="$WORKDIR/$ID_FETUS/${ID_FETUS}.final.cram"
C_K1="$WORKDIR/$ID_FETUS1/${ID_FETUS1}.final.cram"
for f in "$C_F" "$C_M" "$C_K2" "$C_K1"; do [ -f "$f" ] || { echo "❌ 缺 CRAM: $f"; exit 1; }; done

# ---- 1) somalier：亲缘 / 性别 / 指纹（四样本一起）-------------------------
[ -f "$SOMALIER_SITES" ] || { echo "❌ 缺 SOMALIER_SITES=$SOMALIER_SITES（跑 01c）"; exit 1; }
echo "[02bq] somalier extract ×4 → relate …"
rm -rf "$QC/somalier_quad"; mkdir -p "$QC/somalier_quad"
for cram in "$C_F" "$C_M" "$C_K2" "$C_K1"; do
  somalier extract -d "$QC/somalier_quad/" --sites "$SOMALIER_SITES" -f "$REF_FASTA" "$cram"
done
somalier relate --ped "$PED" -o "$QC/somalier_quad/quad" "$QC/somalier_quad/"*.somalier 2>&1 | tail -3 || true
[ -f "$QC/somalier_quad/quad.pairs.tsv" ] || { echo "❌ [02bq] 亲缘门未通过：无 quad.pairs.tsv" >&2; exit 1; }

echo "──── 亲缘系数（亲子应≈0.5，同胞应≈0.5，无关应≈0）────"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
 {printf "  %-8s <-> %-8s  relatedness=%-8s ibs0=%-7s ibs2=%-7s hom_conc=%s\n",$1,$2,$(h["relatedness"]),$(h["ibs0"]),$(h["ibs2"]),$(h["hom_concordance"])}' \
 "$QC/somalier_quad/quad.pairs.tsv"
echo "──── 遗传性别推断 ────"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
 {printf "  %-8s PED性别=%s  推断性别=%s  X杂合数=%s  X深度=%s  Y深度=%s\n",$2,$(h["original_pedigree_sex"]),$(h["sex"]),$(h["X_het"]),$(h["X_depth_mean"]),$(h["Y_depth_mean"])}' \
 "$QC/somalier_quad/quad.samples.tsv"

# ---- 2) MCC —— VerifyBamID2 对【两个胎儿】各做一次 -------------------------
VBDIR=$(ls -d "$CONDA_PREFIX/share/"verifybamid2-* 2>/dev/null | head -1)
SVD="${VERIFYBAMID_RES:-$VBDIR/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat}"
echo "──── 母源细胞污染 MCC（VerifyBamID2 FREEMIX，阈值 $MCC_THRESHOLD）────"
if command -v verifybamid2 >/dev/null && [ -f "$SVD.UD" ]; then
  for pair in "$ID_FETUS:$C_K2" "$ID_FETUS1:$C_K1"; do
    SID="${pair%%:*}"; CR="${pair#*:}"
    ( cd "$QC" && verifybamid2 --SVDPrefix "$SVD" --Reference "$REF_FASTA" \
        --BamFile "$CR" --NumThread "$THREADS" --Output "$QC/vbid_$SID" >/dev/null 2>&1 ) \
      || echo "  ⚠️ [$SID] VerifyBamID2 运行异常"
    if [ -f "$QC/vbid_$SID.selfSM" ]; then
      FM=$(awk 'NR==2{print $7}' "$QC/vbid_$SID.selfSM")
      FLAG=$(python3 -c "print('⚠️ 超阈值' if float('$FM')>float('$MCC_THRESHOLD') else '✅ 可忽略')")
      echo "  $SID  FREEMIX=$FM  $FLAG"
    else
      echo "  ❌ [$SID] 无 selfSM 产出，MCC 未评估"
    fi
  done
else
  echo "  ❌ verifybamid2 或资源缺失（$SVD.UD），MCC 未评估 —— 不得据此下阴性结论"
fi
echo "[02bq] ✅ 完成，产物在 $QC"
