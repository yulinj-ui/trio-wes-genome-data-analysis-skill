#!/usr/bin/env bash
# ============================================================================
# 99_manifest.sh —— 运行溯源清单(每病例跑完后生成)
# 记录:工具版本、参考构建+校验、样本→个体映射、QC汇总(覆盖度/somalier/MCC)、是否降级。
# 便于复现与实验室认可评审。产出 $RESULTDIR/manifest.json。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
OUT="$RESULTDIR/manifest.json"

# 清洗成单行、去引号/制表/回车,避免污染 JSON;$@ 允许传子命令(如 bwa-mem2 version)
ver(){ local v; v=$("$@" 2>&1 | head -1 | tr -d '"\t\r' | sed 's/[[:cntrl:]]//g; s/  */ /g' || true); echo "${v:-NA}"; }
gatk_ver=$(gatk --version 2>&1 | grep -i "Toolkit" | head -1 | tr -d '"\t\r' || echo NA)
snpeff_ver=$(snpEff -version 2>&1 | head -1 | tr -d '"\t\r' || echo NA)
# slivar 版本串里混着 git 报错("fatal: not a git repository")，只取 x.y.z 部分，别让垃圾进溯源文件
slivar_ver=$(slivar 2>&1 | grep -i version | head -1 | sed -n 's/.*version:[[:space:]]*\([0-9.]*\).*/\1/p' || true)
somalier_ver=$(somalier 2>&1 | grep -i version | head -1 | sed -n 's/.*version:[[:space:]]*\([0-9.]*\).*/\1/p' || true)
freemix=$( [ -f "$RESULTDIR/qc/verifybamid_fetus.selfSM" ] && awk 'NR==2{print $7}' "$RESULTDIR/qc/verifybamid_fetus.selfSM" || echo "NA")
dbnsfp_used=$( bcftools view -h "$RESULTDIR/trio.annot.vcf.gz" 2>/dev/null | grep -q "ID=REVEL_score" && echo true || echo false)
spliceai_used=$( bcftools view -h "$RESULTDIR/trio.annot.vcf.gz" 2>/dev/null | grep -q "ID=SpliceAI" && echo true || echo false)
nvar=$(bcftools view -H "$RESULTDIR/trio.annot.vcf.gz" 2>/dev/null | wc -l | tr -d ' ')

cat > "$OUT" <<JSON
{
  "assay": "${ASSAY:-NA}",
  "genome_build": "${GENOME_BUILD:-GRCh38}",
  "reference_fasta": "$REF_FASTA",
  "samples": {
    "father": "$ID_FATHER", "mother": "$ID_MOTHER", "proband": "$ID_FETUS"
  },
  "tools": {
    "bwa_mem2": "$(ver bwa-mem2 version)", "samtools": "$(ver samtools --version)",
    "bcftools": "$(ver bcftools --version)", "gatk": "$gatk_ver", "snpeff": "$snpeff_ver",
    "slivar": "${slivar_ver:-NA}", "somalier": "${somalier_ver:-NA}",
    "verifybamid2": "present", "annovar_dbnsfp": "$DBNSFP_PROTOCOL"
  },
  "run_host": {
    "cpu_cores": ${HOST_NCPU:-0}, "mem_gb": ${HOST_MEM_GB:-0},
    "threads": ${THREADS:-0}, "java_heap_gb": ${MAXMEM_GB:-0}
  },
  "qc": {
    "fetus_FREEMIX": "$freemix", "MCC_THRESHOLD": "$MCC_THRESHOLD"
  },
  "annotation": {
    "dbnsfp_used": $dbnsfp_used, "spliceai_used": $spliceai_used,
    "degraded_note": "$([ "$dbnsfp_used" = false ] && echo 'dbNSFP缺失: PP3/BP4不可用,判读强度下降' || echo '预测注释完整')"
  },
  "variant_count_annotated": $nvar,
  "generated_from": "trio-pipeline scripts 00-07"
}
JSON
echo "[manifest] ✅ $OUT"
cat "$OUT"
