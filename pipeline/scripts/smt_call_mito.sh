#!/usr/bin/env bash
# ============================================================================
# smt_call_mito.sh —— 线粒体专项文库：比对 → Mutect2 线粒体模式 → 异质性变异表
# 适用于厂商单独交付的 mtDNA 富集文库（本例 SEY-Nmtc，实测 chrM 深度 >11000x）。
#
# ⚠️ 已知限度（必须写进报告的「未覆盖清单」）：
#   1. 未使用 GATK 官方的"移位参考(shifted reference)"两遍法 → chrM 控制区跨越
#      线性化断点(约 16024–576)的位点灵敏度下降，该区段结果视为不完整。
#   2. NUMT（核基因组内的线粒体假基因）干扰：本脚本比对到【全基因组】而非单独 chrM，
#      让 NUMT 读段回到它们自己的核位点，这是压 NUMT 假阳性的关键；但低异质性
#      (<1–2%) 位点仍可能是 NUMT 残留，不得据单一低 VAF 位点下结论。
#   3. 大片段 mtDNA 缺失/拷贝数需另做（long-range PCR / 深度断点分析），本脚本不覆盖。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
require_compute_mem "smt_call_mito.sh（bwa-mem2 比对）"
activate_env
_JB="$JAVA_ENV_BIN"; [ -x "$_JB/java" ] && export PATH="$_JB:$PATH"
SID="${ID_MITO:?ID_MITO 未设}"; R1="${FQ_MITO_R1:?}"; R2="${FQ_MITO_R2:?}"
[ -f "$R1" ] && [ -f "$R2" ] || { echo "❌ 线粒体 fastq 不在位: $R1 / $R2" >&2; exit 1; }
O="$WORKDIR/$SID"; mkdir -p "$O" "$RESULTDIR/mito"; LOG="$LOGDIR/smt_mito.log"
exec > >(tee -a "$LOG") 2>&1
echo "===== [$SID] $(date) 线粒体专项分析开始 ====="

CRAM="$O/${SID}.final.cram"
if [ ! -f "$CRAM" ]; then
  fastp -i "$R1" -I "$R2" -o "$O/c1.fq.gz" -O "$O/c2.fq.gz" --detect_adapter_for_pe \
    --thread "$THREADS" --qualified_quality_phred 15 --length_required 50 \
    --json "$O/${SID}.fastp.json" --html "$O/${SID}.fastp.html"
  bwa-mem2 mem -t "$THREADS" -R "@RG\tID:${SID}\tSM:${SID}\tPL:DNBSEQ\tLB:${SID}" \
    "$REF_FASTA" "$O/c1.fq.gz" "$O/c2.fq.gz" \
    | samtools sort -@ "$SORT_THREADS" -m "${SORT_MEM_G}G" -o "$O/${SID}.sorted.bam" -
  samtools index "$O/${SID}.sorted.bam"; rm -f "$O/c1.fq.gz" "$O/c2.fq.gz"
  gatk --java-options "-Xmx${MAXMEM_GB}g" MarkDuplicates -I "$O/${SID}.sorted.bam" \
    -O "$O/${SID}.md.bam" -M "$O/${SID}.dupmetrics.txt" --CREATE_INDEX true
  rm -f "$O/${SID}.sorted.bam" "$O/${SID}.sorted.bam.bai"
  samtools view -@4 -C --output-fmt-option version=3.0 -T "$REF_FASTA" -o "$CRAM" "$O/${SID}.md.bam"
  samtools index "$CRAM"; rm -f "$O/${SID}.md.bam" "$O/${SID}.md.bai"
fi

echo "--- chrM 覆盖度 ---"
samtools depth -a -r chrM "$CRAM" --reference "$REF_FASTA" > "$O/chrM.depth.txt"
awk '{s+=$3; n++; if($3<100) low++} END{printf "chrM 平均深度=%.0f  位点数=%d  深度<100x位点=%d (%.2f%%)\n", s/n, n, low+0, 100*(low+0)/n}' "$O/chrM.depth.txt"
awk '$3<100{print $2}' "$O/chrM.depth.txt" > "$RESULTDIR/mito/chrM_lowcov_positions.txt"

# Mutect2 线粒体模式（高灵敏度、支持低异质性）
MTVCF="$RESULTDIR/mito/${SID}.chrM.vcf.gz"
if [ ! -f "$MTVCF" ]; then
  gatk --java-options "-Xmx${MAXMEM_GB}g" Mutect2 -R "$REF_FASTA" -I "$CRAM" -L chrM \
    --mitochondria-mode -O "$O/${SID}.chrM.raw.vcf.gz"
  gatk --java-options "-Xmx${MAXMEM_GB}g" FilterMutectCalls -R "$REF_FASTA" \
    -V "$O/${SID}.chrM.raw.vcf.gz" --mitochondria-mode -O "$MTVCF"
fi

# 注释：snpEff（chrM 基因）+ ClinVar
# ⛔ 染色体命名坑（2026-08-31 本例实际踩到，一次踩两个注释源）：
#    本流水线参考用 `chrM`，而 snpEff GRCh38.p14 库与 ClinVar VCF 都用 `MT`。
#    snpEff 的自动去 chr 前缀只会把 chrM 变成 M（≠MT），于是【退出码 0、VCF 照产、
#    ANN= 与 CLNSIG= 全部 0 条】—— 与 04_annotate.sh 里 SpliceAI 那个坑同型。
#    故此处先 chrM→MT 再注释，注完改回 chrM，并对两步各加非零断言。
AN="$REFDIR/annot"; export SNPEFF_DB="${SNPEFF_DB:-GRCh38.p14}"; export _JAVA_OPTIONS="-Xmx${MAXMEM_GB}g"
printf 'chrM\tMT\n' > "$O/chrM2MT.txt"; printf 'MT\tchrM\n' > "$O/MT2chrM.txt"
bcftools annotate --rename-chrs "$O/chrM2MT.txt" -Oz -o "$O/mt.MT.vcf.gz" "$MTVCF"; tabix -f -p vcf "$O/mt.MT.vcf.gz"
mt_assert(){ local n; n=$(grep -v '^#' "$1" | grep -c "$2" || true)
  [ "${n:-0}" -gt 0 ] || { echo "❌ [mito] $3 后 $2 命中 0 条 —— 静默失效（多半是染色体命名 chrM vs MT），停。" >&2; exit 1; }
  echo "[mito] ✔ $3: $2 命中 $n 条"; }
snpEff -dataDir "$AN/snpeff_data" -hgvs -noStats "$SNPEFF_DB" "$O/mt.MT.vcf.gz" > "$O/mt.eff.vcf"
mt_assert "$O/mt.eff.vcf" 'ANN=' "snpEff"
SnpSift annotate -info CLNSIG,CLNREVSTAT,CLNDN "$AN/clinvar_GRCh38.vcf.gz" "$O/mt.eff.vcf" > "$O/mt.ann.MT.vcf"
bcftools annotate --rename-chrs "$O/MT2chrM.txt" -o "$O/mt.ann.vcf" "$O/mt.ann.MT.vcf"
SnpSift extractFields -e "." "$O/mt.ann.vcf" CHROM POS REF ALT FILTER \
  "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  CLNSIG CLNREVSTAT CLNDN "GEN[0].AF" "GEN[0].AD" "GEN[0].DP" > "$RESULTDIR/mito/chrM_variants.tsv"
N=$(( $(wc -l < "$RESULTDIR/mito/chrM_variants.tsv") - 1 ))
echo "[mito] ✅ chrM 变异 $N 条 → $RESULTDIR/mito/chrM_variants.tsv"
echo "       AF 列即异质性比例（1.0≈同质）；FILTER=PASS 之外的须人工复核，不得直接丢弃。"
