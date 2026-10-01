#!/usr/bin/env bash
# ============================================================================
# 00b_detect_datatype.sh —— 抽样比对判定 WES vs WGS（★存储/流程决策）
# 只取每样本前 200 万读对，比对到 GRCh38，看外显子 on-target 富集度。
# 需先跑完 00 与 01（至少参考基因组+bwa索引就位）。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
require_compute_mem "00b_detect_datatype.sh（抽样比对，同样要加载完整 bwa-mem2 索引）"
require_fastq
activate_env
TMP="$WORKDIR/_detect"; mkdir -p "$TMP"

echo "[detect] 抽样比对母样本前 200 万读对…"
seqkit head -n 2000000 "$FQ_MOTHER_R1" -o "$TMP/r1.fq.gz"
seqkit head -n 2000000 "$FQ_MOTHER_R2" -o "$TMP/r2.fq.gz"
bwa-mem2 mem -t "$THREADS" "$REF_FASTA" "$TMP/r1.fq.gz" "$TMP/r2.fq.gz" 2>/dev/null \
  | samtools sort -@2 -o "$TMP/sub.bam" - ; samtools index "$TMP/sub.bam"

# 外显子区间上的覆盖占比：WES 通常 >60% reads 落在 CDS±100bp；WGS ~2-3%
mosdepth -t "$THREADS" --by "${CAPTURE_BED:-$REFDIR/GRCh38/exons_cds.bed}" \
  -n "$TMP/md" "$TMP/sub.bam" 2>/dev/null || true
# ⚠️ mosdepth summary 列序为: chrom length bases mean min max —— bases 是第 3 列。
#    2026-08-05 修:原先误取 $5(=min,恒为0),导致 ONTGT 恒为空且被 `|| true` 吞掉,
#    表现为"脚本退出码 0、结论是空的"。改用 $3 并加断言。
ONTGT=$(awk '$1=="total_region"{r=$3} $1=="total"{t=$3} END{if(t>0)printf "%.1f",100*r/t}' "$TMP/md.mosdepth.summary.txt" 2>/dev/null || echo "")
if [ -z "$ONTGT" ]; then
  echo "❌ [detect] on-target 比例计算失败(mosdepth 未产出 summary 或列序变化),中间产物保留在 $TMP 供排查" >&2
  exit 1
fi

echo "----------------------------------------------------------------"
echo "[detect] on-target(外显子区)读段占比 ≈ ${ONTGT}%"
echo "   >50%  → WES  （在 config.sh 里 ASSAY=WES，填 CAPTURE_BED）"
echo "   <10%  → WGS  （ASSAY=WGS，CAPTURE_BED 留空；★磁盘务必用外置盘）"
echo "----------------------------------------------------------------"
rm -rf "$TMP"
