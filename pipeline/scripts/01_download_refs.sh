#!/usr/bin/env bash
# ============================================================================
# 01_download_refs.sh —— 下载参考基因组 + 索引 + 注释资源
# 占磁盘约 40GB（GRCh38 3GB + bwa-mem2 index ~10GB + VEP cache ~25GB + known-sites）
# ★ 若走外置盘：把 config.sh 的 REFDIR 指过去即可，本脚本全部写入 $REFDIR
# 幂等：已存在则跳过。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
G="$REFDIR/GRCh38"; mkdir -p "$G" "$REFDIR/vep"

need_space() { # 参数: 需要GB；不足则报错退出
  local need=$1 avail; avail=$(df_avail_gb "$REFDIR")
  echo "[refs] $REFDIR 可用 ${avail}GB，本步骤约需 ${need}GB"
  [ "$avail" -lt "$need" ] && { echo "❌ 空间不足，请把 REFDIR 指向外置盘"; exit 1; } || true
}

# 1) 参考基因组（GATK 推荐的 no-alt analysis set）
if [ ! -f "$REF_FASTA" ]; then
  need_space 15
  echo "[refs] 下载 GRCh38 …"
  curl -fSL -o "$REF_FASTA.gz" \
   "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
  gunzip "$REF_FASTA.gz"
fi
# .fai / .dict —— 幂等，独立于上面的下载块（重跑时若缺则补建）
[ -f "$REF_FASTA.fai" ] || samtools faidx "$REF_FASTA"
[ -f "${REF_FASTA%.fna}.dict" ] || gatk CreateSequenceDictionary -R "$REF_FASTA" -O "${REF_FASTA%.fna}.dict"

# 2) bwa-mem2 索引（~10GB，最耗时）
if [ ! -f "$REF_FASTA.bwt.2bit.64" ]; then
  need_space 12; echo "[refs] 构建 bwa-mem2 索引（约20-40分钟）…"
  bwa-mem2 index "$REF_FASTA"
fi

# 3) GATK known-sites（BQSR/VQSR 用）
# 注意：macOS 自带 bash 3.2 不支持关联数组(declare -A)，故用普通数组，保持可移植。
cd "$G"
# 桶 gcp-public-data--broad-references（旧 genomics-public-data 已 403）。用 .gz 版更小、自带 .tbi。
KS_BUCKET="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
KS_URLS="\
$KS_BUCKET/Homo_sapiens_assembly38.dbsnp138.vcf.gz
$KS_BUCKET/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
$KS_BUCKET/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
while IFS= read -r url; do
  [ -z "$url" ] && continue
  f="$G/$(basename "$url")"
  if [ ! -f "$f" ]; then
    echo "[refs] 下载 $(basename "$url")"; curl -fSL -o "$f" "$url"
  fi
  if [ ! -f "$f.tbi" ] && [ ! -f "$f.idx" ]; then
    curl -fsSL -o "$f.tbi" "$url.tbi" 2>/dev/null || gatk IndexFeatureFile -I "$f"
  fi
done <<< "$KS_URLS"

# 4) 外显子 CDS 区间（WES on-target 判定 + tier 兜底用）
if [ ! -f "$G/exons_cds.bed" ]; then
  echo "[refs] 生成 GENCODE CDS 区间…"
  curl -fSL "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz" \
   | gunzip | awk -F'\t' '$3=="CDS"{print $1"\t"$4-11"\t"$5+10}' \
   | sort -k1,1 -k2,2n | bedtools merge -i - > "$G/exons_cds.bed" 2>/dev/null \
   || echo "[refs] (bedtools 缺失可后补；非致命)"
fi

# 5) 注释数据库 —— 已改用 snpEff+SnpSift（VEP 在 arm64 段错误），由 01b_download_snpeff.sh 单独下载。
#    本脚本只负责比对/call 所需的参考（基因组+索引+known-sites+CDS bed）。
echo "[refs] ✅ 参考(比对/call用)完成。"
echo "     注释数据库请另跑: bash 01b_download_snpeff.sh"
echo "     下一步: bash 00b_detect_datatype.sh 判定WES/WGS，再 bash 02_align.sh"
