#!/usr/bin/env bash
# ============================================================================
# 01b_download_snpeff.sh —— 注释数据库（替代 VEP：snpEff + ClinVar + dbNSFP）
# 三样注释资源：
#   ① snpEff GRCh38 数据库（基因/转录本/后果预测）—— 自动下载，可靠
#   ② ClinVar GRCh38 VCF（临床致病记录）—— NCBI，curl 直下，~200MB
#   ③ dbNSFP（CADD/REVEL/SIFT/PolyPhen/SpliceAI + gnomAD 频率，编码区）—— 大(~30GB)
#      dbNSFP 官方走学术分发链接，需你手动填 URL（见下），故此步半自动。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
AN="$REFDIR/annot"; mkdir -p "$AN"

# ① snpEff 数据库（★如报"database not found"，先跑 snpEff databases|grep GRCh38 选可用名回填）
export SNPEFF_DB="${SNPEFF_DB:-GRCh38.p14}"
SNPEFF_DATADIR="$AN/snpeff_data"; mkdir -p "$SNPEFF_DATADIR"
if [ ! -d "$SNPEFF_DATADIR/$SNPEFF_DB" ]; then
  echo "[annot] 下载 snpEff 数据库 $SNPEFF_DB …"
  snpEff download -dataDir "$SNPEFF_DATADIR" "$SNPEFF_DB" || {
    echo "❌ 数据库名可能不对。可用列表："; snpEff databases | grep -i "GRCh38" | head; exit 1; }
fi

# ② ClinVar GRCh38
CLINVAR="$AN/clinvar_GRCh38.vcf.gz"
if [ ! -f "$CLINVAR" ]; then
  echo "[annot] 下载 ClinVar GRCh38 …"
  curl -fSL -o "$CLINVAR" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
  curl -fSL -o "$CLINVAR.tbi" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
fi

# ③ dbNSFP（半自动）—— 提供功能预测+gnomAD频率。学术版链接需你填：
#    到 https://www.dbnsfp.org/download 取当前 GRCh38 学术版链接，填入 DBNSFP_URL 环境变量后重跑本脚本。
DBNSFP="$AN/dbNSFP.txt.gz"
if [ ! -f "$DBNSFP" ]; then
  if [ -n "${DBNSFP_URL:-}" ]; then
    echo "[annot] 下载 dbNSFP（大，~30GB，数小时）…"
    curl -fSL -o "$AN/dbNSFP.zip" "$DBNSFP_URL"
    (cd "$AN" && unzip -o dbNSFP.zip && \
     zcat dbNSFP*_variant.chr1.gz | head -1 > h.txt && \
     zcat dbNSFP*_variant.chr{1..22} dbNSFP*_variant.chrX dbNSFP*_variant.chrY 2>/dev/null \
       | grep -v "^#chr" | cat h.txt - | bgzip -@ "$THREADS" > "$DBNSFP" && \
     tabix -s 1 -b 2 -e 2 "$DBNSFP" && rm -f dbNSFP.zip dbNSFP*_variant.chr*.gz h.txt)
  else
    echo "⚠️ [annot] 跳过 dbNSFP：未设 DBNSFP_URL。"
    echo "   去 https://www.dbnsfp.org/download 取 GRCh38 学术版链接，然后："
    echo "   DBNSFP_URL='<链接>' bash 01b_download_snpeff.sh"
    echo "   （没有 dbNSFP 也能跑：04 会退化为 snpEff+ClinVar，仅缺 CADD/REVEL/gnomAD 频率注释，"
    echo "     频率过滤可改用 slivar 自带 gnomAD gnotate 库，见 05 注释）"
  fi
fi
# ④ slivar gnomAD gnotate 库（05 遗传模型的频率过滤用；~5GB，单文件 curl 直下，可靠）
GNZIP="$REFDIR/GRCh38/gnomad.hg38.genomes.v3.fix.zip"
if [ ! -f "$GNZIP" ]; then
  echo "[annot] 下载 slivar gnomAD gnotate 库（~5GB）…"
  curl -fSL -o "$GNZIP" \
    "https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip" || \
    echo "⚠️ slivar gnomAD 库下载失败，可稍后重试；05 频率过滤依赖它"
fi

echo "[annot] ✅ 完成: snpEff=$SNPEFF_DB, ClinVar=$([ -f "$CLINVAR" ]&&echo OK), dbNSFP=$([ -f "$DBNSFP" ]&&echo OK||echo 待补), slivar-gnomAD=$([ -f "$GNZIP" ]&&echo OK||echo 待补)"
