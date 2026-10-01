#!/usr/bin/env bash
# ============================================================================
# 07_benchmark_giab.sh —— 临床级验证:GIAB HG002/3/4 trio 基准(部署/换机/改流程后一次性)
# HG002(子)/HG003(父)/HG004(母) 正好是公开真值三口之家。用本流水线跑一遍,
# rtg vcfeval(替代 arm64 无包的 hap.py)比对 HG002 与 GIAB 真值,报 SNV/indel 灵敏度/精度,
# 并统计 de novo FDR(HG002 相对 003/004——HG002 真 de novo 极少,标出的多为假阳)。
# ⚠️ 重下载(reads 数十GB)+ 完整跑一遍,耗时耗盘;仅在明确要做临床级验证时运行。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
BM="$RESULTDIR/benchmark"; mkdir -p "$BM"; G="$REFDIR/GRCh38"

# ---- 1) GIAB 真值(小,~几十MB)—— HG002 高置信 small-variant VCF + BED --------
TRUTH_VCF="$G/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="$G/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
GIAB="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
[ -f "$TRUTH_VCF" ] || curl -fSL -o "$TRUTH_VCF" "$GIAB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
[ -f "$TRUTH_VCF.tbi" ] || curl -fSL -o "$TRUTH_VCF.tbi" "$GIAB/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" 2>/dev/null || tabix -p vcf "$TRUTH_VCF"
[ -f "$TRUTH_BED" ] || curl -fSL -o "$TRUTH_BED" "$GIAB/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# ---- 2) reads(★重):HG002/3/4 fastq 或 CRAM ------------------------------
# 数十GB级。两种做法:
#   a) 设 GIAB_HG002_R1/R2、HG003_*、HG004_* 指向已下好的 fastq,把它们当普通病例配进 config 跑 02→05;
#   b) 用 GIAB 提供的比对好 CRAM(省比对)直接从 03 起步。
# 本脚本假定你已用本流水线把 HG002/3/4 跑到 trio.dn.vcf.gz(GATK 联合定型+de novo精修)。
QVCF="${GIAB_QUERY_VCF:-$RESULTDIR/trio.dn.vcf.gz}"    # 待评估的 HG002(子)call 结果
[ -f "$QVCF" ] || { echo "❌ 缺待评估 VCF $QVCF——先用本流水线跑完 GIAB trio(见脚本注释a/b)"; exit 1; }

# ---- 3) rtg vcfeval:SNV/indel 灵敏度/精度 ---------------------------------
SDF="$G/rtg_sdf"
[ -d "$SDF" ] || { echo "[07] 构建 RTG SDF…"; rtg format -o "$SDF" "$REF_FASTA"; }
# 取子样本(HG002)那一列做单样本评估
HG002_COL="${GIAB_HG002_SAMPLE:-HG002}"
rm -rf "$BM/vcfeval"
echo "[07] rtg vcfeval(HG002 vs GIAB v4.2.1)…"
rtg vcfeval -b "$TRUTH_VCF" -c "$QVCF" -e "$TRUTH_BED" -t "$SDF" \
  --sample "$HG002_COL" -o "$BM/vcfeval" --threads "$THREADS" 2>&1 | tail -8 || \
  echo "[07] ⚠️ vcfeval 失败:确认 -c 的样本名(--sample)与 GIAB 真值样本一致"
[ -f "$BM/vcfeval/summary.txt" ] && { echo "=== SNV/indel 灵敏度·精度 ==="; cat "$BM/vcfeval/summary.txt"; }

# ---- 4) de novo FDR:HG002 相对 003/004 ------------------------------------
# HG002 真 de novo ~个位数;流水线标出的 hiConfDeNovo 数 ≫ 真值 → 估 FDR。
NDN=$(bcftools view -H -i 'hiConfDeNovo!="."' "$QVCF" 2>/dev/null | wc -l | tr -d ' ')
echo "[07] hiConf de novo 标记数=$NDN（HG002 真 de novo 约个位数;远超即提示 de novo 假阳偏高,需收紧 03b/过滤）"

echo "[07] ✅ 基准报告 → $BM/  （灵敏度/精度见 vcfeval/summary.txt，de novo FDR 见上）"
