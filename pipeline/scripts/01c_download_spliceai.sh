#!/usr/bin/env bash
# ============================================================================
# 01c_download_spliceai.sh —— v2 增强资产:somalier sites + SpliceAI 预计算 VCF
# somalier sites 很小(~260KB)自动下;SpliceAI 预计算 VCF 大(~29GB)且分发需授权,半自动。
# 一次性,跨病例共享。放 $REFDIR/GRCh38/。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
G="$REFDIR/GRCh38"; mkdir -p "$G"

# ---- 1) somalier sites (GRCh38, chr命名) —— 自动下 --------------------------
if [ ! -f "$SOMALIER_SITES" ]; then
  echo "[01c] 下载 somalier sites (GRCh38)…"
  curl -fSL -o "$SOMALIER_SITES" \
    "https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz"
  tabix -f -p vcf "$SOMALIER_SITES" 2>/dev/null || true
  echo "[01c] ✅ somalier sites → $SOMALIER_SITES"
else
  echo "[01c] somalier sites 已存在，跳过。"
fi

# ---- 2) SpliceAI 预计算 VCF (SNV ~29GB + indel ~1GB) —— 半自动 -------------
# 官方走 Illumina BaseSpace(需登录),无稳定 curl 直链。两种拿法:
#   a) 已有直链/镜像 → 设 SPLICEAI_SNV_URL / SPLICEAI_INDEL_URL 环境变量后重跑本脚本;
#   b) 手动从 BaseSpace(spliceai/genome-annotations)下 spliceai_scores.masked.snv.hg38.vcf.gz
#      与 indel 版,连同 .tbi 放到 $G/ 下(文件名匹配 config.sh 的 SPLICEAI_SNV_VCF/INDEL_VCF)。
# 未就绪也不阻塞:04 会自动跳过 SpliceAI 注释,难点手册的剪接挖掘退化为对小范围现算(见 runbook)。
dl_spliceai() {  # $1=url $2=dest
  [ -z "$1" ] && return 1
  echo "[01c] 下载 $(basename "$2")…"
  curl -fSL -C - --retry 20 --retry-all-errors -o "$2" "$1" || return 1
  [ -f "$2.tbi" ] || tabix -f -p vcf "$2" 2>/dev/null || true
}
if [ ! -f "$SPLICEAI_SNV_VCF" ]; then
  if [ -n "${SPLICEAI_SNV_URL:-}" ]; then
    dl_spliceai "${SPLICEAI_SNV_URL}" "$SPLICEAI_SNV_VCF" && \
    dl_spliceai "${SPLICEAI_INDEL_URL:-}" "$SPLICEAI_INDEL_VCF" || echo "[01c] ⚠️ SpliceAI 下载失败，可重试"
  else
    echo "[01c] ⚠️ 未设 SPLICEAI_SNV_URL —— 跳过 SpliceAI(可选,~29GB)。"
    echo "     拿法: 从 Illumina BaseSpace 'spliceai/genome-annotations' 手动下 GRCh38 的"
    echo "     spliceai_scores.masked.snv.hg38.vcf.gz(+.tbi) 与 indel 版，放到 $G/;"
    echo "     或有直链时: SPLICEAI_SNV_URL='...' SPLICEAI_INDEL_URL='...' bash 01c_download_spliceai.sh"
    echo "     下完 bgzip -t 验完整;未下也能跑,04 会自动跳过 SpliceAI 注释。"
  fi
else
  echo "[01c] SpliceAI SNV VCF 已存在。"
fi
echo "[01c] 完成: somalier-sites=$([ -f "$SOMALIER_SITES" ]&&echo OK), SpliceAI=$([ -f "$SPLICEAI_SNV_VCF" ]&&echo OK||echo 待补)"
