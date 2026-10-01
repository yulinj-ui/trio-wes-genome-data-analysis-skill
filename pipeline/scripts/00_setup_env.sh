#!/usr/bin/env bash
# ============================================================================
# 00_setup_env.sh —— 安装 miniforge + 创建生信 conda 环境（arm64 原生）
# 幂等：已装则跳过。预计 10-20 分钟，占磁盘 ~3GB。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
mkdir -p "$LOGDIR"

# 容器内（TRIO_IN_CONTAINER=1）工具已装进镜像，直接跳到版本核对。
if [ "${TRIO_IN_CONTAINER:-0}" != "1" ]; then
# 1) miniforge (conda/mamba) —— 用目录/初始化脚本判断，不能靠 PATH（非交互 shell 里 conda 未必在 PATH）
if [ ! -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]; then
  _os=$(uname -s); _arch=$(uname -m)
  [ "$_os" = "Darwin" ] && _os=MacOSX
  echo "[setup] 安装 miniforge (${_os}-${_arch})…"
  curl -fsSL -o /tmp/mf.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${_os}-${_arch}.sh"
  bash /tmp/mf.sh -b -p "$CONDA_ROOT"
  rm -f /tmp/mf.sh
else
  echo "[setup] miniforge 已存在，跳过安装。"
fi
source "$CONDA_ROOT/etc/profile.d/conda.sh"

# 2) 环境 + 工具（bioconda）—— 版本锁定在仓库根目录 environment.yml（与容器镜像同一份）
if ! conda env list | grep -q "^$ENVNAME "; then
  echo "[setup] 创建环境 $ENVNAME …"
  # ★ 不装 ensembl-vep：它在 arm64 段错误,且会拖 perl-bio-samtools 把 samtools 降到古董版 0.1.19。
  #   注释改用 snpeff+snpsift(基因后果/ClinVar)+ ANNOVAR(dbNSFP,单独装,见 04)。
  #   somalier/verifybamid2(QC+MCC,arm64原生自带资源)/rtg-tools(GIAB基准 vcfeval,替代无arm64包的hap.py)。
  conda env create -y -n "$ENVNAME" -f "$DIR/../../environment.yml"
else
  echo "[setup] 环境 $ENVNAME 已存在，跳过。"
fi
fi

activate_env

# ---- 版本核对 --------------------------------------------------------------
# ⚠️ 这些工具取版本的方式各不相同，统一用 `--version` 会产生假的「缺失!」。
#    2026-07-28 的日志里 bwa-mem2 / slivar / vep 三行都是误报，害人误以为环境装坏了。
#    - bwa-mem2 没有 --version，只有 `bwa-mem2 version`
#    - slivar 的版本串里混着 git 报错，必须抓含 version 的那一行
#    - vep 已按设计弃用（arm64 段错误），不该再出现在核对清单里
echo "[setup] 版本核对："
# ⚠️ 每个取版本命令都必须 `|| true`：脚本开了 set -euo pipefail，而 slivar/verifybamid2
#    这类工具无参调用时退出码非 0，会直接中断整个核对循环（只打印到一半就没了）。
vershow(){
  local t="$1" v=""
  printf "  %-12s " "$t"
  command -v "$t" >/dev/null 2>&1 || { echo "❌ 缺失!"; return 0; }
  case "$t" in
    bwa-mem2)     v=$("$t" version 2>&1 | head -1 || true) ;;
    snpEff)       v=$("$t" -version 2>&1 | head -1 || true) ;;
    slivar|somalier)
                  v=$("$t" 2>&1 | grep -i version | head -1 || true) ;;
    gatk)         v=$("$t" --version 2>&1 | grep -i "Toolkit" | head -1 || true) ;;
    *)            v=$("$t" --version 2>&1 | head -1 || true) ;;
  esac
  echo "${v:-已安装(版本串未解析)}"
}
for t in bwa-mem2 samtools bcftools gatk snpEff mosdepth slivar somalier verifybamid2 rtg bedtools fastp; do
  vershow "$t"
done
echo ""
echo "[setup] 本机: ${HOST_NCPU}核 / ${HOST_MEM_GB}GB 内存 → THREADS=$THREADS, MAXMEM_GB=$MAXMEM_GB"
if [ "${HOST_MEM_GB:-0}" -lt "${MIN_MEM_GB_COMPUTE:-32}" ]; then
  echo "[setup] ⚠️ 内存 < ${MIN_MEM_GB_COMPUTE}GB：本机只能做 04/05 与判读，比对(00b/02)和 call(03) 会被守卫拦下。"
  echo "         详见 README.md「硬件门槛与双机分工」。"
fi
echo "[setup] ✅ 完成。下一步: bash 01_download_refs.sh"
