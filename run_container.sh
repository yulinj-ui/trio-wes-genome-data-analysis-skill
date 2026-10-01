#!/usr/bin/env bash
# ============================================================================
# run_container.sh —— 在容器里跑流水线的某一步（参考库/数据/结果均以卷挂载）
#
# 用法：
#   REFS=/外置盘/genomics/refs  DATA=/外置盘/seqdata/<病例>/raw \
#   OUT=/外置盘/genomics        CASE_ID=<病例ID> \
#   CONFIG=/path/to/config.sh   ./run_container.sh 02_align.sh father
#
# 挂载约定（容器内路径固定）：
#   $REFS   → /refs      参考库（只读）
#   $DATA   → /data      本病例原始 fastq（只读）
#   $OUT    → /out       其下自动建 work/ results/ logs/ 三个子目录（按 CASE_ID 隔离）
#   $TOOLS  → /tools     ANNOVAR 等需自行注册下载的工具（可选，只读）
#   $CASE   → /case      本病例的 HPO / panel 文件目录（可选，只读）
#   $CONFIG → /pipeline/scripts/config.sh（可选；不给则用 config.example.sh）
# 其余环境变量（THREADS、MAXMEM_GB、FAMILY_MODE、ID_* 等）可在命令前 export 后用
#   EXTRA_ENV="THREADS FAMILY_MODE" 透传。
# ============================================================================
set -euo pipefail

IMAGE="${IMAGE:-trio-wes:latest}"
ENGINE="${ENGINE:-docker}"            # 也可 ENGINE=podman
: "${REFS:?请设 REFS=参考库目录}"
: "${OUT:?请设 OUT=输出根目录（其下建 work/results/logs）}"
: "${CASE_ID:?请设 CASE_ID=病例ID（不要用姓名）}"
[ $# -ge 1 ] || { echo "用法: $0 <脚本名> [参数...]，如 $0 00b_detect_datatype.sh" >&2; exit 2; }

mkdir -p "$OUT/work/$CASE_ID" "$OUT/results/$CASE_ID" "$OUT/logs/$CASE_ID"

args="--rm -i"
[ -t 0 ] && args="$args -t"
# 以当前用户身份写文件，避免结果目录变成 root 所有（Linux 主机上尤其重要）
args="$args --user $(id -u):$(id -g)"
args="$args -v $REFS:/refs:ro -v $OUT:/out"
[ -n "${DATA:-}" ]  && args="$args -v $DATA:/data:ro"
[ -n "${TOOLS:-}" ] && args="$args -v $TOOLS:/tools:ro"
[ -n "${CASE:-}" ]  && args="$args -v $CASE:/case:ro"
# 镜像内置 config.sh = config.example.sh；给了 CONFIG 就覆盖挂载
[ -n "${CONFIG:-}" ] && args="$args -v $CONFIG:/pipeline/scripts/config.sh:ro"

envs="-e CASE_ID=$CASE_ID -e REFDIR=/refs -e TOOLS_DIR=/tools -e SEQ_DATA_ROOT=/data -e DATA_ROOT=/data"
envs="$envs -e WORKDIR=/out/work/$CASE_ID -e RESULTDIR=/out/results/$CASE_ID -e LOGDIR=/out/logs/$CASE_ID"
for v in ${EXTRA_ENV:-}; do
  eval "val=\${$v:-}"
  envs="$envs -e $v=$val"
done

script="$1"; shift
# 注意：路径含空格时请先建不含空格的软链接再传入（docker -v 参数按空格拆分）
# shellcheck disable=SC2086
exec "$ENGINE" run $args $envs "$IMAGE" \
  bash -c "cd /pipeline/scripts && bash \"$script\" \"\$@\"" _ "$@"
