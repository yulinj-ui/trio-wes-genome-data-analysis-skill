# 家系 WES/WGS 分析流水线镜像（linux/amd64 + linux/arm64）
# 构建：docker buildx build --platform linux/amd64,linux/arm64 -t trio-wes:latest .
# 参考库、原始数据、结果全部通过卷挂载（见 run_container.sh），不进镜像。
FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes

USER root
RUN mkdir -p /refs /data /out /tools /case /pipeline \
 && chown -R $MAMBA_USER:$MAMBA_USER /refs /data /out /tools /case /pipeline \
 && chmod 1777 /out
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER pipeline /pipeline
# 默认配置 = 模板；真实病例用 run_container.sh 的 CONFIG= 挂载覆盖
RUN cp /pipeline/scripts/config.example.sh /pipeline/scripts/config.sh

# 容器内路径约定（config.sh 用 ${VAR:-默认}，这里给出的值即生效）
ENV TRIO_IN_CONTAINER=1 \
    TRIO_PIPELINE_DIR=/pipeline \
    SSD_ROOT=/out \
    REFDIR=/refs \
    TOOLS_DIR=/tools \
    SEQ_DATA_ROOT=/data \
    JAVA_HOME=/opt/conda/lib/jvm \
    JAVA_ENV_BIN=/opt/conda/lib/jvm/bin \
    PATH=/opt/conda/lib/jvm/bin:/opt/conda/bin:$PATH
# ↑ conda 的 openjdk 包不一定建 $PREFIX/bin/java 软链（macOS 上实测没有），故把 lib/jvm/bin 放进 PATH 最前；
#   缺 java 时 snpEff/SnpSift 会以退出码 0 静默注释 0 条。

# 构建期自检：关键工具必须可执行，否则直接构建失败（不把坏镜像推出去）
RUN java -version && samtools --version | head -1 && bcftools --version | head -1 \
 && bwa-mem2 version && gatk --version | head -3 && snpEff -version \
 && mosdepth --version && slivar 2>&1 | grep -i version && somalier 2>&1 | grep -i version

WORKDIR /pipeline/scripts
# micromamba 镜像的 entrypoint 会激活 base 环境；默认打印用法
CMD ["bash", "-c", "ls /pipeline/scripts/*.sh && echo '用法见仓库 README：run_container.sh <脚本> [参数]'"]
