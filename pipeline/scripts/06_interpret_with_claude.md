# 06 · Claude 判读层（排序 + ACMG 致病性判读 + 报告输出）

计算层（01–05）产出 `results/candidates/` 下的候选变异表后，这一层由 **Claude Code + 你已装的技能**完成。**不需要额外生信工具**——这是 Claude 相对纯管线的增量所在。

## 输入
- `results/candidates/tier1_panel_hits.tsv`（表型 panel 内命中，优先）
- `results/candidates/candidates.tsv`（全外显子兜底）
- `scripts/phenotype_HPO.txt`（胎儿表型/HPO）
- `results/trio.ped`（家系+性别）

## 在 Claude Code 里怎么跑（对话即可）

### Step 1 · 表型驱动排序
把候选表贴给 Claude / 让它读取，指令示例：
> "读取 tier1_panel_hits.tsv，按下列维度对候选变异综合打分排序：①与胎儿HPO表型匹配度 ②遗传模型契合度(de novo/AR-comphet/AR-hom/XL) ③gnomAD popmax 频率 ④CADD/REVEL/SpliceAI 危害预测 ⑤ClinVar 记录。输出 Top15 排序表。"

排序权重建议：表型匹配 > 模型契合 > 频率 > 预测。**复发+一代阴性**语境下，额外加权：
- de novo 且基因为 AD-可生殖腺嵌合（如 GREB1L）→ 上调
- AR 复合杂合，其中一个为剪接/深内含子/LoF → 上调（正是上次 exome 易漏的类型）

### Step 2 · 逐变异 ACMG 判读
对 Top 候选，逐个触发致病性评级技能（本项目已设默认规则：见变异即按规格评级）：
> "对 GREB1L NM_xxx c.xxx 走 gene-variant-pathogenicity 评级"

技能会输出 ACMG 证据链（PVS1/PS/PM/PP/BA1…）、数据库交叉校验、功能预测、中文报告。
必要时对标 Franklin：
> "franklin-variant-lookup 查这个变异在 Genoox 的分类"

### Step 3 · 家系共分离与嵌合复核（本例关键）
- 让 Claude 核对 trio 基因型是否与所报遗传模型自洽（相位、父母来源）。
- **生殖腺嵌合排查**：对"胎儿 de novo、父母血样阴性"的强候选，检查父母 CRAM 在该位点是否有极低 VAF 支持 reads（IGV/`samtools mpileup`）——决定复发风险咨询。

### Step 4 · 输出报告
> "用 prenatal-report-interpreter 风格生成中文家系WES报告：候选变异排序表 + Top变异ACMG证据链 + 表型-基因型吻合度论述 + 遗传咨询要点(复发风险/生殖腺嵌合/PGT-M/再次产前诊断建议)。"

## ⚠️ 必须向报告读者交代的分析盲区
上一胎**同表型 + trio WES/CMA 阴性**，标准外显子流程很可能再次漏诊。报告需显式说明本流程对以下类别的覆盖限度，并给出补充建议：
| 变异类别 | 标准 WES 是否覆盖 | 补充手段 |
|---|---|---|
| CNV / 外显子级缺失重复 | 弱 | 结合本次 CMA；或 ExomeDepth/CNVkit |
| 深内含子 / 非经典剪接 | 弱（除非 SpliceAI 放宽内含子） | SpliceAI 全区间、必要时 RNA/WGS |
| 结构变异(SV)/平衡易位 | 否 | WGS |
| 亲代生殖腺嵌合 | 血样测不到 | 低 VAF 复核 + 再发风险咨询 |
| 印记/UPD | 否 | 甲基化/SNP-array |
| 三核苷酸重复扩增 | 否 | 专项检测 |
