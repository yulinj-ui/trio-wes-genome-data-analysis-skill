#!/usr/bin/env bash
# ============================================================================
# 05_inheritance_filter.sh —— 遗传模型筛选（AD/de novo · AR纯合 · AR复合杂合 · XL）
# 两条互补路径：
#   A) slivar expr：genome-wide 标注 denovo/arhom/xlinked（基因型+gnomAD频率规则）
#   B) panel 内自实现扫描：胎儿携带的罕见中高危害变异，按遗传来源分类
#      （补 A 的盲区——外显不全的显性遗传：变异来自"表面正常"的父/母一方，
#       严格 de novo/AR 模型不会标记，但对本例"复发+一代阴性+疑嵌合/低外显"场景很关键）
# 输出交给 06 判读层做 ACMG 分级。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
activate_env
export _JAVA_OPTIONS="-Xmx${MAXMEM_GB}g"
VCF="$RESULTDIR/trio.annot.vcf.gz"; PED="$RESULTDIR/trio.ped"
GNOMAD_ZIP="$REFDIR/GRCh38/gnomad.hg38.genomes.v3.fix.zip"
OUTDIR="$RESULTDIR/candidates"; mkdir -p "$OUTDIR"

# ---- A) slivar expr: genome-wide denovo/arhom/xlinked 标注 ----------------
DN="INFO.gnomad_popmax_af < $MAX_AF && kid.hom_ref==false && kid.GQ>=20 && kid.AB>0.25 && mom.hom_ref && dad.hom_ref && mom.GQ>=20 && dad.GQ>=20"
AR="INFO.gnomad_popmax_af < 0.005 && kid.hom_alt && kid.GQ>=20 && mom.het && dad.het"
XL="INFO.gnomad_popmax_af < $MAX_AF && kid.alts>0 && kid.GQ>=20"
GNARG=""; [ -f "$GNOMAD_ZIP" ] && GNARG="-g $GNOMAD_ZIP"

echo "[filter] slivar expr —— genome-wide de novo / AR-hom(slivar版,多等位站点需人工复核) / XL"
slivar expr --vcf "$VCF" --ped "$PED" $GNARG \
  --info "INFO.impactful" \
  --trio "denovo:$DN" \
  --trio "arhom:$AR" \
  --trio "xlinked:$XL && variant.CHROM=='chrX'" \
  --out-vcf "$WORKDIR/models.vcf" 2>>"$LOGDIR/05.log"

# ---- 导出全量变异表（genome-wide，含 slivar 标签），供两条路径共用 ---------
# ★ 动态探测 04 是否注入 dbNSFP 字段（取决于下载是否就绪），缺失不请求，避免 SnpSift 报字段不存在崩溃。
# ★ VCF 样本列由 GATK 按名字母序排列 = Fetus_010, GaoHong_011(母), WangQiang_012(父)
#   即 GEN[0]=胎儿 GEN[1]=母 GEN[2]=父 —— 与直觉的"父母胎"顺序不同，务必对应正确。
HDR_FIELDS=$(bcftools view -h "$WORKDIR/models.vcf" 2>/dev/null | grep "^##INFO" | sed -E 's/^##INFO=<ID=([^,]+).*/\1/')
BASE_FIELDS='CHROM POS REF ALT "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" LOF NMD'
OPT_FIELDS=""; OPT_HDR=""
# 动态探测可选注释字段:dbNSFP 预测分(ANNOVAR 裸列名/SnpSift dbNSFP_ 前缀两兼容)、
# 03b 的 de novo 高/低置信标签、SpliceAI 剪接分。只把 VCF 头里真实存在的加入,缺失自动跳过。
# SpliceAI 字段值为 ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|...；DS_max=第3-6项(四个delta)取最大,判读时解析。
for f in CADD_phred REVEL_score AlphaMissense_pred AlphaMissense_score MetaRNN_pred SIFT_pred Polyphen2_HDIV_pred \
         dbNSFP_CADD_phred dbNSFP_REVEL_score dbNSFP_AlphaMissense_pred dbNSFP_MetaRNN_pred \
         hiConfDeNovo loConfDeNovo SpliceAI; do
  if echo "$HDR_FIELDS" | grep -qx "$f"; then OPT_FIELDS="$OPT_FIELDS $f"; OPT_HDR="$OPT_HDR\t${f#dbNSFP_}"; fi
done
[ -z "$OPT_FIELDS" ] && echo "[filter] ⚠️ dbNSFP 尚未就绪，候选表暂缺 CADD/REVEL 等预测分；dbNSFP 下载完成后重跑 04+05 可补全"
TAIL_FIELDS='CLNSIG gnomad_popmax_af "GEN[0].GT" "GEN[1].GT" "GEN[2].GT" denovo arhom xlinked FILTER'

echo "[filter] 导出全量变异表…"
eval SnpSift extractFields -e '"."' "\"\$WORKDIR/models.vcf\"" $BASE_FIELDS $OPT_FIELDS $TAIL_FIELDS 2>/dev/null \
  | tail -n +1 > "$WORKDIR/all_variants_full.tsv"
NCOL_OPT=$(echo "$OPT_FIELDS" | wc -w | tr -d ' ')

# candidates.tsv：genome-wide，仅 slivar 标注出 denovo/arhom/xlinked 的行（严格模型）
# TAIL_FIELDS 末9列固定顺序: CLNSIG,gnomad_popmax_af,GT_fetus,GT_mother,GT_father,denovo,arhom,xlinked,FILTER
# 即 denovo=$(NF-3) arhom=$(NF-2) xlinked=$(NF-1) FILTER=$NF；打印到 GT_father(=$(NF-4)) 为止。
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==1{next}
  { dn=$(NF-3); ah=$(NF-2); xl=$(NF-1);
    if ((dn!="."&&dn!="") || (ah!="."&&ah!="") || (xl!="."&&xl!="")) {
      m=((dn!="."&&dn!="")?"deNovo;":"")((ah!="."&&ah!="")?"AR-hom(slivar,核实多等位);":"")((xl!="."&&xl!="")?"XL;":"");
      s=$1; for(i=2;i<=NF-4;i++) s=s"\t"$i; print s"\t"m
    }
  }' "$WORKDIR/all_variants_full.tsv" > "$OUTDIR/_body.tsv"
HDR="CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tHGVS_c\tHGVS_p\tLOF\tNMD${OPT_HDR}\tCLNSIG\tgnomAD_popmax\tGT_fetus\tGT_mother\tGT_father\tMODEL"
{ echo -e "$HDR"; cat "$OUTDIR/_body.tsv"; } > "$OUTDIR/candidates.tsv"; rm -f "$OUTDIR/_body.tsv"

# ---- B) panel 内自实现扫描：AR-hom(多等位安全) + AR复合杂合(反式) + 遗传来源标注 ----
echo "[filter] panel基因内扫描(胎儿携带+罕见+中高危害, 标注遗传来源含外显不全显性候选)…"
grep -v '^#' "$PANEL_FILE" | cut -f1 | grep -v '^$' > "$WORKDIR/panel_genes.txt"

python3 - "$WORKDIR/all_variants_full.tsv" "$WORKDIR/panel_genes.txt" "$OUTDIR/tier1_panel_hits.tsv" "$NCOL_OPT" <<'PYEOF'
import csv, sys
from collections import defaultdict

full_tsv, panel_file, out_file, nopt = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
panel_genes = set(l.strip() for l in open(panel_file) if l.strip())

def gt_alleles(gt):
    if gt in ('.', './.', '.|.', ''): return None
    gt = gt.replace('|', '/')
    try: return [int(p) for p in gt.split('/')]
    except: return None

def gt_class(gt):
    a = gt_alleles(gt)
    if a is None: return 'missing'
    if len(set(a)) == 1: return 'hom_ref' if a[0] == 0 else 'hom_alt'
    return 'het'

with open(full_tsv) as f:
    header = f.readline().rstrip('\n').split('\t')
rows = []
with open(full_tsv) as f:
    r = csv.DictReader(f, delimiter='\t', fieldnames=header)
    next(r)
    for row in r:
        if row['ANN[0].GENE'] not in panel_genes: continue
        if gt_class(row['GEN[0].GT']) not in ('het', 'hom_alt'): continue   # 胎儿必须携带
        try: af = float(row['gnomad_popmax_af'])
        except: af = 0.0   # gnotate '.'=不在gnomAD=按极罕见处理
        if af >= 0.01: continue
        if row['ANN[0].IMPACT'] not in ('HIGH', 'MODERATE'): continue
        rows.append(row)

by_gene = defaultdict(list)
for row in rows: by_gene[row['ANN[0].GENE']].append(row)

results = []
for gene, vs in by_gene.items():
    het_variants = [v for v in vs if gt_class(v['GEN[0].GT']) == 'het']
    comphet_partner = {}
    for i in range(len(het_variants)):
        for j in range(i + 1, len(het_variants)):
            a, b = het_variants[i], het_variants[j]
            a_dad, a_mom = gt_class(a['GEN[2].GT']), gt_class(a['GEN[1].GT'])
            b_dad, b_mom = gt_class(b['GEN[2].GT']), gt_class(b['GEN[1].GT'])
            if (a_dad == 'het' and a_mom == 'hom_ref' and b_mom == 'het' and b_dad == 'hom_ref') or \
               (a_mom == 'het' and a_dad == 'hom_ref' and b_dad == 'het' and b_mom == 'hom_ref'):
                comphet_partner[a['POS']] = b['POS']; comphet_partner[b['POS']] = a['POS']
    for row in vs:
        fetus, mom, dad = gt_class(row['GEN[0].GT']), gt_class(row['GEN[1].GT']), gt_class(row['GEN[2].GT'])
        if row['POS'] in comphet_partner:
            pattern = f"AR复合杂合(配对POS={comphet_partner[row['POS']]})"
        elif fetus == 'hom_alt' and mom == 'het' and dad == 'het':
            pattern = "AR纯合(父母各携带一份)"
        elif fetus == 'hom_alt' and (mom == 'hom_ref' or dad == 'hom_ref'):
            pattern = "hom_alt但父母基因型不支持典型AR-需人工复核(可能denovo/UPD/测序误差)"
        elif mom == 'hom_ref' and dad == 'hom_ref':
            pattern = "denovo(父母均hom_ref)"
        elif dad == 'het' and mom == 'hom_ref':
            pattern = "遗传自父(单杂合,若为AD基因=外显不全候选)"
        elif mom == 'het' and dad == 'hom_ref':
            pattern = "遗传自母(单杂合,若为AD基因=外显不全候选)"
        elif mom == 'het' and dad == 'het':
            pattern = "父母均携带(需查是否同一变异/AR纯合已覆盖)"
        else:
            pattern = "其他(需人工核对基因型)"
        results.append(row | {'INHERITANCE_PATTERN': pattern})

opt_cols = header[11:11+nopt] if nopt else []   # BASE_FIELDS共11列(CHROM..NMD),之后紧跟可选dbNSFP列
out_cols = ['CHROM','POS','REF','ALT','ANN[0].GENE','ANN[0].EFFECT','ANN[0].IMPACT',
            'ANN[0].HGVS_C','ANN[0].HGVS_P','LOF','NMD'] + opt_cols + \
           ['CLNSIG','gnomad_popmax_af','GEN[0].GT','GEN[1].GT','GEN[2].GT','INHERITANCE_PATTERN']
disp_hdr = ['CHROM','POS','REF','ALT','GENE','EFFECT','IMPACT','HGVS_c','HGVS_p','LOF','NMD'] + \
           [c.replace('dbNSFP_','') for c in opt_cols] + \
           ['CLNSIG','gnomAD_popmax','GT_fetus','GT_mother','GT_father','INHERITANCE_PATTERN']

with open(out_file, 'w') as out:
    out.write('\t'.join(disp_hdr) + '\n')
    for row in sorted(results, key=lambda r: (r['ANN[0].GENE'], int(r['POS']))):
        out.write('\t'.join(row.get(c, '.') for c in out_cols) + '\n')

print(f"[filter] panel扫描: {len(results)} 条记录, 涉及 {len(by_gene)} 个基因", file=sys.stderr)
PYEOF

echo "[filter] ✅ 候选: $OUTDIR/tier1_panel_hits.tsv (panel内,含遗传来源判断) / candidates.tsv (genome-wide严格模型)"
echo "下一步: 把候选表交给 Claude 判读层 → 见 06_interpret_with_claude.md"
