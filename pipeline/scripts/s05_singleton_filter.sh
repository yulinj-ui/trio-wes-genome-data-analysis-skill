#!/usr/bin/env bash
# ============================================================================
# s05_singleton_filter.sh —— 单人（singleton）候选变异筛选
#
# ⛔ 与 05_inheritance_filter.sh 的本质区别：没有父母基因型。
#    因此本脚本【不产出】de novo 判定，也【不宣称】复合杂合成立：
#      · AD 候选     = 罕见杂合（来源不明：可能 de novo，也可能遗传自表型正常/轻表现的亲代）
#      · AR-hom      = 纯合（无父母佐证，须防"半合于缺失之上"与技术性假纯合）
#      · XL-hemi     = 男性 chrX 非 PAR 区半合
#      · comphet?    = 同基因 ≥2 个罕见杂合 —— 仅"疑似"，相位未验证（相位门）
#      · ClinVar-P/LP= 无论后果等级一律捞出（捕获同义/剪接/内含子等被 IMPACT 漏掉的已知致病位点）
# 输出交给阶段⑤⑥ 判读层。
# ============================================================================
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "$DIR/config.sh"
[ "${FAMILY_MODE:-trio}" = "singleton" ] || { echo "❌ 本脚本仅用于 singleton" >&2; exit 1; }
activate_env
_JB="$JAVA_ENV_BIN"; [ -x "$_JB/java" ] && export PATH="$_JB:$PATH"
export _JAVA_OPTIONS="-Xmx${MAXMEM_GB}g"
VP="${VCF_PREFIX:-proband}"
VCF="$RESULTDIR/${VP}.annot.vcf.gz"
[ -f "$VCF" ] || { echo "❌ 找不到 $VCF，先跑 04_annotate.sh" >&2; exit 1; }
GNOMAD_ZIP="$REFDIR/GRCh38/gnomad.hg38.genomes.v3.fix.zip"
OUTDIR="$RESULTDIR/candidates"; mkdir -p "$OUTDIR"; LOG="$LOGDIR/s05.log"

# ---- A) slivar gnotate：加 gnomad_popmax_af（不做任何基因型过滤，保留全部变异）----
echo "[filter] slivar gnotate 加 gnomAD popmax…"
GN=""; [ -f "$GNOMAD_ZIP" ] && GN="-g $GNOMAD_ZIP" || echo "⚠️ gnomAD gnotate 库缺失，频率列将全为 . —— 频率过滤失效，须在报告中声明" >&2
slivar expr --vcf "$VCF" $GN --info "true" --out-vcf "$WORKDIR/${VP}.gn.vcf" 2>>"$LOG"
NGN=$(grep -vc '^#' "$WORKDIR/${VP}.gn.vcf" || true)
[ "${NGN:-0}" -gt 0 ] || { echo "❌ slivar 输出 0 条变异 —— 静默失效，停" >&2; exit 1; }
echo "[filter] slivar 输出 $NGN 条变异"

# ---- B) 导出全量变异表（单样本：只有 GEN[0]）----
HDR_FIELDS=$(bcftools view -h "$WORKDIR/${VP}.gn.vcf" | grep "^##INFO" | sed -E 's/^##INFO=<ID=([^,]+).*/\1/')
BASE='CHROM POS REF ALT "ANN[0].GENE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" LOF NMD'
OPT=""; OPTH=""
for f in CADD_phred REVEL_score AlphaMissense_pred AlphaMissense_score MetaRNN_pred SIFT_pred Polyphen2_HDIV_pred \
         dbNSFP_CADD_phred dbNSFP_REVEL_score dbNSFP_AlphaMissense_pred dbNSFP_MetaRNN_pred SpliceAI; do
  echo "$HDR_FIELDS" | grep -qx "$f" && { OPT="$OPT $f"; OPTH="$OPTH\t${f#dbNSFP_}"; }
done
TAIL='CLNSIG CLNREVSTAT CLNDN gnomad_popmax_af "GEN[0].GT" "GEN[0].AD" "GEN[0].DP" "GEN[0].GQ" FILTER'
echo "[filter] 导出全量变异表…"
eval SnpSift extractFields -e '"."' "\"\$WORKDIR/\${VP}.gn.vcf\"" $BASE $OPT $TAIL 2>>"$LOG" > "$WORKDIR/all_variants_full.tsv"
NROW=$(( $(wc -l < "$WORKDIR/all_variants_full.tsv") - 1 ))
echo "[filter] 全量表 $NROW 行"
[ "$NROW" -ge 20000 ] || echo "⚠️ 全量表仅 $NROW 行，低于 WES 经验下界，先怀疑静默失效" >&2

# ---- C) 单人遗传模型分类 ----
python3 - "$WORKDIR/all_variants_full.tsv" "$PANEL_FILE" "${PHENO_PANEL:-}" "$OUTDIR" "$SEX_FETUS" <<'PYEOF'
# -*- coding: utf-8 -*-
import csv,sys,collections
full,panel_file,pheno_file,outdir,sex = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]

pheno={}   # gene -> 匹配的 HPO 术语数（tier-1 排序用）
if pheno_file:
    for l in open(pheno_file,encoding='utf-8'):
        if l.startswith('#') or not l.strip(): continue
        p=l.rstrip('\n').split('\t')
        tags=[t for t in (p[1].split(',') if len(p)>1 else []) if t.startswith('HP:')]
        pheno[p[0]]=len(tags)

def gtclass(gt):
    if gt in ('.','./.','.|.',''): return 'missing'
    a=gt.replace('|','/').split('/')
    try: a=[int(x) for x in a]
    except: return 'missing'
    if len(a)==1: return 'hemi_ref' if a[0]==0 else 'hemi_alt'
    if len(set(a))==1: return 'hom_ref' if a[0]==0 else 'hom_alt'
    return 'het'

def af(v):
    try: return float(v)
    except: return 0.0        # gnotate '.' = 不在 gnomAD = 按极罕见处理

def spliceai_max(v):
    if not v or v=='.': return 0.0
    best=0.0
    for rec in v.split(','):
        f=rec.split('|')
        if len(f)>=6:
            for x in f[2:6]:
                try: best=max(best,float(x))
                except: pass
    return best

rows=list(csv.DictReader(open(full,encoding='utf-8'),delimiter='\t'))
G='ANN[0].GENE'; IMP='ANN[0].IMPACT'; EFF='ANN[0].EFFECT'
PAR_X=[(10001,2781479),(155701383,156030895)]   # GRCh38 chrX PAR1/PAR2

def in_par(pos): return any(a<=pos<=b for a,b in PAR_X)

hits=[]
gene_het=collections.defaultdict(list)
for r in rows:
    gt=gtclass(r['GEN[0].GT'])
    if gt in ('missing','hom_ref','hemi_ref'): continue
    a=af(r['gnomad_popmax_af']); imp=r[IMP]; sp=spliceai_max(r.get('SpliceAI','.'))
    cln=(r.get('CLNSIG') or '.')
    plp = ('Pathogenic' in cln or 'Likely_pathogenic' in cln) and 'Conflicting' not in cln
    damaging = imp in ('HIGH','MODERATE') or sp>=0.20 or plp
    if not damaging: continue
    models=[]
    if gt=='het' and a<0.0001: models.append('AD?罕见杂合(来源不明:无父母数据)')
    if gt=='hom_alt' and a<0.005: models.append('AR-hom(纯合;无父母佐证,须排除半合/假纯合)')
    if r['CHROM']=='chrX' and sex=='1' and not in_par(int(r['POS'])) and gt in ('hom_alt','hemi_alt','het') and a<0.0001:
        models.append('XL-hemi(男性chrX非PAR)')
    if plp and a<0.01: models.append(f'ClinVar-{cln}')
    if gt=='het' and a<0.01: gene_het[r[G]].append(r['POS'])
    if models:
        r=dict(r); r['MODEL']=';'.join(models); r['SpliceAI_DSmax']=f"{sp:.2f}"
        r['PHENO_HPO_N']=str(pheno.get(r[G],0)); hits.append(r)

# 疑似复合杂合（相位未验证）：同基因 ≥2 个罕见杂合
comphet={g:v for g,v in gene_het.items() if len(v)>=2}
for r in hits:
    if r[G] in comphet and 'het' in gtclass(r['GEN[0].GT']):
        r['MODEL']=r['MODEL']+f";疑似复合杂合(同基因{len(comphet[r[G]])}处杂合,⚠️相位未验证)"
# 补：只因 comphet 才入选的（单条频率未过 AD 阈值）
seen={(r['CHROM'],r['POS']) for r in hits}
for r in rows:
    if r[G] not in comphet: continue
    if (r['CHROM'],r['POS']) in seen: continue
    if gtclass(r['GEN[0].GT'])!='het': continue
    a=af(r['gnomad_popmax_af']); sp=spliceai_max(r.get('SpliceAI','.'))
    if a>=0.01: continue
    if not (r[IMP] in ('HIGH','MODERATE') or sp>=0.20): continue
    r=dict(r); r['MODEL']=f"疑似复合杂合(同基因{len(comphet[r[G]])}处杂合,⚠️相位未验证)"
    r['SpliceAI_DSmax']=f"{sp:.2f}"; r['PHENO_HPO_N']=str(pheno.get(r[G],0)); hits.append(r)

cols=[c for c in rows[0].keys()]+['SpliceAI_DSmax','PHENO_HPO_N','MODEL']
def write(fn,data):
    with open(fn,'w',encoding='utf-8') as f:
        f.write('\t'.join(cols)+'\n')
        for r in data: f.write('\t'.join(str(r.get(c,'.')) for c in cols)+'\n')

hits.sort(key=lambda r:(-int(r['PHENO_HPO_N']), r[G], int(r['POS'])))
write(f"{outdir}/candidates_singleton.tsv",hits)
tier1=[r for r in hits if int(r['PHENO_HPO_N'])>0]
write(f"{outdir}/tier1_pheno_hits.tsv",tier1)
mc=collections.Counter(m for r in hits for m in r['MODEL'].split(';'))
print(f"[filter] 候选总数 {len(hits)}；表型 panel 内(tier-1) {len(tier1)}",file=sys.stderr)
for k,v in mc.most_common(): print(f"          {k}: {v}",file=sys.stderr)
PYEOF
echo "[filter] ✅ $OUTDIR/candidates_singleton.tsv （全基因）/ tier1_pheno_hits.tsv （表型 panel 内）"
