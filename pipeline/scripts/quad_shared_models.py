#!/usr/bin/env python3
"""四人家系交集分析 —— 复发病例的核心判据。
两个子代均受累、父母表型正常时，找【两个子代共有而父母不共有】的致病候选。
六种模型: ①X连锁隐性共有 ②共有de novo(提示生殖腺嵌合) ③AR纯合共有
          ④母源杂合共有(显性外显不全) ⑤父源杂合共有 ⑥反式复合杂合共有
用法: quad_shared_models.py <all_variants_full.tsv> <quad_gt.tsv> <out_prefix>
  all_variants_full.tsv 由 05_inheritance_filter.sh 产出（trio 三列基因型 + 注释）
  quad_gt.tsv 由 bcftools query 从 quad.filtered.vcf.gz 取四列基因型，列序须为
  先证者/母/父/第二子代（即 config 的 ID_FETUS,ID_MOTHER,ID_FATHER,ID_FETUS1 字母序）
⛔ 输出的每一条在写进报告前都必须过伪影筛（回原始 BAM 复核 AD），见 runbook。"""
import sys, csv
from collections import defaultdict

full, quadgt, outp = sys.argv[1], sys.argv[2], sys.argv[3]
PAR = [(10001, 2781479), (155701383, 156030895)]          # GRCh38 chrX 拟常染色体区
def in_par(pos): return any(a <= pos <= b for a, b in PAR)

def cls(gt):
    if gt in ('./.', '.|.', '.', ''): return 'miss'
    a = gt.replace('|', '/').split('/')
    try: a = [int(x) for x in a]
    except ValueError: return 'miss'
    if len(a) != 2: return 'miss'
    if a[0] == a[1]: return 'homref' if a[0] == 0 else 'homalt'
    return 'het'

def spliceai_max(s):
    if not s or s == '.': return None
    best = None
    for rec in s.split(','):
        f = rec.split('|')
        if len(f) < 7: continue
        for v in f[2:6]:
            try: v = float(v)
            except ValueError: continue
            best = v if best is None else max(best, v)
    return best

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None

# 胎1(P4) 基因型 + AD/DP/GQ
q = {}
for line in open(quadgt):
    f = line.rstrip('\n').split('\t')
    key = (f[0], f[1], f[2], f[3])
    q[key] = {'F2': f[4], 'MO': f[5], 'FA': f[6], 'F1': f[7]}

rows = []
with open(full) as fh:
    rd = csv.DictReader(fh, delimiter='\t')
    for r in rd:
        key = (r['CHROM'], r['POS'], r['REF'], r['ALT'])
        qq = q.get(key)
        if not qq: continue
        r['_F1'] = qq['F1'].split(':')[0]
        r['_F1_full'] = qq['F1']
        r['_F2_full'] = qq['F2']; r['_MO_full'] = qq['MO']; r['_FA_full'] = qq['FA']
        rows.append(r)
print(f"[join] all_variants_full={sum(1 for _ in open(full))-1}  成功匹配胎1基因型={len(rows)}")

def keep(r, maxaf=0.01):
    af = fnum(r['gnomad_popmax_af'])
    if af is not None and af >= 0 and af >= maxaf: return False
    sa = spliceai_max(r['SpliceAI'])
    return r['ANN[0].IMPACT'] in ('HIGH', 'MODERATE') or (sa is not None and sa >= 0.5)

def gts(r):
    return cls(r['GEN[0].GT']), cls(r['GEN[1].GT']), cls(r['GEN[2].GT']), cls(r['_F1'])

models = defaultdict(list)
for r in rows:
    if not keep(r): continue
    f2, mo, fa, f1 = gts(r)
    if 'miss' in (f2, mo, fa, f1): continue
    chrom, pos = r['CHROM'], int(r['POS'])
    carry2 = f2 in ('het', 'homalt'); carry1 = f1 in ('het', 'homalt')

    # 1) X连锁隐性：两胎均携带(男胎半合应为homalt)，母杂合，父不携带
    if chrom == 'chrX' and not in_par(pos):
        if carry2 and carry1 and mo == 'het' and fa == 'homref':
            models['XL_shared'].append(r)
    # 2) 两胎共有 de novo（父母外周血均阴性 → 提示亲代生殖腺嵌合）
    if carry2 and carry1 and mo == 'homref' and fa == 'homref':
        models['sharedDeNovo'].append(r)
    # 3) 两胎共有 AR 纯合
    if f2 == 'homalt' and f1 == 'homalt' and mo == 'het' and fa == 'het' and chrom != 'chrX':
        models['AR_hom_shared'].append(r)
    # 4) 两胎共有的母源杂合（显性外显不全 / 母方嵌合）
    if f2 == 'het' and f1 == 'het' and mo == 'het' and fa == 'homref' and chrom not in ('chrX', 'chrY'):
        models['matDominant_shared'].append(r)
    # 5) 两胎共有的父源杂合
    if f2 == 'het' and f1 == 'het' and fa == 'het' and mo == 'homref' and chrom not in ('chrX', 'chrY'):
        models['patDominant_shared'].append(r)

# 6) 两胎共有的复合杂合（反式）
bygene = defaultdict(list)
for r in rows:
    if not keep(r): continue
    f2, mo, fa, f1 = gts(r)
    if 'miss' in (f2, mo, fa, f1): continue
    if f2 == 'het' and f1 == 'het':
        bygene[r['ANN[0].GENE']].append((r, mo, fa))
for g, vs in bygene.items():
    mat = [v for v in vs if v[1] == 'het' and v[2] == 'homref']
    pat = [v for v in vs if v[2] == 'het' and v[1] == 'homref']
    for a in mat:
        for b in pat:
            models['compHet_shared'].append(a[0]); models['compHet_shared'].append(b[0])

hdrs = ['CHROM','POS','REF','ALT','ANN[0].GENE','ANN[0].EFFECT','ANN[0].IMPACT','ANN[0].HGVS_C',
        'ANN[0].HGVS_P','CADD_phred','REVEL_score','AlphaMissense_pred','SpliceAI','CLNSIG',
        'gnomad_popmax_af','hiConfDeNovo','_F2_full','_MO_full','_FA_full','_F1_full']
name = {'XL_shared':'① X连锁隐性·两胎共有','sharedDeNovo':'② 两胎共有de novo(提示生殖腺嵌合)',
        'AR_hom_shared':'③ AR纯合·两胎共有','matDominant_shared':'④ 母源杂合·两胎共有(显性外显不全)',
        'patDominant_shared':'⑤ 父源杂合·两胎共有','compHet_shared':'⑥ 复合杂合(反式)·两胎共有'}
with open(outp + '.tsv', 'w') as out:
    out.write('MODEL\t' + '\t'.join(hdrs) + '\n')
    for m in ['XL_shared','sharedDeNovo','AR_hom_shared','compHet_shared','matDominant_shared','patDominant_shared']:
        vs = models.get(m, [])
        seen = set(); uniq = []
        for r in vs:
            k = (r['CHROM'], r['POS'], r['REF'], r['ALT'])
            if k in seen: continue
            seen.add(k); uniq.append(r)
        print(f"\n{'='*70}\n{name[m]}  —— {len(uniq)} 条")
        for r in uniq:
            out.write(m + '\t' + '\t'.join(str(r.get(h, '.')) for h in hdrs) + '\n')
            sa = spliceai_max(r['SpliceAI'])
            print(f"  {r['ANN[0].GENE']:<12} {r['CHROM']}:{r['POS']} {r['REF']}>{r['ALT']}  {r['ANN[0].EFFECT'][:30]}")
            print(f"      {r['ANN[0].HGVS_C']} {r['ANN[0].HGVS_P']}  popmax={r['gnomad_popmax_af']} CADD={r['CADD_phred']} REVEL={r['REVEL_score']} SpliceAI_max={sa}")
            print(f"      GT(胎2/母/父/胎1)= {r['_F2_full']} | {r['_MO_full']} | {r['_FA_full']} | {r['_F1_full']}")
            if r['CLNSIG'] != '.': print(f"      ClinVar={r['CLNSIG']}")
print(f"\n已写出: {outp}.tsv")
