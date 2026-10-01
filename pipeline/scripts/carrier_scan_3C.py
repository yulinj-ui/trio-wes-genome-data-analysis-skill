#!/usr/bin/env python3
"""3-C 夫妇共同携带扫描 —— 【不看子代基因型】的独立扫描。
主流程一切筛选都隐含"先证者携带"前提，故 3-C 无法从已有结果推导，必须单独跑。
条件: popmax<=5%, 亲代 GQ>=30 & DP>=20, (HIGH impact | ClinVar P/LP | SpliceAI>=0.5)
形态: (a) 双方携带同一变异  (b) 同基因内各携带不同变异"""
import sys, csv
from collections import defaultdict
full, quadgt = sys.argv[1], sys.argv[2]
def cls(gt):
    if gt in ('./.','.|.','.',''): return 'miss'
    a=gt.replace('|','/').split('/')
    try: a=[int(x) for x in a]
    except ValueError: return 'miss'
    if len(a)!=2: return 'miss'
    return ('homref' if a[0]==0 else 'homalt') if a[0]==a[1] else 'het'
def sai(s):
    if not s or s=='.': return None
    b=None
    for rec in s.split(','):
        f=rec.split('|')
        if len(f)<7: continue
        for v in f[2:6]:
            try: v=float(v)
            except ValueError: continue
            b=v if b is None else max(b,v)
    return b
def fnum(x):
    try: return float(x)
    except (TypeError,ValueError): return None
def ok_qual(fld):
    p=fld.split(':')
    try: return int(p[2])>=20 and int(p[3])>=30
    except (IndexError,ValueError): return False

q={}
for line in open(quadgt):
    f=line.rstrip('\n').split('\t'); q[(f[0],f[1],f[2],f[3])]=f[4:8]
rows=[]
for r in csv.DictReader(open(full),delimiter='\t'):
    qq=q.get((r['CHROM'],r['POS'],r['REF'],r['ALT']))
    if not qq: continue
    r['_F2f'],r['_MOf'],r['_FAf'],r['_F1f']=qq
    rows.append(r)

cand=[]
for r in rows:
    af=fnum(r['gnomad_popmax_af'])
    if af is not None and af>=0.05: continue
    if not (ok_qual(r['_MOf']) and ok_qual(r['_FAf'])): continue
    c=r['CLNSIG']; s=sai(r['SpliceAI'])
    plp = c and (c.startswith('Pathogenic') or c.startswith('Likely_pathogenic')) and 'Conflicting' not in c
    if not (r['ANN[0].IMPACT']=='HIGH' or plp or (s is not None and s>=0.5)): continue
    cand.append(r)
print(f"[3C] 通过条件的变异 {len(cand)} 条（不看先证者基因型）")

bygene=defaultdict(lambda: {'mo':[], 'fa':[], 'both':[]})
for r in cand:
    mo,fa=cls(r['GEN[1].GT']),cls(r['GEN[2].GT'])
    g=r['ANN[0].GENE']
    mo_c = mo in ('het','homalt'); fa_c = fa in ('het','homalt')
    if mo_c and fa_c: bygene[g]['both'].append(r)
    elif mo_c: bygene[g]['mo'].append(r)
    elif fa_c: bygene[g]['fa'].append(r)

def line(r,who):
    f2=cls(r['GEN[0].GT']); f1=cls(r['_F1f'].split(':')[0])
    inh=[]
    if f2 in('het','homalt'): inh.append('胎2:'+f2)
    if f1 in('het','homalt'): inh.append('胎1:'+f1)
    print(f"     {who} {r['CHROM']}:{r['POS']} {r['REF']}>{r['ALT'][:14]} {r['ANN[0].EFFECT'][:26]} {r['ANN[0].HGVS_C']} {r['ANN[0].HGVS_P']}")
    print(f"        popmax={r['gnomad_popmax_af']} ClinVar={r['CLNSIG'][:40]} SpliceAI={sai(r['SpliceAI'])}")
    print(f"        子代遗传: {'、'.join(inh) if inh else '两胎均未遗传'}   AD 母|父 = {r['_MOf']} | {r['_FAf']}")

print("\n──── (a) 夫妇携带【同一】变异 ────")
na=0
for g,d in sorted(bygene.items()):
    for r in d['both']:
        print(f"  ★ {g}"); line(r,'双方'); na+=1
print(f"  → {na} 条" if na else "  → 本形态为空")
print("\n──── (b) 同一基因内夫妇【各携带不同】变异 ────")
nb=0
for g,d in sorted(bygene.items()):
    if d['mo'] and d['fa']:
        print(f"  ★ {g}")
        for r in d['mo']: line(r,'母源')
        for r in d['fa']: line(r,'父源')
        nb+=1
print(f"  → {nb} 个基因" if nb else "  → 本形态为空")
