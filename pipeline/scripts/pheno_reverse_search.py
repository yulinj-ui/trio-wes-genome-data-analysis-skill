#!/usr/bin/env python3
"""表型优先反向搜索 + 全基因组 ClinVar P/LP + HIGH impact 扫描（四人）"""
import sys, csv
from collections import defaultdict
full, quadgt, panelf = sys.argv[1], sys.argv[2], sys.argv[3]

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

panel=set(l.strip() for l in open(panelf) if l.strip() and not l.startswith('#'))
q={}
for line in open(quadgt):
    f=line.rstrip('\n').split('\t'); q[(f[0],f[1],f[2],f[3])]=f[4:8]
rows=[]
for r in csv.DictReader(open(full),delimiter='\t'):
    k=(r['CHROM'],r['POS'],r['REF'],r['ALT']); qq=q.get(k)
    if not qq: continue
    r['_F2f'],r['_MOf'],r['_FAf'],r['_F1f']=qq
    r['_F1']=qq[3].split(':')[0]
    rows.append(r)

def origin(r):
    f2,mo,fa,f1=cls(r['GEN[0].GT']),cls(r['GEN[1].GT']),cls(r['GEN[2].GT']),cls(r['_F1'])
    if mo in('het','homalt') and fa=='homref': o='母源'
    elif fa in('het','homalt') and mo=='homref': o='父源'
    elif mo in('het','homalt') and fa in('het','homalt'): o='双亲均携带'
    elif mo=='homref' and fa=='homref': o='★de novo'
    else: o='?'
    return f2,mo,fa,f1,o
def show(r,tag):
    f2,mo,fa,f1,o=origin(r)
    who=[]
    if f2 in('het','homalt'): who.append('胎2')
    if f1 in('het','homalt'): who.append('胎1')
    print(f"  [{tag}] {r['ANN[0].GENE']:<12} {r['CHROM']}:{r['POS']} {r['REF']}>{r['ALT']}  {r['ANN[0].EFFECT'][:34]}")
    print(f"       {r['ANN[0].HGVS_C']} {r['ANN[0].HGVS_P']}  携带者={'+'.join(who)}  来源={o}")
    print(f"       popmax={r['gnomad_popmax_af']} CADD={r['CADD_phred']} REVEL={r['REVEL_score']} AM={r['AlphaMissense_pred']} SpliceAI={sai(r['SpliceAI'])} ClinVar={r['CLNSIG'][:45]}")
    print(f"       AD 胎2|母|父|胎1 = {r['_F2f']} | {r['_MOf']} | {r['_FAf']} | {r['_F1f']}")

print("="*80); print("【A】表型基因集内 · 胎儿携带 · popmax≤2% · 不限后果类型")
n=0
for r in rows:
    if r['ANN[0].GENE'] not in panel: continue
    f2,mo,fa,f1,o=origin(r)
    if f2 not in('het','homalt') and f1 not in('het','homalt'): continue
    af=fnum(r['gnomad_popmax_af'])
    if af is not None and af>=0.02: continue
    s=sai(r['SpliceAI'])
    imp=r['ANN[0].IMPACT']
    if imp in ('HIGH',) or (s is not None and s>=0.2) or \
       (imp=='MODERATE' and ((fnum(r['REVEL_score']) or 0)>=0.5 or (fnum(r['CADD_phred']) or 0)>=25)) or \
       ('athogenic' in r['CLNSIG']):
        show(r,'表型集'); n+=1
print(f"  → 共 {n} 条\n")

print("="*80); print("【B】全基因组 · ClinVar Pathogenic/Likely_pathogenic · 胎儿携带")
n=0
for r in rows:
    c=r['CLNSIG']
    if not c or c=='.': continue
    if 'Conflicting' in c: continue
    if not (c.startswith('Pathogenic') or c.startswith('Likely_pathogenic') or c=='Pathogenic/Likely_pathogenic'): continue
    f2,mo,fa,f1,o=origin(r)
    if f2 not in('het','homalt') and f1 not in('het','homalt'): continue
    af=fnum(r['gnomad_popmax_af'])
    if af is not None and af>=0.02: continue
    show(r,'ClinVar'); n+=1
print(f"  → 共 {n} 条\n")

print("="*80); print("【C】全基因组 · HIGH impact(LoF) · popmax≤0.1% · 胎儿携带 · 非父母共有多态")
n=0
for r in rows:
    if r['ANN[0].IMPACT']!='HIGH': continue
    f2,mo,fa,f1,o=origin(r)
    if f2 not in('het','homalt') and f1 not in('het','homalt'): continue
    af=fnum(r['gnomad_popmax_af'])
    if af is not None and af>=0.001: continue
    if o=='双亲均携带': continue
    show(r,o); n+=1
print(f"  → 共 {n} 条")
