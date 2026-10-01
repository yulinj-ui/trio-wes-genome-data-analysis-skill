#!/usr/bin/env python3
"""KING-robust 亲缘系数。用法: king_kinship.py <gt.tsv> <sample1,sample2,...>
gt.tsv 格式: CHROM POS [GT DP GQ]×N
φ_ij = [N(Aa,Aa) - 2*N(AA,aa)] / [N(Aa_i) + N(Aa_j)]   不依赖等位基因频率面板。"""
import sys, itertools
tsv = sys.argv[1]; S = [x for x in sys.argv[2].split(',') if x]; N = len(S)

def code(g):
    if g in ('./.', '.|.', '.'): return None
    a = g.replace('|', '/').split('/')
    try: a = [int(x) for x in a]
    except ValueError: return None
    return sum(a) if len(a) == 2 else None

rows = []
for line in open(tsv):
    f = line.rstrip('\n').split('\t'); gts = []; ok = True
    for i in range(N):
        g, dp, gq = f[2 + i*3], f[3 + i*3], f[4 + i*3]
        try:
            if int(dp) < 20 or int(gq) < 30: ok = False; break
        except ValueError: ok = False; break
        c = code(g)
        if c is None: ok = False; break
        gts.append(c)
    if ok: rows.append(gts)

n = len(rows)
if n < 2000:
    sys.exit(f"❌ 有效位点仅 {n} 个 (<2000)，KING 估计不可靠，拒绝输出")
print(f"KING-robust 亲缘系数  (双等位常染色体 SNV, DP>=20, GQ>=30)  有效位点 n={n}")
print(f"{'配对':<22}{'kinship':>10}{'IBS0率':>10}  判定")
for i, j in itertools.combinations(range(N), 2):
    Nhh  = sum(1 for r in rows if r[i] == 1 and r[j] == 1)
    Nopp = sum(1 for r in rows if (r[i] == 0 and r[j] == 2) or (r[i] == 2 and r[j] == 0))
    Ni   = sum(1 for r in rows if r[i] == 1); Nj = sum(1 for r in rows if r[j] == 1)
    phi  = (Nhh - 2*Nopp) / (Ni + Nj); ibs0 = Nopp / n
    if   phi > 0.354:  rel = "同一人/单卵双胎"
    elif phi > 0.177:  rel = "一级亲·亲子" if ibs0 < 0.002 else "一级亲·全同胞"
    elif phi > 0.0884: rel = "二级亲"
    elif phi > 0.0442: rel = "三级亲(一级表亲)"
    elif phi > 0.0221: rel = "四级亲(二级表亲)"
    else:              rel = "无可检出血缘关系"
    print(f"{S[i] + ' - ' + S[j]:<22}{phi:>10.4f}{ibs0:>10.4f}  {rel}")
print("\n判定阈值(KING标准): 一级>0.177  二级>0.0884  三级(一级表亲)>0.0442  四级>0.0221")
print("\n各样本 杂合/非ref纯合 比 (近亲子代应显著低于父母):")
for k in range(N):
    het = sum(1 for r in rows if r[k] == 1); hom = sum(1 for r in rows if r[k] == 2)
    print(f"  {S[k]:<10} het={het:<7} homalt={hom:<7} ratio={het/hom:.3f}")
