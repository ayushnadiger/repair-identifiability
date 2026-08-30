"""Symbolic verification of the two-context wheel certificate formulas.

Criterion (quotient edge labels): for context beta on G, let Q = F_2^V / row(M_beta).
Edge label L(uv) = span{[e_u],[e_v]} in Q. Two contexts diagnose iff the ordered
pairs (L1(e), L2(e)) are pairwise distinct over e in E and never (0,0).
"""
import itertools

def wheel(n):
    m = n - 1
    E = [(0, i) for i in range(1, n)] + [(i, (i % m) + 1) for i in range(1, m + 1)]
    return n, sorted(set(tuple(sorted(e)) for e in E))

def rowspace_basis(rows, n):
    basis=[]
    for r in rows:
        cur=r
        for b in basis:
            if cur ^ b < cur: cur ^= b
        if cur:
            basis.append(cur); basis.sort(reverse=True)
    return basis

def reduce_vec(v,basis):
    for b in basis:
        if v ^ b < v: v ^= b
    return v

def context_rows(n,E,beta):
    Gamma=[0]*n
    for u,v in E:
        Gamma[u]|=1<<v; Gamma[v]|=1<<u
    rows=[]
    for v in range(n):
        if beta[v]==0: rows.append(Gamma[v])
        elif beta[v]==1: rows.append(Gamma[v]^(1<<v))
        else: rows.append(1<<v)
    return rows

def edge_labels(n,E,beta):
    basis=rowspace_basis(context_rows(n,E,beta),n)
    cls=[reduce_vec(1<<v,basis) for v in range(n)]
    return {e:tuple(sorted({cls[e[0]],cls[e[1]],cls[e[0]]^cls[e[1]]}-{0})) for e in E}

def pair_diagnoses(n,E,b1,b2):
    L1,L2=edge_labels(n,E,b1),edge_labels(n,E,b2); seen=set()
    for e in E:
        key=(L1[e],L2[e])
        if key==((),()) or key in seen:return False
        seen.add(key)
    return True

D={'X':0,'Y':1,'Z':2}
def mk(hub,rim):return tuple(D[c] for c in hub+rim)
def formulas(n):
    m=n-1
    if n%6==2:
        t=(m-1)//6; return mk('X','ZYY'*(2*t)+'Z'),mk('X','YYZ'*(2*t)+'Z')
    if n%6==4:
        s=m//3;base='YYZ'*s;return mk('Y',base),mk('Y',base[1:]+base[0])
    raise ValueError(n)

if __name__=='__main__':
    from cdiag_verify import CdiagSolver,wheel as bw
    import random
    random.seed(3)
    for nn in (6,8,10):
        n,E=bw(nn);S=CdiagSolver(n,E);full=(1<<S.npairs)-1
        for _ in range(300):
            b1=tuple(random.randrange(3) for _ in range(n));b2=tuple(random.randrange(3) for _ in range(n))
            assert ((S.sep_mask(b1)|S.sep_mask(b2))==full)==pair_diagnoses(n,E,b1,b2)
        print(f'label criterion == brute force on W_{nn} (300 random pairs)')
    ok=[]
    for n in range(8,65,2):
        if n%6==0:continue
        nn,E=wheel(n);b1,b2=formulas(n);assert pair_diagnoses(n,E,b1,b2);ok.append(n)
    print('certificate formulas verified for W_n, n =',ok)
