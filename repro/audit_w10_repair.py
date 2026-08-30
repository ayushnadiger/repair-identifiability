"""Exact raw-repair obstruction for W10.
Expected: c_raw^rep(W10)=3, 23698 distinct nonzero masks.
"""
import multiprocessing as mp, random, sys
from audit_raw_repair_fast import wheel, make_graph, context_mask
from raw_repair_context import sep_mask as slow_sep_mask
from cdiag_verify import CdiagSolver

N=10
E=wheel(N)

def trits(x):
    b=[0]*N
    for i in range(N-1,-1,-1):
        b[i]=x%3; x//=3
    return tuple(b)

def chunk(args):
    lo,hi=args
    E0=wheel(N); T=make_graph(N,E0); C=[make_graph(N,[f for f in E0 if f!=e]) for e in E0]
    return [context_mask(N,E0,T,C,trits(i),True)[0] for i in range(lo,hi)]

def no_two_cover(masks,H=37):
    masks=list(dict.fromkeys(masks))
    npairs=H*(H-1)//2; full=(1<<npairs)-1
    B=(len(masks)+7)//8
    sep=[bytearray(B) for _ in range(npairs)]
    for i,M in enumerate(masks):
        by=i>>3; bt=1<<(i&7); x=M
        while x:
            lb=x&-x; p=lb.bit_length()-1; sep[p][by]|=bt; x^=lb
    sepI=[int.from_bytes(b,'little') for b in sep]
    counts=[x.bit_count() for x in sepI]
    allbits=(1<<len(masks))-1
    for M in masks:
        miss=full^M; ps=[];x=miss
        while x:
            lb=x&-x;p=lb.bit_length()-1;ps.append(p);x^=lb
        ps.sort(key=lambda p:counts[p])
        cand=allbits
        for p in ps:
            cand &= sepI[p]
            if not cand: break
        if cand:
            return False
    return True

def main():
    total=3**N
    bounds=[(0,15000),(15000,30000),(30000,45000),(45000,total)]
    with mp.Pool(4) as P:
        parts=P.map(chunk,bounds)
    masks=[m for part in parts for m in part]
    distinct=set(masks)
    nz={m for m in distinct if m}
    print('contexts',len(masks),'distinct_all',len(distinct),'distinct_nonzero',len(nz))
    ct=CdiagSolver(N,E).cdiag()
    print('c_targ^cut(W10)',ct)
    assert ct==2
    assert len(masks)==59049
    assert len(nz)==23698
    assert no_two_cover(masks), 'found a two-context repair schedule'
    print('exact two-cover obstruction: PASS')
    T=make_graph(N,E); C=[make_graph(N,[f for f in E if f!=e]) for e in E]
    rng=random.Random(9173)
    for i in range(120):
        b=tuple(rng.randrange(3) for _ in range(N))
        fast,H1=context_mask(N,E,T,C,b,True)
        slow,H2=slow_sep_mask(N,E,b,include_null=True)
        assert fast==slow and H1==H2, (i,b)
    print('120-context independent cross-check: PASS')
    print('RESULT c_raw^rep(W10)=3 (lower bound exact; upper bound is the 3-color rim certificate)')

if __name__=='__main__': main()
