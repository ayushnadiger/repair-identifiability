import itertools, time
from functools import lru_cache
from raw_repair_context import wheel, graph_stab

def kernel_basis(rows,n):
    rows=[r for r in rows if r]
    piv=[]; basis=[]
    r=0
    rows=rows[:]
    for c in range(n):
        p=next((i for i in range(r,len(rows)) if (rows[i]>>c)&1),None)
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and ((rows[i]>>c)&1): rows[i]^=rows[r]
        piv.append(c); r+=1
        if r==len(rows): break
    free=[c for c in range(n) if c not in piv]
    out=[]
    for f in free:
        x=1<<f
        for rr,c in zip(rows[:len(piv)],piv):
            if (rr>>f)&1: x|=1<<c
        out.append(x)
    return out

def canon_signed(gens):
    rows=[q | (s<<64) for q,s in gens if q]
    r=0
    for c in range(63):
        p=next((i for i in range(r,len(rows)) if (rows[i]>>c)&1),None)
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and ((rows[i]>>c)&1): rows[i]^=rows[r]
        r+=1
        if r==len(rows): break
    for x in rows:
        if (x & ((1<<64)-1))==0 and (x>>64)&1:
            raise RuntimeError('inconsistent signed subgroup')
    rows=[x for x in rows if x & ((1<<64)-1)]
    rows.sort(key=lambda x:(x & -x).bit_length())
    return tuple(rows)

def make_graph(n,E):
    nbr=[0]*n
    for u,v in E:
        nbr[u]|=1<<v; nbr[v]|=1<<u
    stab=graph_stab(n,E)
    return nbr,stab

def ctx_rows(n,nbr,b):
    rows=[]
    for v in range(n):
        if b[v]==0: rows.append(nbr[v])
        elif b[v]==1: rows.append(nbr[v]^(1<<v))
        else: rows.append(1<<v)
    return rows

def signature(n,nbr,stab,b,extra_zero=()):
    rows=ctx_rows(n,nbr,b)+[1<<v for v in extra_zero]
    bas=kernel_basis(rows,n)
    gens=[]
    for a in bas:
        z=0
        aa=a
        while aa:
            lb=aa&-aa; v=lb.bit_length()-1; z^=nbr[v]; aa^=lb
        q=a|z
        sign=0 if stab[(a,z)]==1 else 1
        gens.append((q,sign))
    return canon_signed(gens)

def context_mask(n,E,target,cuts,b,include_null=True):
    labels=[]; sig=[]
    if include_null:
        labels.append(('N',None)); sig.append(signature(n,*target,b))
    for j,e in enumerate(E):
        labels.append(('C',j)); sig.append(signature(n,*cuts[j],b))
        labels.append(('D',j)); sig.append(signature(n,*target,b,extra_zero=e))
    M=0;k=0
    for i in range(len(sig)):
        si=sig[i]
        for j in range(i+1,len(sig)):
            if si!=sig[j]: M|=1<<k
            k+=1
    return M,len(sig)

def solve_wheel(n):
    E=wheel(n); target=make_graph(n,E)
    cuts=[make_graph(n,[f for f in E if f!=e]) for e in E]
    d={}; H=None; t=time.time()
    for ii,b in enumerate(itertools.product(range(3),repeat=n)):
        M,H=context_mask(n,E,target,cuts,b,True)
        if M and M not in d: d[M]=b
        if ii and ii%10000==0: print('enum',ii,'distinct',len(d),'sec',round(time.time()-t,1),flush=True)
    full=(1<<(H*(H-1)//2))-1
    print('enum done contexts',3**n,'distinct',len(d),'sec',time.time()-t)
    if full in d: return 1,[d[full]],len(d),None
    ml=sorted(d,key=int.bit_count,reverse=True); mx=[]
    for idx,M in enumerate(ml):
        if not any((M|Q)==Q for Q,_ in mx): mx.append((M,d[M]))
        if idx and idx%5000==0: print('prune',idx,'max',len(mx),flush=True)
    print('maximal',len(mx))
    np=H*(H-1)//2
    cand=[[] for _ in range(np)]
    for j,(M,b) in enumerate(mx):
        x=M
        while x:
            lb=x&-x; cand[lb.bit_length()-1].append(j);x^=lb
    for i,(M,b) in enumerate(mx):
        miss=full^M
        if not miss: return 1,[b],len(d),len(mx)
        bits=[];x=miss
        while x:
            lb=x&-x;r=lb.bit_length()-1;bits.append(r);x^=lb
        r=min(bits,key=lambda r:len(cand[r]))
        for j in cand[r]:
            if (M|mx[j][0])==full:
                return 2,[b,mx[j][1]],len(d),len(mx)
        if i and i%1000==0: print('2cover',i,flush=True)
    return 3,None,len(d),len(mx)

def fmt(b):return ''.join('XYZ'[x] for x in b)
if __name__=='__main__':
    import sys
    for n in map(int,sys.argv[1:] or [8]):
        k,bs,nd,nm=solve_wheel(n)
        print('RESULT W',n,'k',k,'cert',None if not bs else [fmt(x) for x in bs],'distinct',nd,'maximal',nm)
