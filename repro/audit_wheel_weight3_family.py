import itertools

def wheel(n):
    m=n-1
    E=[(0,i) for i in range(1,n)]
    E += [(i,1+(i % m)) for i in range(1,n)]
    return [tuple(sorted(e)) for e in E]

def stab_map(n,E):
    N=[0]*n
    for u,v in E:
        N[u]|=1<<v; N[v]|=1<<u
    d={}
    for a in range(1<<n):
        z=0
        for v in range(n):
            if ((N[v]&a).bit_count()&1): z|=1<<v
        q=0
        for u,v in E:
            if (a>>u)&1 and (a>>v)&1: q^=1
        d[(a,z)] = -1 if q else 1
    return d

def build_context_masks(n,w=3):
    E=wheel(n); T=stab_map(n,E)
    cuts=[stab_map(n,[f for f in E if f!=e]) for e in E]
    H=[T]+cuts
    pairs=list(itertools.combinations(range(len(H)),2)); full=(1<<len(pairs))-1
    by_req={}
    for r in range(1,w+1):
        for supp in itertools.combinations(range(n),r):
            for labs in itertools.product((1,2,3),repeat=r):
                x=z=0
                for v,l in zip(supp,labs):
                    if l in (1,2): x|=1<<v
                    if l in (2,3): z|=1<<v
                vals=[S.get((x,z),0) for S in H]
                mask=0
                for k,(i,j) in enumerate(pairs):
                    if vals[i]!=vals[j]: mask|=1<<k
                if mask:
                    req=tuple(zip(supp,labs))
                    by_req[req]=by_req.get(req,0)|mask
    masks={}
    for req,pmask in by_req.items():
        fixed={v:l for v,l in req}; free=[v for v in range(n) if v not in fixed]
        for rest in itertools.product((1,2,3), repeat=len(free)):
            b=[0]*n
            for v,l in fixed.items(): b[v]=l
            for v,l in zip(free,rest): b[v]=l
            bt=tuple(b)
            masks[bt]=masks.get(bt,0)|pmask
    return masks,full

def prune(masks):
    rev={}
    for b,M in masks.items():
        if M: rev.setdefault(M,b)
    items=sorted(rev.items(),key=lambda kv:kv[0].bit_count(),reverse=True)
    maximal=[]
    for M,b in items:
        if not any((M|Q)==Q for Q,_ in maximal):
            maximal.append((M,b))
    return rev,maximal

def min_cover(masks,full):
    rev,maximal=prune(masks)
    for M,b in maximal:
        if M==full: return 1,[b],len(rev),len(maximal)
    for i,(M,b) in enumerate(maximal):
        need=full & ~M
        for Q,bq in maximal[i:]:
            if need & ~Q==0:
                return 2,[b,bq],len(rev),len(maximal)
    return 3,None,len(rev),len(maximal)

def fmt(b): return ''.join('XYZ'[x-1] for x in b)

for n in range(6,14):
    masks,full=build_context_masks(n,3)
    k,cert,nd,nm=min_cover(masks,full)
    print(f'W_{n}: c_raw_cut_<=3={k}; chi_rim={2 if (n-1)%2==0 else 3}; distinct={nd}; maximal={nm}; cert={[fmt(x) for x in cert] if cert else None}')
