from itertools import combinations

def wheel(n):
    m=n-1
    E=[(0,i) for i in range(1,n)] + [(i,1+(i % m)) for i in range(1,n)]
    return [tuple(sorted(e)) for e in E]

def neigh(n,E):
    N=[set() for _ in range(n)]
    for u,v in E:N[u].add(v);N[v].add(u)
    return N

def support(N,A):
    A=set(A); odd={v for v in range(len(N)) if sum(u in A for u in N[v])%2}
    return A|odd,odd

def target_min(n):
    E=wheel(n);N=neigh(n,E);best=n+1;As=[]
    for r in range(1,min(4,n)+1):
      for A in combinations(range(n),r):
        w=len(support(N,A)[0])
        if w<best:best=w;As=[A]
        elif w==best:As.append(A)
    return best,As

def spoke_low(n):
    E=wheel(n); e=(0,1); H=[f for f in E if f!=e]
    NG=neigh(n,E);NH=neigh(n,H);out=[]
    for r in range(1,4):
      for A in combinations(range(n),r):
        sup,zH=support(NH,A); _,zG=support(NG,A)
        if len(sup)<=3 and zH!=zG: out.append((A,len(sup)))
    return out

for n in range(5,21):
    dt,As=target_min(n)
    sl=spoke_low(n) if n>=6 else []
    print(f'W_{n}: d_targ={dt}; minA={As[:4]}; spoke_non_target_w<=3={sl}')
