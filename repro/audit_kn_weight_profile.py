from itertools import combinations

def complete(n):
    return [(i,j) for i in range(n) for j in range(i+1,n)]

def neigh(n,E):
    N=[set() for _ in range(n)]
    for u,v in E:N[u].add(v);N[v].add(u)
    return N

def stab_support(N,A):
    A=set(A); odd={v for v in range(len(N)) if sum(u in A for u in N[v])%2}
    return A|odd,odd

def mins(n):
    E=complete(n);NG=neigh(n,E)
    dt=n+1
    for r in range(1,n+1):
      for A in combinations(range(n),r):dt=min(dt,len(stab_support(NG,A)[0]))
    dfs=[]
    for e in E:
      H=[f for f in E if f!=e];NH=neigh(n,H);de=n+1
      for r in range(1,n+1):
        for A in combinations(range(n),r):
          sup,zH=stab_support(NH,A);_,zG=stab_support(NG,A)
          if zH!=zG:de=min(de,len(sup))
      dfs.append(de)
    return dt,min(dfs),max(dfs)

for n in range(4,9):
    print(f'K_{n}: d_targ, min d_fail, max d_fail = {mins(n)}')
