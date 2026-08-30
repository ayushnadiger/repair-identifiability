import itertools, math, networkx as nx, collections
from functools import lru_cache

def graph_stab(n,E):
    nbr=[0]*n
    for u,v in E:nbr[u]|=1<<v;nbr[v]|=1<<u
    elems={(0,0):0}
    for v in range(n):
        gx=1<<v;gz=nbr[v];new={}
        for (x,z),p in list(elems.items()):
            q=(p+2*((z&gx).bit_count()&1))%4
            new[(x^gx,z^gz)]=q
        elems.update(new)
    out={}
    for (x,z),p in elems.items():
        k=(x&z).bit_count();r=(p-k)%4;assert r in(0,2)
        out[(x,z)]=1 if r==0 else -1
    return out

def wheel(n):
 m=n-1;E=[(0,i) for i in range(1,n)]+[(i,(i%m)+1) for i in range(1,n)];return sorted(set(tuple(sorted(e)) for e in E))

def context_paulis(n,b):
  vals=[]
  for q in range(1<<n):
    x=z=0
    for v in range(n):
      if (q>>v)&1:
        if b[v]==0:x|=1<<v
        elif b[v]==1:x|=1<<v;z|=1<<v
        else:z|=1<<v
    vals.append((x,z))
  return vals

def signatures_context(n,E,b,include_null=False):
  E=sorted(set(tuple(sorted(e)) for e in E));target=graph_stab(n,E)
  cuts=[graph_stab(n,[f for f in E if f!=e]) for e in E]
  Ps=context_paulis(n,b)
  H=[]; sig=[]
  if include_null:
    H.append(('N',None));sig.append(tuple(target.get(p,0) for p in Ps))
  for j,e in enumerate(E):
    H.append(('C',j));sig.append(tuple(cuts[j].get(p,0) for p in Ps))
    H.append(('D',j));
    sig.append(tuple((0 if any((x>>v)&1 for v in e) else target.get((x,z),0)) for x,z in Ps))
  return H,sig

def sep_mask(n,E,b,include_null=False):
  H,sig=signatures_context(n,E,b,include_null)
  M=0
  for k,(i,j) in enumerate(itertools.combinations(range(len(H)),2)):
    if sig[i]!=sig[j]:M|=1<<k
  return M,len(H)

def minrawrepair(n,E,include_null=False,maxk=5):
  d={};Hn=None
  for b in itertools.product(range(3),repeat=n):
    M,Hn=sep_mask(n,E,b,include_null)
    if M:d.setdefault(M,b)
  full=(1<<(Hn*(Hn-1)//2))-1
  ml=sorted(d,key=int.bit_count,reverse=True);maximal=[]
  for M in ml:
    if not any((M|Q)==Q for Q,_ in maximal):maximal.append((M,d[M]))
  for M,b in maximal:
    if M==full:return 1,[b],len(d),len(maximal)
  np=Hn*(Hn-1)//2;cand=[[] for _ in range(np)]
  for idx,(M,b) in enumerate(maximal):
    x=M
    while x:
      lb=x&-x;r=lb.bit_length()-1;cand[r].append(idx);x^=lb
  @lru_cache(None)
  def dfs(cov,dep):
    if cov==full:return ()
    if dep==0:return None
    unc=full^cov
    best=None;x=unc
    while x:
      lb=x&-x;r=lb.bit_length()-1
      us=[i for i in cand[r] if maximal[i][0]&unc]
      if best is None or len(us)<len(best[1]):best=(r,us)
      x^=lb
    if best is None or not best[1]:return None
    mg=max((maximal[i][0]&unc).bit_count() for i in best[1])
    if (unc.bit_count()+mg-1)//mg>dep:return None
    seen=set()
    for i in sorted(best[1],key=lambda i:(maximal[i][0]&unc).bit_count(),reverse=True):
      nc=cov|maximal[i][0]
      if nc in seen:continue
      seen.add(nc)
      r=dfs(nc,dep-1)
      if r is not None:return (i,)+r
    return None
  for k in range(2,maxk+1):
    dfs.cache_clear();a=dfs(0,k)
    if a is not None:return k,[maximal[i][1] for i in a],len(d),len(maximal)
  return f'>{maxk}',None,len(d),len(maximal)

def fmt(b):return ''.join('XYZ'[x] for x in b)

if __name__=='__main__':
 for n in range(3,8):
  fam=[(f'P{n}',[(i,i+1) for i in range(n-1)]),(f'S{n}',[(0,i) for i in range(1,n)]),(f'W{n}',wheel(n))]
  if n<=6:fam.append((f'K{n}',[(i,j) for i in range(n) for j in range(i+1,n)]))
  for name,E in fam:
    if len(E)>12:continue
    k,bs,nd,nm=minrawrepair(n,E,maxk=4)
    print(name,'m',len(E),'c_raw_repair',k,[fmt(b) for b in bs] if bs else None,'masks',nd,'max',nm)
