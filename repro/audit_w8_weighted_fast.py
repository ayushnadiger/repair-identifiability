import itertools,sys
from raw_repair_context import graph_stab,wheel,fmt
n=8; E=wheel(n); E=sorted(E); T=graph_stab(n,E); cuts=[graph_stab(n,[f for f in E if f!=e]) for e in E]

def prep(mode):
    H=[]
    if mode=='cut':
        H=[('S',T)]+[('S',S) for S in cuts]
    else:
        H=[('S',T)]
        for j,e in enumerate(E):H += [('S',cuts[j]),('D',(T,e))]
    pairs=list(itertools.combinations(range(len(H)),2)); full=(1<<len(pairs))-1
    sep=[0]*(1<<(2*n)); wt=[0]*(1<<(2*n))
    for x in range(1<<n):
      for z in range(1<<n):
        code=x|(z<<n); wt[code]=(x|z).bit_count()
        vals=[]
        for typ,obj in H:
          if typ=='S':vals.append(obj.get((x,z),0))
          else:
            target,e=obj
            vals.append(0 if ((x>>e[0])&1 or (x>>e[1])&1) else target.get((x,z),0))
        m=0
        for k,(i,j) in enumerate(pairs):
          if vals[i]!=vals[j]:m|=1<<k
        sep[code]=m
    return sep,wt,full,len(H)

def context_codes(b):
    # all subset products; local X=0,Y=1,Z=2
    gens=[]
    for v,t in enumerate(b):
      if t==0: x,z=1<<v,0
      elif t==1:x,z=1<<v,1<<v
      else:x,z=0,1<<v
      gens.append(x|(z<<n))
    codes=[0]
    for g in gens:codes += [q^g for q in list(codes)]
    return codes

def masks_all_contexts(mode):
  sep,wt,full,H=prep(mode); out=[{} for _ in range(n+1)]
  for b in itertools.product(range(3),repeat=n):
    acc=[0]*(n+1)
    bucket=[0]*(n+1)
    for code in context_codes(b)[1:]:
      w=wt[code]
      bucket[w]|=sep[code]
    cur=0
    for w in range(1,n+1):
      cur |= bucket[w]
      if cur: out[w].setdefault(cur,b)
  return out,full

def prune(d):
  items=sorted(d.items(),key=lambda kv:kv[0].bit_count(),reverse=True);maxi=[]
  for M,b in items:
    if not any((M|Q)==Q for Q,_ in maxi):maxi.append((M,b))
  return maxi

def mincover(d,full,maxk=3):
  maxi=prune(d)
  for M,b in maxi:
    if M==full:return 1,[b],len(d),len(maxi)
  # exact two cover
  for i,(M,b) in enumerate(maxi):
    need=full^M
    for j in range(i,len(maxi)):
      Q,bq=maxi[j]
      if need & ~Q ==0:return 2,[b,bq],len(d),len(maxi)
  if maxk<3:return '>2',None,len(d),len(maxi)
  # exact three via uncovered-pair candidate lists recursive
  np=full.bit_length(); cand=[[] for _ in range(np)]
  for idx,(M,b) in enumerate(maxi):
    x=M
    while x:
      lb=x&-x;r=lb.bit_length()-1;cand[r].append(idx);x^=lb
  from functools import lru_cache
  @lru_cache(None)
  def dfs(cov,dep):
    if cov==full:return ()
    if dep==0:return None
    unc=full^cov; best=None; x=unc
    while x:
      lb=x&-x;r=lb.bit_length()-1
      us=[i for i in cand[r] if maxi[i][0]&unc]
      if best is None or len(us)<len(best):best=us
      x^=lb
    if not best:return None
    for i in sorted(best,key=lambda i:(maxi[i][0]&unc).bit_count(),reverse=True):
      r=dfs(cov|maxi[i][0],dep-1)
      if r is not None:return (i,)+r
    return None
  ans=dfs(0,3)
  if ans is None:return '>3',None,len(d),len(maxi)
  return 3,[maxi[i][1] for i in ans],len(d),len(maxi)

for mode in ['cut','repair']:
  print('building',mode,flush=True); ds,full=masks_all_contexts(mode)
  for w in range(1,n+1):
    k,cert,nd,nm=mincover(ds[w],full,3)
    print(mode,'w',w,'k',k,'cert',[fmt(b) for b in cert] if cert else None,'distinct',nd,'max',nm,flush=True)
