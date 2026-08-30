"""Independent verifier for the diagnostic context number c_diag(G)."""
import itertools

def f2_kernel_basis(rows,n):
    ech=[];pivs=[]
    for r in rows:
        cur=r
        for p,er in zip(pivs,ech):
            if (cur>>p)&1:cur^=er
        if cur:
            p=(cur&-cur).bit_length()-1
            for i,er in enumerate(ech):
                if (er>>p)&1:ech[i]^=cur
            pivs.append(p);ech.append(cur)
    free=[c for c in range(n) if c not in set(pivs)];basis=[]
    for fc in free:
        vec=1<<fc
        for p,er in zip(pivs,ech):
            if (er>>fc)&1:vec|=1<<p
        basis.append(vec)
    return basis

def codewords_from_basis(basis):
    words=[0]
    for b in basis:words += [w^b for w in words]
    return words

class CdiagSolver:
    def __init__(self,n,edges):
        self.n=n;self.edges=[tuple(sorted(e)) for e in edges];self.m=len(self.edges);self.Gamma=[0]*n
        for u,v in self.edges:self.Gamma[u]|=1<<v;self.Gamma[v]|=1<<u
        self.pairs=list(itertools.combinations(range(self.m+1),2));self.npairs=len(self.pairs);self.edge_masks=[(1<<u)|(1<<v) for u,v in self.edges]
    def context_code(self,beta):
        rows=[]
        for v in range(self.n):
            if beta[v]==0:rows.append(self.Gamma[v])
            elif beta[v]==1:rows.append(self.Gamma[v]^(1<<v))
            else:rows.append(1<<v)
        return codewords_from_basis(f2_kernel_basis(rows,self.n))
    def sep_mask(self,beta):
        mask=0
        for A in self.context_code(beta):
            if A==0:continue
            hits=0
            for i,em in enumerate(self.edge_masks):
                if A&em:hits|=1<<i
            for idx,(i,j) in enumerate(self.pairs):
                hi=(hits>>i)&1 if i<self.m else 0;hj=(hits>>j)&1 if j<self.m else 0
                if hi!=hj:mask|=1<<idx
        return mask

def wheel(n):
    m=n-1;E=[(0,i) for i in range(1,n)]+[(i,i%m+1) for i in range(1,n)]
    return n,sorted(set(tuple(sorted(e)) for e in E))
