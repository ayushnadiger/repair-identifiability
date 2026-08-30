"""Small exact audit of the single-edge blindness/cure theorem."""
import itertools, random
import numpy as np

rng=random.Random(20260817)

def graph_state(n, edges):
    N=1<<n
    psi=np.ones(N,dtype=complex)/np.sqrt(N)
    for x in range(N):
        phase=1
        for u,v in edges:
            if ((x>>u)&1) and ((x>>v)&1): phase*=-1
        psi[x]*=phase
    return psi

def pauli_expect(psi, A, edges):
    # target stabilizer K_A: X on A, Z on odd neighborhood Gamma A; include exact stabilizer phase by applying generators sequentially
    n=int(np.log2(len(psi)))
    # apply product of K_a directly to basis
    out=psi.copy()
    # matrix-free sequential K_a applications
    for a in A:
        new=np.zeros_like(out)
        neigh=[v if u==a else u for u,v in edges if u==a or v==a]
        for x,amp in enumerate(out):
            zphase=(-1)**sum((x>>w)&1 for w in neigh)
            new[x^(1<<a)] += zphase*amp
        out=new
    return float(np.real_if_close(np.vdot(psi,out)))

def dephased_target_expect(A,u,v):
    # target stabilizer expectation on ideal is 1; independent Z dephasing kills iff A contains endpoint
    return 1.0 if u not in A and v not in A else 0.0

def boundary_expect(psi, u, neigh):
    n=int(np.log2(len(psi)))
    out=np.zeros_like(psi)
    for x,amp in enumerate(psi):
        phase=(-1)**sum((x>>w)&1 for w in neigh)
        out[x^(1<<u)] += phase*amp
    return float(np.real_if_close(np.vdot(psi,out)))

cases=0
for n in range(2,7):
    verts=list(range(n))
    allE=list(itertools.combinations(verts,2))
    for _ in range(12):
        E=[e for e in allE if rng.random()<0.35]
        if not E: continue
        target=graph_state(n,E)
        for e in E[:min(3,len(E))]:
            u,v=e; Ep=[f for f in E if f!=e]
            cut=graph_state(n,Ep)
            for mask in range(1<<n):
                A=[i for i in range(n) if (mask>>i)&1]
                x=pauli_expect(cut,A,E)
                y=dephased_target_expect(A,u,v)
                if abs(x-y)>1e-8:
                    raise AssertionError((n,E,e,A,x,y))
            Nu=[w for w in verts if tuple(sorted((u,w))) in Ep]
            Nv=[w for w in verts if tuple(sorted((v,w))) in Ep]
            if abs(boundary_expect(cut,u,Nu)-1)>1e-8 or abs(boundary_expect(cut,v,Nv)-1)>1e-8:
                raise AssertionError('boundary cure failed on cut state')
            cases+=1
print(f'PASS single-edge blindness/cure on {cases} edge instances through n=6')
