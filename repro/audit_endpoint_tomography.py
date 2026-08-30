"""Exhaustive audit of the endpoint-pair tomography corollary.

Checks every edge of every connected NetworkX graph-atlas graph on 2--7
vertices.  The direct stabilizer-support condition is compared against the
claimed graph criterion: a leaf endpoint or equal external neighborhoods.
"""
import networkx as nx
from raw_repair_context import graph_stab

instances=0
for G in nx.graph_atlas_g():
    n=G.number_of_nodes()
    if not (2 <= n <= 7) or not nx.is_connected(G):
        continue
    E=sorted(tuple(sorted(e)) for e in G.edges())
    target=graph_stab(n,E)
    for e in E:
        u,v=e
        H=graph_stab(n,[f for f in E if f!=e])
        pairmask=(1<<u)|(1<<v)
        direct=False
        for (x,z),sgn in H.items():
            if (x,z) in target or (x,z)==(0,0):
                continue
            support=x|z
            if support & ~pairmask == 0:
                direct=True; break
        Nu=set(G.neighbors(u))-{v}
        Nv=set(G.neighbors(v))-{u}
        predicted=(len(Nu)==0 or len(Nv)==0 or Nu==Nv)
        assert direct==predicted, (n,E,e,direct,predicted,Nu,Nv)
        instances += 1
print(f'PASS endpoint-pair tomography criterion on {instances} edge instances through n=7')
