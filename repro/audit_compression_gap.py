"""Exact raw-record compression-gap checks for K3 and W6.

The target-stabilizer invariant can exceed unrestricted raw product-context
complexity because compatible non-target Pauli correlators are present in the
same local measurement record.
"""
from raw_repair_context import graph_stab, wheel


def compatible_signature(n, stab, beta):
    vals=[]
    for q in range(1<<n):
        x=z=0
        for v in range(n):
            if (q>>v)&1:
                if beta[v]==0: x |= 1<<v
                elif beta[v]==1: x |= 1<<v; z |= 1<<v
                else: z |= 1<<v
        vals.append(stab.get((x,z),0))
    return tuple(vals)


def cut_catalogue_distinct(n,E,beta):
    E=sorted(set(tuple(sorted(e)) for e in E))
    states=[graph_stab(n,E)]
    states += [graph_stab(n,[f for f in E if f!=e]) for e in E]
    sig=[compatible_signature(n,S,beta) for S in states]
    return len(sig)==len(set(sig)),sig

# K3: XXX.
E3=[(0,1),(0,2),(1,2)]
ok,_=cut_catalogue_distinct(3,E3,(0,0,0))
assert ok
print('K3: XXX separates intact + all three cuts using full raw correlations.')

# W6: XYYYYY (hub first).
E6=wheel(6)
ok,_=cut_catalogue_distinct(6,E6,(0,1,1,1,1,1))
assert ok
print('W6: XYYYYY separates intact + all ten cuts using full raw correlations.')
