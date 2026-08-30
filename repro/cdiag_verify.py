"""Independent verifier for the diagnostic context number c_diag(G).

Model (matching the paper):
  - Context beta: V -> {X,Y,Z}. Support code C_beta = ker M_beta over F2, where
    row at v is Gamma_v (X), Gamma_v + e_v (Y), e_v (Z).
  - Measurable stabilizers in beta are exactly s_A with 1_A in C_beta.
  - Signature of s_A under cut at e: sigma_e = 1[A cap e = empty]; under null: 1.
  - Hypotheses: all edges + null. A pair is separated by beta if some codeword's
    hit indicators differ on the pair (for null: some codeword hits the edge).
  - cdiag = min #contexts whose separated-pair sets union to all pairs.
"""
import itertools
import numpy as np


def f2_kernel_basis(rows, n):
    rows = [r for r in rows]
    pivots = {}
    ech = []
    for r in rows:
        cur = r
        for c, er in pivots.items():
            if (cur >> c) & 1:
                cur ^= ech[er]
        if cur:
            c = (cur & -cur).bit_length() - 1
            pivots[c] = len(ech)
            ech.append(cur)
    free_cols = [c for c in range(n) if c not in set(pivots.keys())]
    basis = []
    for fc in free_cols:
        basis.append(1 << fc)
    if not free_cols:
        return []
    ech2 = []
    piv2 = []
    for r in rows:
        cur = r
        for p, er in zip(piv2, ech2):
            if (cur >> p) & 1:
                cur ^= er
        if cur:
            p = (cur & -cur).bit_length() - 1
            ech2_new = []
            for q, er in zip(piv2, ech2):
                if (er >> p) & 1:
                    er ^= cur
                ech2_new.append(er)
            ech2 = ech2_new
            piv2.append(p)
            ech2.append(cur)
    pivset = set(piv2)
    free = [c for c in range(n) if c not in pivset]
    basis = []
    for fc in free:
        vec = 1 << fc
        for p, er in zip(piv2, ech2):
            if (er >> fc) & 1:
                vec |= 1 << p
        basis.append(vec)
    return basis


def codewords_from_basis(basis, cap=1 << 14):
    d = len(basis)
    if (1 << d) > cap:
        raise RuntimeError(f"kernel too large: dim {d}")
    words = [0]
    for b in basis:
        words += [w ^ b for w in words]
    return words


class CdiagSolver:
    def __init__(self, n, edges):
        self.n = n
        self.edges = [tuple(sorted(e)) for e in edges]
        self.m = len(self.edges)
        self.Gamma = [0] * n
        for (u, v) in self.edges:
            self.Gamma[u] |= 1 << v
            self.Gamma[v] |= 1 << u
        self.H = self.m + 1
        self.pairs = list(itertools.combinations(range(self.H), 2))
        self.npairs = len(self.pairs)
        self.pair_index = {p: i for i, p in enumerate(self.pairs)}
        self.edge_masks = [(1 << u) | (1 << v) for (u, v) in self.edges]

    def context_code(self, beta):
        rows = []
        for v in range(self.n):
            if beta[v] == 0:
                rows.append(self.Gamma[v])
            elif beta[v] == 1:
                rows.append(self.Gamma[v] | (1 << v) if not ((self.Gamma[v] >> v) & 1) else self.Gamma[v] ^ (1 << v))
            else:
                rows.append(1 << v)
        basis = f2_kernel_basis(rows, self.n)
        return codewords_from_basis(basis)

    def sep_mask(self, beta):
        words = self.context_code(beta)
        mask = 0
        for A in words:
            if A == 0:
                continue
            hits = 0
            for i, em in enumerate(self.edge_masks):
                if A & em:
                    hits |= 1 << i
            if hits == 0:
                continue
            for idx, (i, j) in enumerate(self.pairs):
                hi = (hits >> i) & 1 if i < self.m else 0
                hj = (hits >> j) & 1 if j < self.m else 0
                if hi != hj:
                    mask |= 1 << idx
        return mask

    def all_masks(self, verbose=False):
        masks = set()
        for beta in itertools.product((0, 1, 2), repeat=self.n):
            masks.add(self.sep_mask(beta))
        return masks

    def cdiag(self, max_k=4):
        full = (1 << self.npairs) - 1
        masks = self.all_masks()
        masks.discard(0)
        ml = sorted(masks, key=lambda x: -bin(x).count('1'))
        maximal = []
        for m_ in ml:
            if not any((m_ | M) == M for M in maximal):
                maximal.append(m_)
        if any(m_ == full for m_ in maximal):
            return 1
        for a in range(len(maximal)):
            for b in range(a + 1, len(maximal)):
                if (maximal[a] | maximal[b]) == full:
                    return 2
        for a in range(len(maximal)):
            for b in range(a + 1, len(maximal)):
                ab = maximal[a] | maximal[b]
                for c in range(b + 1, len(maximal)):
                    if (ab | maximal[c]) == full:
                        return 3
        union = 0
        for m_ in maximal:
            union |= m_
        if union != full:
            return None
        return ">=4"


def wheel(n):
    m = n - 1
    edges = [(0, i) for i in range(1, n)]
    edges += [(i, i % m + 1) for i in range(1, n)]
    return n, sorted(set(tuple(sorted(e)) for e in edges))


def complete_bipartite(a, b):
    edges = [(i, a + j) for i in range(a) for j in range(b)]
    return a + b, edges


if __name__ == "__main__":
    tests = [
        ("K_3", 3, [(0, 1), (0, 2), (1, 2)], 2),
        ("P_4", 4, [(0, 1), (1, 2), (2, 3)], 2),
        ("K_4", 4, [(i, j) for i in range(4) for j in range(i + 1, 4)], 1),
        ("K_5", 5, [(i, j) for i in range(5) for j in range(i + 1, 5)], 1),
        ("star K_{1,4}", 5, [(0, i) for i in range(1, 5)], 1),
        ("C_5", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)], 2),
    ]
    for name, n, E, expect in tests:
        s = CdiagSolver(n, E)
        got = s.cdiag()
        print(f"{name}: cdiag = {got} (expected {expect}) {'OK' if got == expect else '**MISMATCH**'}")
