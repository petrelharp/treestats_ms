"""
Direct implementations of the algorithms in the paper.
"""
import msprime
import tskit
import numpy as np

def diffs_branch_algorithm(ts, f, S, w):

    tau = ts.tables.nodes.time
    pi = np.zeros(ts.num_nodes, dtype=np.int32) - 1
    beta = np.zeros(ts.num_nodes)
    x = np.zeros(ts.num_nodes)
    F = np.zeros(ts.num_nodes)
    sigma = 0
    s = 0

    for j, u in enumerate(S):
        x[u] = w[j]
        F[u] = f(x[u])

    for (t_left, t_right), edges_out, edges_in in ts.edge_diffs():
        for edge in edges_out:
            u, v = edge.child, edge.parent
            s -= beta[u] * F[u]
            pi[u] = -1
            beta[u] = 0

            while v != -1:
                s -= beta[v] * F[v]
                x[v] -= x[u]
                F[v] = f(x[v])
                s += beta[v] * F[v]
                v = pi[v]

        for edge in edges_in:
            u, v = edge.child, edge.parent
            pi[u] = v
            beta[u] = tau[v] - tau[u]
            s += beta[u] * F[u]

            while v != -1:
                s -= beta[v] * F[v]
                x[v] += x[u]
                F[v] = f(x[v])
                s += beta[v] * F[v]
                v = pi[v]

        sigma += (t_right - t_left) * s

    return sigma / ts.sequence_length


def tables_branch_algorithm(ts, f, S, w):

    N = ts.tables.nodes
    E = ts.tables.edges
    L = ts.sequence_length

    tau = N.time
    pi = np.zeros(len(N), dtype=np.int32) - 1
    beta = np.zeros(len(N))
    x = np.zeros(len(N))
    F = np.zeros(len(N))
    sigma = 0
    s = 0

    for j, u in enumerate(S):
        x[u] = w[j]
        F[u] = f(x[u])

    I = sorted(
        range(len(E)), key=lambda j: (
            E[j].left, tau[E[j].parent], E[j].parent, E[j].child))
    O = sorted(
        range(len(E)), key=lambda j: (
            E[j].right, -tau[E[j].parent], -E[j].parent, -E[j].child))
    j = 0
    k = 0
    t_l = 0

    while j < len(E) or t_l < L:
        while k < len(E) and E[O[k]].right == t_l:
            u = E[O[k]].child
            v = E[O[k]].parent
            k += 1

            s -= beta[u] * F[u]
            pi[u] = -1
            beta[u] = 0

            while v != -1:
                s -= beta[v] * F[v]
                x[v] -= x[u]
                F[v] = f(x[v])
                s += beta[v] * F[v]
                v = pi[v]

        while j < len(E) and E[I[j]].left == t_l:
            u = E[I[j]].child
            v = E[I[j]].parent
            j += 1

            pi[u] = v
            beta[u] = tau[v] - tau[u]
            s += beta[u] * F[u]

            while v != -1:
                s -= beta[v] * F[v]
                x[v] += x[u]
                F[v] = f(x[v])
                s += beta[v] * F[v]
                v = pi[v]

        t_r = L
        if j < len(E):
            t_r = min(t_r, E[I[j]].left)
        if k < len(E):
            t_r = min(t_r, E[O[k]].right)

        sigma += (t_r - t_l) * s
        t_l = t_r

    return sigma / L


def check_algorithms(ts, f):
    # print(ts.draw_text())
    weights = np.ones((ts.num_samples, 1))

    sigma1 = ts.general_stat(
        weights, f, 1, strict=False, polarised=True, mode="branch")[0]

    w = weights.reshape(ts.num_samples)
    sigma2 = diffs_branch_algorithm(ts, f, ts.samples(), w)
    sigma3 = tables_branch_algorithm(ts, f, ts.samples(), w)
    print(sigma1, sigma2, sigma3)

if __name__ == "__main__":
    ts = msprime.simulate(30, recombination_rate=10.0, random_seed=2)
    check_algorithms(ts, lambda x: x + 1)



