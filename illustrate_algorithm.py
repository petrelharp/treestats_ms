"""
Draw an illustration of the branch stats algorithm.
"""

import msprime
import tskit
import numpy as np
import subprocess


def convert_svg(svg_file):
    pdf_file = svg_file.replace(".svg", ".pdf")
    print("converting", svg_file, pdf_file)
    subprocess.check_call(f"rsvg-convert {svg_file} -o {pdf_file}", shell=True)

def print_latex_table(pi, beta, x, F):
    print("\\begin{tabular}{c|cccc}")
    print("& $\pi$ & $\\beta$ & $x$ & $F$\\\\")
    for u in range(len(pi)):
        print(f"{u} & {pi[u]} & {beta[u]} & {x[u]} & {F[u]}\\\\")
    print("\\end{tabular}")


def branch_stat_algorithm(ts, f, S, w):

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

    tree = tskit.Tree(ts)
    # TODO Better SVG for not drawing this.
    edge_attrs = {u: {"stroke": "white"} for u in range(ts.num_nodes)}

    for (t_left, t_right), edges_out, edges_in in ts.edge_diffs():
        for edge_num, edge in enumerate(edges_out):
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

            print("out", edge_num)
            print_latex_table(pi, beta, x, F)
            edge_attrs[edge.child] = {"stroke": "white"}
            svg_file = f"figures/tree_{tree.index}_out_{edge_num}.svg"
            tree.draw_svg(
                svg_file,
                max_tree_height="ts",
                edge_attrs=edge_attrs)
            convert_svg(svg_file)

        tree.next()

        for edge_num, edge in enumerate(edges_in):
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

            print("in", edge_num, edge)
            print_latex_table(pi, beta, x, F)

            edge_attrs[edge.child] = {"stroke": "black"}
            svg_file = f"figures/tree_{tree.index}_in_{edge_num}.svg"
            tree.draw_svg(
                svg_file,
                max_tree_height="ts",
                edge_attrs=edge_attrs)
            convert_svg(svg_file)

        sigma += (t_right - t_left) * s

    return sigma / ts.sequence_length

def illustrate(ts):
    n = ts.num_samples
    w = np.ones(ts.num_samples)
    f = lambda x: x * (n - x) / (n * (n - 1))
    branch_stat_algorithm(ts, f, ts.samples(), w)



if __name__ == "__main__":
    # ts = msprime.simulate(8, recombination_rate=0.2, random_seed=42)
    tables = tskit.TableCollection(2)
    for j in range(5):
        tables.nodes.add_row(flags=1, time=0)
    tables.nodes.add_row(time=1)
    tables.nodes.add_row(time=1)
    tables.nodes.add_row(time=2)
    tables.nodes.add_row(time=3)

    tables.edges.add_row(0, 2, 6, 2)
    tables.edges.add_row(0, 2, 6, 1)
    tables.edges.add_row(0, 1, 6, 3)
    tables.edges.add_row(0, 1, 7, 0)
    tables.edges.add_row(0, 2, 7, 6)
    tables.edges.add_row(0, 2, 8, 4)
    tables.edges.add_row(0, 2, 8, 7)

    tables.edges.add_row(1, 2, 5, 0)
    tables.edges.add_row(1, 2, 5, 3)
    tables.edges.add_row(1, 2, 7, 5)

    tables.sort()
    ts = tables.tree_sequence()

    print(ts.draw_text())
    illustrate(ts)
