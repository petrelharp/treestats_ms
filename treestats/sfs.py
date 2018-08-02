# Experimental code for computing functions of the joint SFS along a sequence
# and related things

import msprime
import numpy as np

def get_derived_state(site, mut_id):
    """
    Find the derived state of the mutation with id `mut_id` at site `site`.
    """
    if mut_id == -1:
        state = site.ancestral_state
    else:
        for m in site.mutations:
            if m.id == mut_id:
                state = m.derived_state
    return state

def tracked_samples_joint_branch_sfs(tree_sequence, pop0, pop1, begin=0.0, end=None):
    if end is None:
        end = tree_sequence.sequence_length
    n0 = len(pop0)
    n1 = len(pop1)
    S = np.zeros((n0 + 1, n1 + 1))
    tt0 = tree_sequence.trees(tracked_samples=pop0,
                                   sample_counts=True)
    tt1 = tree_sequence.trees(tracked_samples=pop1,
                                   sample_counts=True)
    for t0, t1 in zip(tt0, tt1):
        root = t0.root
        tr_len = min(end, t0.interval[1]) - max(begin, t0.interval[0])
        if tr_len > 0:
            for node in t0.nodes():
                if node != root:
                    x0 = t0.num_tracked_samples(node)
                    x1 = t1.num_tracked_samples(node)
                    if (x0 > 0) or (x1 > 0):
                        S[x0, x1] += t0.branch_length(node) * tr_len
    S /= (end-begin)
    return S


def variants_joint_sfs(tree_sequence, pop0, pop1, begin=0.0, end=None):
    '''
    '''
    if end is None:
        end = tree_sequence.sequence_length
    n0 = len(pop0)
    n1 = len(pop1)
    S = np.zeros((n0 + 1, n1 + 1))
    for v in tree_sequence.variants():
        if v.position >= end:
            break
        if v.position >= begin:
            for a in set(v.genotypes):
                x0 = sum(v.genotypes[pop0] == a)
                x1 = sum(v.genotypes[pop1] == a)
                if (x0 > 0) or (x1 > 0):
                    S[x0, x1] += 1.0
    S /= (end - begin)
    return S


def joint_site_frequency_spectrum(tree_sequence, pop0, pop1, windows=None):
    '''
    Computes the joint folded site frequency spectrum between pop0 and pop1,
    independently in windows.

    :param list pop0: A list of IDs of samples of length n0.
    :param list pop1: A list of IDs of samples of length n1.
    :param iterable windows: The breakpoints of the windows (including start
        and end, so has one more entry than number of windows).
    :return: A list of matrices of dimension n0 x n1, one for each window,
        whose kth entry gives the number of mutations in that window at which a
        mutation is seen by exactly k of the samples, divided by the window length.
    '''
    if windows is None:
        windows = (0, tree_sequence.sequence_length)
    if ((not isinstance(pop0, list)) or
         (not isinstance(pop1, list)) or
         (len(pop0) + len(pop1) != len(set(pop0+pop1)))):
        raise ValueError(
            "pop0 and pop1 must not contain repeated elements.")
    if (len(pop0) == 0) or (len(pop1) == 0):
        raise ValueError("sample_set cannot be empty.")
    for u in pop0 + pop1:
        if not tree_sequence.node(u).is_sample():
            raise ValueError("Not all elements of pop0 and pop1 are samples.")
    num_windows = len(windows) - 1
    if windows[0] != 0.0:
        raise ValueError(
            "Windows must start at the start of the sequence (at 0.0).")
    if windows[-1] != tree_sequence.sequence_length:
        raise ValueError("Windows must extend to the end of the sequence.")
    for k in range(num_windows):
        if windows[k + 1] <= windows[k]:
            raise ValueError("Windows must be increasing.")
    num_sites = tree_sequence.num_sites
    n0 = len(pop0)
    n1 = len(pop1)
    # we store the final answers here
    S = np.zeros((num_windows, n0 + 1, n1 + 1))
    if num_sites == 0:
        return S
    N = tree_sequence.num_nodes
    # initialize: with no tree, each node is either in a sample set or not
    X0 = [int(u in pop0) for u in range(N)]
    X1 = [int(u in pop1) for u in range(N)]
    # we will construct the tree here
    pi = [-1 for j in range(N)]
    # keep track of which site we're looking at
    sites = tree_sequence.sites()
    ns = 0  # this will record number of sites seen so far
    s = next(sites)
    # index of *left-hand* end of the current window
    window_num = 0
    while s.position > windows[window_num + 1]:
        window_num += 1
    for interval, records_out, records_in in tree_sequence.edge_diffs():
        # if we've done all the sites then stop
        if ns == num_sites:
            break
        # update the tree
        for sign, records in ((-1, records_out), (+1, records_in)):
            for edge in records:
                dx1 = 0
                dx2 = 0
                if sign == +1:
                    pi[edge.child] = edge.parent
                dx1+= sign * X0[edge.child]
                dx2+= sign * X1[edge.child]
                if sign == -1:
                    pi[edge.child] = -1
                X0[edge.parent] += dx1
                X1[edge.parent] += dx2
                # propagate change up the tree
                u = pi[edge.parent]
                if u != -1:
                    next_u = pi[u]
                    while u != -1:
                        X0[u] += dx1
                        X1[u] += dx2
                        u = next_u
                        next_u = pi[next_u]
        # loop over sites in this tree
        while s.position < interval[1]:
            if s.position > windows[window_num + 1]:
                # finalize this window and move to the next
                window_length = windows[window_num + 1] - windows[window_num]
                S[window_num] /= window_length
                # may need to advance through empty windows
                while s.position > windows[window_num + 1]:
                    window_num += 1
            nm = len(s.mutations)
            if nm > 0:
                U = {s.ancestral_state: [n0, n1]}
                for mut in s.mutations:
                    if mut.derived_state not in U:
                        U[mut.derived_state] = [0, 0]
                    U[mut.derived_state][0] += X0[mut.node]
                    U[mut.derived_state][1] += X1[mut.node]
                    parent_state = get_derived_state(s, mut.parent)
                    if parent_state not in U:
                        U[parent_state] = [0, 0]
                    U[parent_state][0] -= X0[mut.node]
                    U[parent_state][1] -= X1[mut.node]
                for a in U:
                    if max(U[a]) > 0:
                        S[window_num][U[a][0], U[a][1]] += 1.0
            ns += 1
            if ns == num_sites:
                break
            s = next(sites)
    # wrap up the final window
    window_length = windows[window_num + 1] - windows[window_num]
    S[window_num] /= window_length
    return S

def joint_branch_frequency_spectrum(tree_sequence, pop0, pop1, windows=None):
    '''
    Computes the expected *derived* (unfolded) joint site frequency spectrum,
    between pop0 and pop1, based on tree lengths, separately in each window.

    :param list pop0: A list of IDs of samples of length n0.
    :param list pop1: A list of IDs of samples of length n1.
    :param iterable windows: The breakpoints of the windows (including start
        and end, so has one more entry than number of windows).
    :return: A list of lists of length n, one for each window, whose kth
        entry gives the total length of any branches in the marginal trees
        over that window that are ancestral to exactly k of the samples,
        divided by the length of the window.
    '''
    if windows is None:
        windows = (0, tree_sequence.sequence_length)
    if ((not isinstance(pop0, list)) or
         (not isinstance(pop1, list)) or
         (len(pop0) + len(pop1) != len(set(pop0+pop1)))):
        raise ValueError(
            "pop0 and pop1 must not contain repeated elements.")
    if (len(pop0) == 0) or (len(pop1) == 0):
        raise ValueError("sample_set cannot be empty.")
    for u in pop0 + pop1:
        if not tree_sequence.node(u).is_sample():
            raise ValueError("Not all elements of pop0 and pop1 are samples.")
    num_windows = len(windows) - 1
    if windows[0] != 0.0:
        raise ValueError(
            "Windows must start at the start of the sequence (at 0.0).")
    if windows[-1] != tree_sequence.sequence_length:
        raise ValueError("Windows must extend to the end of the sequence.")
    for k in range(num_windows):
        if windows[k + 1] <= windows[k]:
            raise ValueError("Windows must be increasing.")
    n0 = len(pop0)
    n1 = len(pop1)
    S = np.zeros((num_windows, n0 + 1, n1 + 1))
    L = np.zeros((n0 + 1, n1 + 1))
    N = tree_sequence.num_nodes
    X0 = [int(u in pop0) for u in range(N)]
    X1 = [int(u in pop1) for u in range(N)]
    # we will essentially construct the tree
    pi = [-1 for j in range(N)]
    node_time = [tree_sequence.node(u).time for u in range(N)]
    # keep track of where we are for the windows
    chrom_pos = 0.0
    # index of *left-hand* end of the current window
    window_num = 0
    for interval, records_out, records_in in tree_sequence.edge_diffs():
        length = interval[1] - interval[0]
        for sign, records in ((-1, records_out), (+1, records_in)):
            for edge in records:
                dx0 = 0
                dx1 = 0
                if sign == +1:
                    pi[edge.child] = edge.parent
                dx0 += sign * X0[edge.child]
                dx1 += sign * X1[edge.child]
                dt = (node_time[pi[edge.child]] - node_time[edge.child])
                if (X0[edge.child] > 0) or (X1[edge.child] > 0):
                    L[X0[edge.child], X1[edge.child]] += sign * dt
                if sign == -1:
                    pi[edge.child] = -1
                old_X0 = X0[edge.parent]
                old_X1 = X1[edge.parent]
                X0[edge.parent] += dx0
                X1[edge.parent] += dx1
                if pi[edge.parent] != -1:
                    dt = (node_time[pi[edge.parent]] - node_time[edge.parent])
                    if (X0[edge.parent] > 0) or (X1[edge.parent] > 0):
                        L[X0[edge.parent], X1[edge.parent]] += dt
                    if (old_X0 > 0) or (old_X1 > 0):
                        L[old_X0, old_X1] -= dt
                # propagate change up the tree
                u = pi[edge.parent]
                if u != -1:
                    next_u = pi[u]
                    while u != -1:
                        old_X0 = X0[u]
                        old_X1 = X1[u]
                        X0[u] += dx0
                        X1[u] += dx1
                        # need to update X for the root,
                        # but the root does not have a branch length
                        if next_u != -1:
                            dt = (node_time[pi[u]] - node_time[u])
                            if (X0[u] > 0) or (X1[u] > 0):
                                L[X0[u], X1[u]] += dt
                            if (old_X0 > 0) or (old_X1 > 0):
                                L[old_X0, old_X1] -= dt
                        u = next_u
                        next_u = pi[next_u]
        while chrom_pos + length >= windows[window_num + 1]:
            # wrap up the last window
            this_length = windows[window_num + 1] - chrom_pos
            window_length = windows[window_num + 1] - windows[window_num]
            S[window_num] += L * this_length
            S[window_num] /= window_length
            length -= this_length
            # start the next
            if window_num < num_windows - 1:
                window_num += 1
                chrom_pos = windows[window_num]
            else:
                # skips the else statement below
                break
        else:
            S[window_num] += L * length
            chrom_pos += length
    return S

