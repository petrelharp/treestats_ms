import msprime
import itertools

def weighted_trees(ts, sample_weight_list, node_fun=sum):
    '''
    Here ``sample_weight_list`` is a list of lists of weights, each of the same
    length as the samples in the tree sequence ``ts``. This returns an iterator
    over the trees in ``ts`` that is identical to ``ts.trees()`` except that
    each tree ``t`` has the additional method `t.node_weights()` which returns
    an iterator over the "weights" for each node in the tree, in the same order
    as ``t.nodes()``.

    Each node has one weight, computed separately for each set of weights in
    ``sample_weight_list``. Each such weight is defined for a particular list
    of ``sample_weights`` recursively:

    1. First define ``all_weights[ts.samples()[j]] = sample_weights[j]``
        and ``all_weights[k] = 0`` otherwise.
    2. The weight for a node ``j`` with children ``u1, u2, ..., un`` is
        ``node_fun([all_weights[j], weight[u1], ..., weight[un]])``.

    For instance, if ``sample_weights`` is a vector of all ``1``s, and
    ``node_fun`` is ``sum``, then the weight for each node in each tree
    is the number of samples below it, equivalent to ``t.num_samples(j)``.

    To do this, we need to only recurse upwards from the parent of each
    added or removed edge, updating the weights.
    '''
    samples = ts.samples()
    num_weights = len(sample_weight_list)
    # initialize the weights
    base_X = [[0.0 for _ in range(num_weights)] for _ in range(ts.num_nodes)]
    X = [[0.0 for _ in range(num_weights)] for _ in range(ts.num_nodes)]
    for j, u in enumerate(samples):
        for k in range(num_weights):
            X[u][k] = sample_weight_list[k][j]
            base_X[u][k] = sample_weight_list[k][j]

    for t, (interval, records_out, records_in) in zip(ts.trees(), ts.edge_diffs()):
        for edge in itertools.chain(records_out, records_in):
            u = edge.parent
            while u != -1:
                for k in range(num_weights):
                    U = [base_X[u][k]] + [X[u][k] for u in t.children(u)]
                    X[u][k] = node_fun(U)
                u = t.parent(u)

        def the_node_weights(self):
            for u in self.nodes():
                yield X[u]

        # magic that uses "descriptor protocol"
        t.node_weights = the_node_weights.__get__(t, msprime.trees.SparseTree)
        yield t



