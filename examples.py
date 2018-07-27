import msprime
import numpy as np
from treestats import *

## Assign angles in [0,2*pi] to each sample, then compute the average angle
## for each internal node of each tree of all samples below each node.

ts = msprime.simulate(10, recombination_rate=1.0)

sample_angles = np.random.uniform(0.0, 2*np.pi, ts.num_samples)

def angle_mean(a):
    # gives the angle of the average vector in the unit circle
    # https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    return np.math.atan2(sum(np.sin(a)), sum(np.cos(a)))

wt = weighted_trees(ts, [sample_angles], angle_mean)

for t in wt:
    print("Tree ", t.index)
    for u, w in zip(t.nodes(), t.node_weights()):
        print("   node ", u, ": angle ", w)
